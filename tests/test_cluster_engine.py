"""Synthetic numerical tests and mocked external adapters (not biological validation)."""
import argparse
import csv
import json
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from cluster_engine import data, evidence, scoring, methods, cli


class ClusterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.genes = ['a', 'b', 'c', 'd']
        self.seqs = {g: 'M' * 20 for g in self.genes}
        self.fa = self.root/'target.fa'
        data.write_fasta(self.fa, self.seqs)
        self.map = self.root/'map.tsv'
        data.tsv(self.map, [dict(gene_ID=g, assembly_ID='asm'+g) for g in self.genes], ['gene_ID', 'assembly_ID'])
        self.good = data.partition([['a', 'b'], ['c', 'd']], self.genes)
        self.all = data.partition([self.genes], self.genes)
        self.single = data.partition([[g] for g in self.genes], self.genes)
        self.distance = np.array([[0,.1,.9,.9],[.1,0,.9,.9],[.9,.9,0,.1],[.9,.9,.1,0]])

    def constraints(self):
        return [dict(gene_a=a, gene_b=b, relation=r, weight=1, source='independent', block_id=block,
                     role='evaluation') for a,b,r,block in [('a','b','same','s1'), ('c','d','same','s2'),
                        ('a','c','different','d1'), ('b','d','different','d2')]]

    def args(self, extra=()):
        p = argparse.ArgumentParser()
        cli.add_parser(p.add_subparsers())
        return p.parse_args(['cluster', '-i', str(self.fa), '--gene-map', str(self.map),
                             '-o', str(self.root/'out')]+list(extra))

    def test_exact_membership_and_mapping(self):
        for groups in [[['a','b'],['b','c','d']], [['a','b','c']]]:
            with self.assertRaises(ValueError):
                data.partition(groups, self.genes)
        self.map.write_text('gene_ID\tassembly_ID\na\tA\na\tB\n')
        with self.assertRaises(ValueError):
            data.mapping(self.map, self.genes)
        self.map.write_text('gene_ID\tassembly_ID\na\tA\n')
        with self.assertRaises(ValueError):
            data.mapping(self.map, self.genes)

    def test_boundary_all_merged_split_missing(self):
        self.assertEqual(evidence.boundary(self.good, self.constraints())['B'], 1)
        self.assertEqual(evidence.boundary(self.all, self.constraints())['B'], 0)
        self.assertEqual(evidence.boundary(self.single, self.constraints())['B'], 0)
        self.assertIsNone(evidence.boundary(self.good, self.constraints()[:2])['B'])

    def test_block_balance(self):
        rows = self.constraints()
        rows.append(dict(rows[0], gene_a='a', gene_b='c', block_id='s3'))
        baseline = evidence.boundary(self.good, rows)['T_positive']
        rows += [dict(rows[-1]) for _ in range(100)]
        self.assertEqual(evidence.boundary(self.good, rows)['T_positive'], baseline)

    def test_conflicts_transitive_and_leakage(self):
        file = self.root/'constraints.tsv'
        rows = self.constraints()
        rows = [dict(rows[0], gene_a='a', gene_b='b'), dict(rows[0], gene_a='b',gene_b='c'),
                dict(rows[2], gene_a='a',gene_b='c')]
        data.tsv(file, rows, list(rows[0]))
        clean, conflicts = evidence.read_constraints(file, self.genes, 'hog', 'hog')
        self.assertEqual(clean, [])
        self.assertEqual(len(conflicts), 3)
        rows = self.constraints()
        rows[0]['role'] = 'construction'
        data.tsv(file, rows, list(rows[0]))
        clean, conflicts = evidence.read_constraints(file, self.genes, 'hog', 'hog')
        self.assertTrue(all(r['role']=='construction' for r in clean))
        self.assertEqual(len(conflicts), 3)
        with self.assertRaises(ValueError):
            evidence.read_constraints(file, self.genes, 'locus', 'hog')

    def test_quality_degenerate_negative_and_singleton(self):
        for labels in (self.all, self.single):
            q, _ = scoring.quality(labels, self.genes, self.distance)
            self.assertIsNone(q['Q'])
        wrong = data.partition([['a','c'],['b','d']], self.genes)
        q, _ = scoring.quality(wrong, self.genes, self.distance)
        self.assertLess(q['silhouette_mean'], 0)
        self.assertGreater(q['Q'], 0)
        mixed = data.partition([['a','b'],['c'],['d']], self.genes)
        q, sil = scoring.quality(mixed, self.genes, self.distance)
        self.assertEqual(sil['c'], 0)
        self.assertEqual(q['singleton_gene_fraction'], .5)
        q, _ = scoring.quality(self.good, self.genes, None)
        self.assertIsNone(q['Q'])

    def test_common_downgrade_zero_and_formal(self):
        def row(**kw):
            return dict(B=1, R=1, A=1, Q=1, positive_blocks=2, negative_blocks=2,
                        positive_gene_coverage=1, negative_gene_coverage=1, unresolved_fraction=0, **kw)
        rows = [row(), row()]
        scoring.totals(rows, scoring.WEIGHTS)
        self.assertEqual(rows[0]['status'], 'formal')
        rows[0]['R'] = None
        rows[1]['B'] = 0
        common, omitted = scoring.totals(rows, scoring.WEIGHTS)
        self.assertEqual(omitted, ['R'])
        self.assertEqual(rows[0]['score_config'], rows[1]['score_config'])
        self.assertEqual(rows[1]['total_score'], 0)
        self.assertEqual(rows[0]['status'], 'provisional')

    def test_ari_categories_no_parameter_votes(self):
        self.assertEqual(scoring.ari(self.good, self.good), 1)
        wrong = data.partition([['a','c'],['b','d']], self.genes)
        self.assertLess(scoring.ari(self.good, wrong), 0)
        c = dict(id='q', method='mmseqs', labels=self.good)
        a = dict(id='a', method='orthofinder', labels=self.good)
        b = dict(id='b', method='similarity-mcl', labels=wrong)
        score, raw = scoring.agreement(c, [c,a,b])
        self.assertEqual(score, .5)
        score2, _ = scoring.agreement(c, [c,a]+[dict(b,id=str(i)) for i in range(10)])
        self.assertEqual(score, score2)
        self.assertIsNone(scoring.agreement(c, [c,dict(c,id='other')])[0])
        with self.assertRaises(ValueError):
            scoring.ari(self.good, {'a':'1'})

    def test_no_full_proteome_no_orthofinder(self):
        context = dict(args=self.args(), constraints=[])
        with patch.object(methods.shutil, 'which', return_value='/fake'):
            for m in ['orthofinder','orthofinder-mmseqs','orthofinder-cdhit']:
                self.assertIn('complete proteomes', methods.applicable(m, context))

    def test_combination_never_crosses_hog(self):
        context = dict(args=self.args(), seqs=self.seqs)
        with patch.object(methods, 'hogs', return_value=([['a','b'],['c','d']], [])), \
             patch.object(methods, 'sequence_groups', side_effect=lambda m,g,p,c,w: ([g],[])):
            labels, _, _ = methods.execute('orthofinder-mmseqs', {}, context, self.root)
        self.assertEqual(labels['a'], labels['b'])
        self.assertNotEqual(labels['a'], labels['c'])

    def test_similarity_symmetric_inclusive_and_isolates(self):
        sim = self.root/'hits.tsv'
        sim.write_text('b\ta\t80\t20\t0\t0\t1\t20\t1\t20\t0\t90\n'
                       'a\tb\t50\t20\t0\t0\t1\t20\t1\t20\t0\t50\n')
        args = self.args(['--similarity', str(sim)])
        class FakeRunner:
            def run(inner, cmd, work):
                self.assertIn('a\ta\t1', (work/'graph.abc').read_text())
                self.assertIn('d\td\t1', (work/'graph.abc').read_text())
                (work/'mcl.txt').write_text('a\tb\nc\nd\n')
        context = dict(args=args, seqs=self.seqs, constraints=[], config=cli.DEFAULTS,
                       manifest={}, runner=FakeRunner())
        self.assertEqual(methods.similarity_edges(context)[('a','b')][0], .8)
        groups, isolated, _ = methods.graph_groups('similarity-mcl', {'identity':.5,'coverage':.8,'inflation':1.5}, context,self.root)
        self.assertEqual(isolated, ['c','d'])
        self.assertEqual(set(data.partition(groups,self.genes)), set(self.genes))

    def test_missing_distance_not_maximum_and_msa(self):
        file = self.root/'dist.tsv'
        file.write_text('gene_a\tgene_b\tdistance\na\tb\t0.1\n')
        d, info = data.distance_matrix(self.genes,self.seqs,str(file))
        self.assertIsNone(d)
        self.assertFalse(info['complete'])
        d, info = data.distance_matrix(self.genes,self.seqs,alignment=str(self.fa))
        self.assertTrue(info['complete'])
        self.assertEqual(float(d.sum()), 0)

    def test_weighted_graph_unknown_synteny_neutral(self):
        args = self.args()
        pairfile = self.root/'synteny.tsv'
        pairfile.write_text('a\tb\n')
        args.collinear_pairs = str(pairfile)
        class FakeRunner:
            def run(inner, cmd, work):
                (work/'mcl.txt').write_text('a\tb\nc\nd\n')
        context = dict(args=args, seqs=self.seqs, constraints=[], config=cli.DEFAULTS,
                       manifest={}, runner=FakeRunner(), edges={('a','b'):(.8,.8,1),('a','c'):(.8,.8,1)})
        params = dict(identity=.5,coverage=.8,inflation=1.5)
        methods.graph_groups('weighted-mcl', params, context, self.root)
        graph = (self.root/'graph.abc').read_text()
        self.assertIn('a\tb\t1.6', graph)
        self.assertIn('a\tc\t0.8', graph)
        methods.graph_groups('similarity-mcl', params, context, self.root)
        self.assertIn('a\tb\t0.8', (self.root/'graph.abc').read_text())

    def test_adapter_parameters_and_missing_member_rejected(self):
        recorded = []
        class FakeRunner:
            def run(inner, cmd, work):
                recorded.append(cmd)
                (work/'cluster_cluster.tsv').write_text('a\ta\na\tb\na\tc\n')
        context=dict(args=self.args(),seqs=self.seqs,runner=FakeRunner())
        with self.assertRaises(ValueError):
            methods.sequence_groups('mmseqs',self.genes,dict(identity=.8,coverage=.8),context,self.root)
        cmd=recorded[0]
        self.assertEqual(cmd[cmd.index('--cov-mode')+1],'0')
        self.assertEqual(cmd[cmd.index('--alignment-mode')+1],'3')

    def test_runner_real_subprocess_logging(self):
        # Real Python subprocess validates process/log plumbing, not an MCL executable.
        manifest={'commands':[]}
        runner=methods.Runner(manifest,self.root,cli.DEFAULTS)
        runner.run([sys.executable,'-c','print("integration probe")'],self.root)
        self.assertEqual(manifest['commands'][0]['returncode'],0)
        self.assertIn('integration probe',pathlib.Path(manifest['commands'][0]['log']).read_text())

    def test_grid_round_robin_and_duplicate_candidates(self):
        config=dict(cli.DEFAULTS,grid={'mmseqs':[{'identity':.8},{'identity':.9},{'identity':.8}]})
        plans=cli.candidates_for(['mmseqs','similarity-mcl'],self.args(),config)
        self.assertEqual([m for m,p in plans],['mmseqs','similarity-mcl','mmseqs'])

    def test_scheduler_failure_resume_reproducibility(self):
        def fake(method, params, context, work):
            if method == 'cdhit':
                raise RuntimeError('simulated tool failure')
            return self.good, [], ['mock adapter']
        args = self.args()
        with patch.object(cli, 'tool_versions', return_value={}), \
             patch.object(cli, 'applicable', side_effect=lambda m,c: None if m in ('mmseqs','cdhit','similarity-mcl') else 'not available'), \
             patch.object(cli, 'execute', side_effect=fake) as run:
            first = cli.run(args)
            self.assertEqual(run.call_count, 3)
            saved = (self.root/'out'/'selected_clusters.tsv').read_bytes()
            status = data.read_tsv(self.root/'out'/'method_status.tsv')
            self.assertEqual(next(r['status'] for r in status if r['method']=='cdhit'), 'failed')
            self.assertEqual(first['score_status'], 'provisional')
            self.assertEqual(len(first['tied_candidates']), 2)
            self.assertIsNone(first['selected_candidate'])
            self.assertEqual(data.read_tsv(self.root/'out'/'selected_clusters.tsv'), [])
            run.reset_mock()
            args.resume = True
            second = cli.run(args)
            self.assertEqual(run.call_count, 1)  # retry failed adapter only
            self.assertEqual(saved, (self.root/'out'/'selected_clusters.tsv').read_bytes())
            self.assertEqual(first['tied_candidates'], second['tied_candidates'])
            self.fa.write_text(self.fa.read_text()+'\n')
            with self.assertRaisesRegex(ValueError, 'fingerprint changed'):
                cli.run(args)

    def test_actual_hog_format_import_integration(self):
        # Real Python HOG reader, fixture data; no simulated external process here.
        result = self.root/'Results'
        (result/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
        (result/'Log.txt').write_text('OrthoFinder run completed')
        (result/'Phylogenetic_Hierarchical_Orthogroups'/'N0.tsv').write_text(
            'HOG\tOG\tGene Tree Parent Clade\tasma\tasmb\tasmc\tasmd\n'
            'N0.HOG1\tOG1\tn1\ta\tb\t\t\nN0.HOG2\tOG1\tn2\t\t\tc\t\n')
        args = self.args(['-M','orthofinder','--orthofinder-results',str(result)])
        with patch.object(cli, 'tool_versions', return_value={}):
            selection = cli.run(args)
        table = data.read_tsv(self.root/'out'/'selected_clusters.tsv')
        self.assertEqual({r['gene_ID'] for r in table}, set(self.genes))
        self.assertEqual(selection['status'], 'provisional')
        manifest=json.loads((self.root/'out'/'manifest.json').read_text())
        self.assertEqual(manifest['commands'], [])
        self.assertEqual(manifest['orthofinder_source']['unresolved_targets'], ['d'])


if __name__ == '__main__':
    unittest.main()
