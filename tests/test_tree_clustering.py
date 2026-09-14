"""Distance, internal-validity and CLI contracts with real Newick/SciPy parsing."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data, methods, scoring
from cluster_engine.tree import tree_distances, complete_linkage_groups


class TreeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.tree = self.root/'family.nwk'
        self.tree.write_text('((a:0.05,b:0.05)99:0.4,(c:0.05,d:0.05)0.9:0.4);')
        self.genes = list('abcd')
        self.good = data.partition([['a','b'], ['c','d']], self.genes)
        self.wrong = data.partition([['a','c'], ['b','d']], self.genes)
        self.scope, self.d, self.info = tree_distances(self.tree, self.genes)
        self.fa, self.map = self.root/'family.fa', self.root/'map.tsv'
        data.write_fasta(self.fa, {g:'M'*20 for g in 'abcde'})
        data.tsv(self.map, [dict(gene_ID=g, assembly_ID='asm'+g) for g in 'abcde'], ['gene_ID','assembly_ID'])

    def args(self, extra=()):
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        return parser.parse_args(['cluster','-i',str(self.fa),'--gene-map',str(self.map),
            '--gene-tree',str(self.tree),'-o',str(self.root/'out')]+list(extra))

    def candidates(self, pairs):
        return [dict(id=str(i), method=method, labels=labels, params={}, unsupported=[], runtime_seconds=0)
                for i,(method,labels) in enumerate(pairs)]

    def test_expected_patristic_distances_and_rerooting(self):
        expected = np.array([[0,.1,.9,.9],[.1,0,.9,.9],[.9,.9,0,.1],[.9,.9,.1,0]])
        np.testing.assert_allclose(self.d, expected)
        # Same unrooted tree represented with a trifurcation and arbitrary stem.
        self.tree.write_text('(b:0.05,(d:0.05,c:0.05):0.8,a:0.05):42;')
        _, rerooted, _ = tree_distances(self.tree, self.genes)
        np.testing.assert_allclose(rerooted, self.d)

    def test_quotes_comments_supports_zero_edges_and_exact_ids(self):
        self.tree.write_text("[&U](('a_b':0,'NA':0,'it''s':0)87[comment]:0.2,z:3e-2);")
        scope, d, info = tree_distances(self.tree, ['NA', 'a_b', "it's", 'z', 'missing'])
        self.assertEqual(scope, ['NA', 'a_b', "it's", 'z'])
        self.assertEqual(d[0,1], 0)
        self.assertAlmostEqual(d[0,3], .23)
        self.assertEqual(info['coverage'], .8)
        self.assertEqual(info['missing_target_genes'], ['missing'])
        q,_ = scoring.quality(dict(zip(scope,[0,0,0,1])), scope, d)
        self.assertEqual(q['silhouette_mean'], .75)  # singleton contributes zero
        self.assertEqual(q['dunn_status'], 'unbounded_positive_separation')
        self.assertEqual(q['dunn_bounded'], 1)

    def test_invalid_trees_and_budget(self):
        for newick in ['(a:1,a:2);', '(a,b:2);', '(a:-1,b:2);',
                       '(a:1,b:2);(c:1,d:2);', '(a:nan,b:1);']:
            self.tree.write_text(newick)
            with self.subTest(newick=newick), self.assertRaises(ValueError):
                tree_distances(self.tree, self.genes)
        self.tree.write_text('(a:1,b:1,c:1,d:1);')
        scope, d, info = tree_distances(self.tree, self.genes, max_genes=3)
        self.assertIsNone(d)
        self.assertIn('budget',info['reason'])
        with self.assertRaisesRegex(ValueError, 'no leaf IDs matching'):
            tree_distances(self.tree, ['other'])

    def test_complete_linkage_cut_and_order_invariance(self):
        groups,z = complete_linkage_groups(self.genes, self.d, .2)
        self.assertEqual(groups, [['a','b'],['c','d']])
        perm = [3,2,0,1]
        reverse,_ = complete_linkage_groups([self.genes[i] for i in perm],self.d[np.ix_(perm,perm)], .2)
        self.assertEqual(groups,reverse)
        self.assertEqual(len(complete_linkage_groups(self.genes,self.d,1,z)[0]),1)
        self.assertEqual(len(complete_linkage_groups(self.genes,self.d,0,z)[0]),4)
        for g in groups:
            ix = [self.genes.index(x) for x in g]
            self.assertLessEqual(self.d[np.ix_(ix,ix)].max(), .2)

    def test_hand_calculated_indices_and_scale_invariance(self):
        q,_ = scoring.quality(self.good, self.genes, self.d)
        self.assertAlmostEqual(q['silhouette_mean'], 8/9)
        self.assertAlmostEqual(q['dunn_index'], 9)
        self.assertAlmostEqual(q['dunn_bounded'], .9)
        scaled,_ = scoring.quality(self.good,self.genes,20*self.d)
        for k in ['silhouette_mean','dunn_index','dunn_bounded']:
            self.assertAlmostEqual(q[k],scaled[k])
        wrong,_ = scoring.quality(self.wrong,self.genes,self.d)
        self.assertLess(wrong['silhouette_mean'],0)
        for labels in [data.partition([self.genes], self.genes),dict(zip(self.genes,self.genes))]:
            q,_=scoring.quality(labels,self.genes,self.d)
            self.assertIsNone(q['silhouette_mean'])
            self.assertIsNone(q['dunn_bounded'])
        q,_=scoring.quality(self.good,self.genes,np.zeros_like(self.d))
        self.assertEqual(q['silhouette_mean'],0)
        self.assertIsNone(q['dunn_bounded'])
        self.assertEqual(q['dunn_status'],'undefined_zero_over_zero')

    def test_invalid_evaluation_matrix(self):
        for d in [self.d[:3,:3], self.d+np.eye(4), -self.d,
                  np.full((4,4),np.nan), np.triu(self.d)]:
            with self.assertRaises(ValueError):
                scoring.quality(self.good,self.genes,d)

    def test_identical_partitions_equal_scores_and_no_agreement_fallback(self):
        renamed={g:'label'+str(self.good[g]) for g in self.genes}
        pairs=[('mmseqs',self.good),('similarity-mcl',renamed),('cdhit',self.wrong),
               ('tree',dict(zip(self.genes,self.genes)))]
        c=self.candidates(pairs)
        rows,*_=scoring.evaluate(c,self.genes,self.d,[],cli.DEFAULTS)
        self.assertEqual(rows[0]['total_score'],rows[1]['total_score'])
        self.assertEqual(rows[0]['partition_signature'],rows[1]['partition_signature'])
        self.assertGreater(rows[0]['total_score'],rows[2]['total_score'])
        self.assertIsNone(rows[3]['total_score'])
        self.assertTrue(rows[0]['ranking_eligible'])
        more,*_=scoring.evaluate(c+c[:1],self.genes,self.d,[],cli.DEFAULTS)
        self.assertEqual(rows[0]['total_score'],more[0]['total_score'])
        unavailable,*_=scoring.evaluate(c,self.genes,None,[],cli.DEFAULTS,[])
        self.assertTrue(all(r['total_score'] is None for r in unavailable))
        dunn,*_=scoring.evaluate(c,self.genes,self.d,[],dict(cli.DEFAULTS,ranking_metric='dunn'))
        self.assertAlmostEqual(dunn[0]['total_score'],90)

    def test_fixed_scope_coverage_and_singleton_convention(self):
        all_genes = self.genes+['e','f']
        a = dict(self.good,e='extra',f='extra')
        b = dict(self.good,e=self.good['a'],f='unique')
        c=self.candidates([('tree',a),('mmseqs',b)])
        rows,sil,*_=scoring.evaluate(c,all_genes,self.d,[],cli.DEFAULTS,self.genes)
        self.assertEqual(rows[0]['total_score'],rows[1]['total_score'])
        self.assertNotEqual(rows[0]['partition_signature'],rows[1]['partition_signature'])
        self.assertEqual(rows[0]['evaluation_partition_signature'],rows[1]['evaluation_partition_signature'])
        self.assertEqual(rows[0]['evaluation_coverage'],4/6)
        self.assertEqual(set(sil['0']),set(self.genes))
        blocked,*_=scoring.evaluate(c,all_genes,self.d,[],dict(cli.DEFAULTS,min_evaluation_coverage=.8),self.genes)
        self.assertFalse(blocked[0]['ranking_eligible'])
        self.assertIsNotNone(blocked[0]['total_score'])

    def test_real_tree_cli_outputs_missing_genes_and_resume(self):
        args=self.args(['-M','tree','--tree-threshold','.2'])
        with patch.object(cli,'tool_versions',return_value={}):
            selection=cli.run(args)
            self.assertEqual(selection['status'],'selected')
            self.assertAlmostEqual(selection['best_score'],50*(1+8/9))
            self.assertEqual(selection['evaluation_scope']['coverage'],.8)
            table=data.read_tsv(self.root/'out'/'selected_clusters.tsv')
            self.assertEqual(len(table),5)
            self.assertEqual(len({r['cluster_ID'] for r in table}),3)
            self.assertEqual(data.read_tsv(self.root/'out'/'scores.tsv')[0]['assignment_coverage'],'0.8')
            args.resume=True
            self.assertEqual(cli.run(args)['best_score'],selection['best_score'])
            self.tree.write_text(self.tree.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'fingerprint changed'):
                cli.run(args)

    def test_auto_tree_without_proteomes_and_equivalent_methods(self):
        labels=dict(self.good,e='extra')
        def execute(method,params,context,work):
            if method=='tree':
                return methods.execute(method,params,context,work)
            return labels,[],['mock sequence adapter']
        with patch.object(cli,'tool_versions',return_value={}), \
             patch.object(cli,'applicable',side_effect=lambda m,c: None if m in ('tree','mmseqs') else 'unavailable'), \
             patch.object(cli,'execute',side_effect=execute):
            selection=cli.run(self.args(['--tree-threshold','.2']))
        self.assertEqual(selection['status'],'equivalent_best')
        self.assertEqual(selection['recommended_methods'],['mmseqs','tree'])
        self.assertEqual(len(data.read_tsv(self.root/'out'/'selected_clusters.tsv')),5)

    def test_tie_on_evaluation_scope_does_not_hide_different_full_partitions(self):
        def execute(method,params,context,work):
            return dict(self.good,e='extra' if method=='mmseqs' else self.good['a']),[],[]
        with patch.object(cli,'tool_versions',return_value={}), \
             patch.object(cli,'applicable',side_effect=lambda m,c: None if m in ('cdhit','mmseqs') else 'unavailable'), \
             patch.object(cli,'execute',side_effect=execute):
            selection=cli.run(self.args())
        self.assertEqual(selection['status'],'ambiguous')
        self.assertIsNone(selection['selected_candidate'])
        self.assertEqual(data.read_tsv(self.root/'out'/'selected_clusters.tsv'),[])

    def test_explicit_evaluation_distances_override_tree_and_tree_is_unchanged(self):
        path=self.root/'dist.tsv'
        genes=list('abcde')
        data.tsv(path,[dict(gene_a=a,gene_b=b,distance=.5) for i,a in enumerate(genes) for b in genes[:i]],
                 ['gene_a','gene_b','distance'])
        with patch.object(cli,'tool_versions',return_value={}):
            result=cli.run(self.args(['-M','tree','--tree-threshold','.2',
                '--evaluation-distances',str(path),'--distance-provenance','synthetic independent fixture']))
        self.assertEqual(result['evaluation_scope']['coverage'],1)
        self.assertEqual(result['best_score'],50)
        manifest=json.loads((self.root/'out'/'manifest.json').read_text())
        self.assertEqual(manifest['distance']['source'],'external-distances')
        rows=data.read_tsv(self.root/'out'/'selected_clusters.tsv')
        self.assertEqual(len({r['cluster_ID'] for r in rows}),3)


if __name__=='__main__':
    unittest.main()
