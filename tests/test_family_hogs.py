"""Family HOG inference contracts using real tree reconciliation and no executables."""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data, hog_tree


def memberships(groups):
    return {frozenset(group) for group in groups}


class FamilyHogTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.gene_tree = self.root / 'family.nwk'
        self.species_tree = self.root / 'species.nwk'
        self.output = self.root / 'hogs'
        self.mapping = {'o': 'O', 'a1': 'A', 'b1': 'B', 'a2': 'A', 'b2': 'B'}
        self.gene_tree.write_text('(o:1,((a1:1,b1:1):1,(a2:1,b2:1):1):1);', encoding='utf-8')
        self.species_tree.write_text('(O:1,(A:1,B:1):1);', encoding='utf-8')

    def infer(self, **kwargs):
        params = dict(gene_tree=self.gene_tree, species_tree=self.species_tree,
                      gene_to_assembly=self.mapping, outgroup='O', output_dir=self.output)
        params.update(kwargs)
        return hog_tree.infer_hogs(**params)

    def assert_complete(self, result, mapping=None):
        """Every target is retained exactly once, even without selected-node support."""
        target = set(self.mapping if mapping is None else mapping)
        grouped = [gene for group in result['groups'] for gene in group]
        self.assertEqual(len(grouped), len(set(grouped)))
        self.assertEqual(len(result['unsupported']), len(set(result['unsupported'])))
        self.assertEqual(set(grouped), target)
        self.assertTrue(set(result['unsupported']).issubset(target))
        for group in result['groups']:
            if set(group) & set(result['unsupported']):
                self.assertEqual(len(group), 1)

    def test_ancestral_duplication_splits_root_hogs(self):
        self.gene_tree.write_text('((A1:1,B1:1):1,(A2:1,B2:1):1);', encoding='utf-8')
        self.species_tree.write_text('(A:1,B:1);', encoding='utf-8')
        mapping = {'A1': 'A', 'B1': 'B', 'A2': 'A', 'B2': 'B'}
        result = self.infer(gene_to_assembly=mapping, outgroup='A')
        self.assertEqual(memberships(result['groups']), memberships([['A1', 'B1'], ['A2', 'B2']]))
        self.assertEqual(result['unsupported'], [])
        self.assertEqual(result['selected_level'], 'N0')
        self.assertEqual(set(result['selected_scope']), {'A', 'B'})
        self.assert_complete(result, mapping)

    def test_duplication_changes_membership_at_descendant_level(self):
        root_result = self.infer()
        self.assertEqual(memberships(root_result['groups']), {frozenset(self.mapping)})
        self.assertEqual(root_result['unsupported'], [])
        descendant = self.infer(hog_level='N1', output_dir=self.root / 'descendant')
        supported = [group for group in descendant['groups'] if set(group) - set(descendant['unsupported'])]
        self.assertEqual(memberships(supported), memberships([['a1', 'b1'], ['a2', 'b2']]))
        self.assertEqual(descendant['unsupported'], ['o'])
        self.assertEqual(set(descendant['selected_scope']), {'A', 'B'})
        self.assert_complete(descendant)

    def test_terminal_duplication_remains_in_ancestral_hog(self):
        self.gene_tree.write_text('(o:1,((a1:1,a2:1):1,b:1):1);', encoding='utf-8')
        mapping = {'o': 'O', 'a1': 'A', 'a2': 'A', 'b': 'B'}
        result = self.infer(gene_to_assembly=mapping, hog_level='N1')
        supported = [group for group in result['groups'] if set(group) - set(result['unsupported'])]
        self.assertEqual(memberships(supported), memberships([['a1', 'a2', 'b']]))
        self.assertEqual(result['unsupported'], ['o'])
        self.assert_complete(result, mapping)

    def test_two_gene_family_and_outgroup_absent_from_gene_tree(self):
        self.gene_tree.write_text('(a:1,b:1);', encoding='utf-8')
        mapping = {'a': 'A', 'b': 'B'}
        result = self.infer(gene_to_assembly=mapping, hog_level='N1')
        self.assertEqual(memberships(result['groups']), memberships([['a', 'b']]))
        self.assertEqual(result['unsupported'], [])
        self.assertEqual(result['outgroup_gene_count'], 0)
        self.assert_complete(result, mapping)

    def test_optional_paralog_split_changes_hogs_above_missing_outgroup(self):
        self.gene_tree.write_text('((a1:1,b1:1):1,(a2:1,b2:1):1);', encoding='utf-8')
        mapping = {gene: assembly for gene, assembly in self.mapping.items() if gene != 'o'}
        default = self.infer(gene_to_assembly=mapping)
        split = self.infer(gene_to_assembly=mapping, split_paralogous_clades=True,
                           output_dir=self.root / 'extra_split')
        self.assertEqual(memberships(default['groups']), {frozenset(mapping)})
        self.assertEqual(memberships(split['groups']), memberships([['a1', 'b1'], ['a2', 'b2']]))
        self.assertFalse(default['split_paralogous_clades'])
        self.assertTrue(split['split_paralogous_clades'])
        self.assert_complete(split, mapping)

    def test_unrooted_binary_representation_accepted_but_polytomy_rejected(self):
        self.gene_tree.write_text('(o:1,a:1,b:1);', encoding='utf-8')
        mapping = {'o': 'O', 'a': 'A', 'b': 'B'}
        result = self.infer(gene_to_assembly=mapping)
        self.assertEqual(memberships(result['groups']), {frozenset(mapping)})
        self.gene_tree.write_text('(o:1,(a1:1,b1:1,a2:1,b2:1):1);', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'polytom'):
            self.infer(output_dir=self.root / 'polytomy')

    def test_equivalent_unrooted_trees_and_input_order_preserve_hogs(self):
        reference = self.infer(hog_level='N1')
        self.gene_tree.write_text('((b2:1,a2:1):1,o:2,(b1:1,a1:1):1);', encoding='utf-8')
        self.species_tree.write_text('(B:1,O:2,A:1);', encoding='utf-8')
        reordered = dict(reversed(list(self.mapping.items())))
        result = self.infer(gene_to_assembly=reordered, hog_level='N1', output_dir=self.root / 'reordered')
        self.assertEqual(memberships(result['groups']), memberships(reference['groups']))
        self.assertEqual(result['unsupported'], reference['unsupported'])
        self.assertEqual(result['selected_scope'], reference['selected_scope'])
        self.assert_complete(result)

    def test_gene_names_use_explicit_map_and_preserve_punctuation(self):
        mapping = {'unrelated|alpha.1': 'Assembly-A', 'second_gene-v2': 'Assembly-B',
                   'outside|isoform.3': 'Outgroup-assembly'}
        self.gene_tree.write_text('(outside|isoform.3:1,(unrelated|alpha.1:1,second_gene-v2:1):1);', encoding='utf-8')
        self.species_tree.write_text('(Outgroup-assembly,(Assembly-A,Assembly-B));', encoding='utf-8')
        before = dict(mapping)
        result = self.infer(gene_to_assembly=mapping, outgroup='Outgroup-assembly')
        self.assertEqual(mapping, before)
        self.assertEqual(memberships(result['groups']), {frozenset(mapping)})
        tree = Phylo.read(str(self.output / 'GeneTree_rooted.nwk'), 'newick')
        self.assertEqual({tip.name for tip in tree.get_terminals()}, set(mapping))
        self.assert_complete(result, mapping)

    def test_iqtree_and_fasttree_support_labels(self):
        for label in ['99.8/100', '0.999', '100']:
            with self.subTest(label=label):
                self.gene_tree.write_text('(o:1,((a1:1,b1:1){}:1,(a2:1,b2:1){}:1){}:1);'.format(label, label, label), encoding='utf-8')
                self.species_tree.write_text('(O:1,(A:1,B:1){}:1);'.format(label), encoding='utf-8')
                result = self.infer(output_dir=self.root / ('support_' + label.replace('/', '_')))
                self.assertEqual(memberships(result['groups']), {frozenset(self.mapping)})
                self.assert_complete(result)

    def test_original_gene_confidence_is_audited_on_original_clades(self):
        self.gene_tree.write_text('(o:1,((a1:1,b1:1)97:1,(a2:1,b2:1)94:1)86:1)83;', encoding='utf-8')
        self.infer()
        rows = data.read_tsv(self.output / 'input_gene_node_labels.tsv')
        labels = {frozenset(json.loads(row['genes'])): float(row['input_label']) for row in rows}
        self.assertEqual(labels, {frozenset(self.mapping): 83,
                                  frozenset(['a1', 'b1', 'a2', 'b2']): 86,
                                  frozenset(['a1', 'b1']): 97,
                                  frozenset(['a2', 'b2']): 94})
        for name in ['GeneTree_rooted.nwk', 'GeneTree_reconciled.nwk']:
            tree = Phylo.read(str(self.output / name), 'newick')
            self.assertTrue(all(node.confidence is None for node in tree.get_nonterminals()))
        reconciled = Phylo.read(str(self.output / 'GeneTree_reconciled.nwk'), 'newick')
        self.assertTrue(all(node.name.startswith('n') for node in reconciled.get_nonterminals()))

    def test_species_rerooting_does_not_reassign_old_clade_confidence(self):
        self.species_tree.write_text('((O:1,A:1)97:1,B:1)55;', encoding='utf-8')
        self.infer()
        node_map = {row['node']: row for row in data.read_tsv(self.output / 'species_node_map.tsv')}
        self.assertEqual(set(json.loads(node_map['N1']['assemblies'])), {'A', 'B'})
        self.assertEqual(node_map['N1']['input_label'], '')
        original = data.read_tsv(self.output / 'input_species_node_labels.tsv')
        old_clade = next(row for row in original if set(json.loads(row['assemblies'])) == {'O', 'A'})
        self.assertEqual(float(old_clade['input_label']), 97)
        original_root = next(row for row in original if set(json.loads(row['assemblies'])) == {'O', 'A', 'B'})
        self.assertEqual(float(original_root['input_label']), 55)

    def test_reroots_multispecies_outgroup_on_complementary_edge(self):
        self.species_tree.write_text('(O1:1,((A:1,B:1):1,O2:1):1);', encoding='utf-8')
        self.gene_tree.write_text('((a:1,b:1):1,(o1:1,o2:1):1);', encoding='utf-8')
        mapping = {'a': 'A', 'b': 'B', 'o1': 'O1', 'o2': 'O2'}
        result = self.infer(gene_to_assembly=mapping, outgroup=['O1', 'O2'])
        rooted = Phylo.read(str(self.output / 'SpeciesTree_rooted_node_labels.txt'), 'newick')
        splits = [{tip.name for tip in clade.get_terminals()} for clade in rooted.root.clades]
        self.assertCountEqual(splits, [{'O1', 'O2'}, {'A', 'B'}])
        ingroup = next(clade for clade in rooted.root.clades if {tip.name for tip in clade.get_terminals()} == {'A', 'B'})
        self.assertEqual(rooted.root.name, 'N0')
        self.assertEqual(ingroup.name, 'N1')
        self.assert_complete(result, mapping)

    def test_required_known_nonexhaustive_monophyletic_outgroup(self):
        for outgroup in [None, '', [], 'unknown', ['O', 'A', 'B']]:
            with self.subTest(outgroup=outgroup), self.assertRaises(ValueError):
                self.infer(outgroup=outgroup)
        self.species_tree.write_text('((O1:1,A:1):1,(O2:1,B:1):1);', encoding='utf-8')
        mapping = {'a1': 'A', 'b1': 'B', 'a2': 'A', 'b2': 'B', 'o': 'O1'}
        with self.assertRaises(ValueError):
            self.infer(gene_to_assembly=mapping, outgroup=['O1', 'O2'])

    def test_mapping_and_tree_tips_must_match_exactly(self):
        for mapping in [dict(self.mapping, absent='A'), {k: v for k, v in self.mapping.items() if k != 'a1'},
                        dict(self.mapping, a1='unknown')]:
            with self.subTest(mapping=mapping), self.assertRaises(ValueError):
                self.infer(gene_to_assembly=mapping)

    def test_duplicate_gene_and_species_tips_rejected(self):
        self.gene_tree.write_text('(o:1,((a1:1,b1:1):1,(a1:1,b2:1):1):1);', encoding='utf-8')
        with self.assertRaises(ValueError):
            self.infer()
        self.gene_tree.write_text('(o:1,((a1:1,b1:1):1,(a2:1,b2:1):1):1);', encoding='utf-8')
        self.species_tree.write_text('(O:1,(A:1,B:1,A:1):1);', encoding='utf-8')
        with self.assertRaises(ValueError):
            self.infer()

    def test_invalid_gene_lengths_and_multiple_trees_rejected(self):
        for newick in ['(o,(a1:1,b1:1,a2:1,b2:1):1);',
                       '(o:-1,(a1:1,b1:1,a2:1,b2:1):1);',
                       '(o:nan,(a1:1,b1:1,a2:1,b2:1):1);',
                       '(o:1,(a1:1,b1:1,a2:1,b2:1):1);(extra:1,more:1);']:
            self.gene_tree.write_text(newick, encoding='utf-8')
            with self.subTest(newick=newick), self.assertRaises(ValueError):
                self.infer()

    def test_unknown_hog_level_rejected(self):
        for level in ['N99', 'A', '']:
            with self.subTest(level=level), self.assertRaises(ValueError):
                self.infer(hog_level=level)

    def test_budget_callback_interrupts_inference(self):
        calls = []

        def exhausted():
            calls.append(True)
            raise TimeoutError('test budget exhausted')

        with self.assertRaisesRegex(TimeoutError, 'test budget exhausted'):
            self.infer(check_budget=exhausted)
        self.assertEqual(calls, [True])

    def test_artifacts_retain_original_identifiers_and_hierarchy(self):
        result = self.infer()
        expected = ['SpeciesTree_rooted_node_labels.txt', 'GeneTree_rooted.nwk',
                    'GeneTree_reconciled.nwk', 'species_node_map.tsv', 'gene_id_map.tsv',
                    'gene_events.tsv', 'hog_members.tsv', 'inference.json',
                    'Phylogenetic_Hierarchical_Orthogroups/N0.tsv',
                    'Phylogenetic_Hierarchical_Orthogroups/N1.tsv']
        for name in expected:
            with self.subTest(artifact=name):
                path = self.output / name
                self.assertTrue(path.is_file(), str(path))
                self.assertGreater(path.stat().st_size, 0)
        with (self.output / 'Phylogenetic_Hierarchical_Orthogroups/N0.tsv').open(encoding='utf-8', newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 1)
        self.assertEqual(set(rows[0]) - {'HOG', 'OG', 'Gene Tree Parent Clade'}, {'A', 'B', 'O'})
        exported = {gene.strip() for assembly in ('A', 'B', 'O') for gene in rows[0][assembly].split(',') if gene.strip()}
        self.assertEqual(exported, set(self.mapping))
        metadata = json.loads((self.output / 'inference.json').read_text(encoding='utf-8'))
        self.assertEqual(metadata['selected_level'], result['selected_level'])
        self.assertEqual(set(metadata['selected_scope']), set(result['selected_scope']))
        events = (self.output / 'gene_events.tsv').read_text(encoding='utf-8')
        self.assertIn('N1', events)
        self.assert_complete(result)

    def test_cli_accepts_hog_tree_with_species_tree_alias(self):
        for species_option in ['--tree', '--species-tree']:
            parser = argparse.ArgumentParser()
            cli.add_parser(parser.add_subparsers())
            args = parser.parse_args(['cluster', '-M', 'hog-tree', '-i', 'family.fa',
                                     '--gene-map', 'gene_map.tsv', '--gene-tree', 'family.nwk',
                                     species_option, 'species.treefile', '--outgroup', 'O',
                                     '--hog-level', 'N1', '-o', 'report'])
            self.assertEqual(args.method, 'hog-tree')
            self.assertEqual(args.tree, 'species.treefile')
            self.assertEqual(args.gene_tree, 'family.nwk')
            self.assertEqual(args.hog_level, 'N1')
            self.assertIn('O', args.outgroup)

    def test_cli_requires_explicit_species_tree_gene_tree_and_outgroup(self):
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        required = ['--tree', 'species.treefile', '--gene-tree', 'family.nwk', '--outgroup', 'O']
        for index in range(0, len(required), 2):
            args = parser.parse_args(['cluster', '-M', 'hog-tree', '-i', 'family.fa',
                                      '--gene-map', 'map.tsv', '-o', 'report'] +
                                     required[:index] + required[index + 2:])
            with self.subTest(missing=required[index]), self.assertRaisesRegex(ValueError, 'outgroup'):
                cli.configure(args)

    def test_real_cluster_worker_reports_membership_and_reuses_cache(self):
        fasta = self.root / 'family.fa'
        mapping = self.root / 'map.tsv'
        report = self.root / 'report'
        data.write_fasta(fasta, {gene: 'MPEPTIDE' for gene in self.mapping})
        data.tsv(mapping, [dict(gene_ID=gene, assembly_ID=assembly)
                           for gene, assembly in self.mapping.items()], ['gene_ID', 'assembly_ID'])
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        args = parser.parse_args(['cluster', '-M', 'hog-tree', '-i', str(fasta),
                                 '--gene-map', str(mapping), '--gene-tree', str(self.gene_tree),
                                 '--species-tree', str(self.species_tree), '--outgroup', 'O',
                                 '--hog-level', 'N1', '-o', str(report)])
        with patch.object(cli, 'tool_versions', return_value={}):
            selection = cli.run(args)
            self.assertEqual(selection['status'], 'selected')
            rows = data.read_tsv(report / 'selected_clusters.tsv')
            groups = {}
            for row in rows:
                groups.setdefault(row['cluster_ID'], []).append(row['gene_ID'])
            self.assertEqual(memberships(groups.values()), memberships([['a1', 'b1'], ['a2', 'b2'], ['o']]))
            manifest = json.loads((report / 'manifest.json').read_text(encoding='utf-8'))
            self.assertEqual(manifest['candidate_status'][0]['status'], 'success')
            self.assertEqual(len(manifest['family_hog_runs']), 1)
            self.assertEqual(manifest['family_hog_runs'][0]['outgroup'], ['O'])
            self.assertTrue(any('hog_tree_worker.py' in arg for command in manifest['commands'] for arg in command['argv']))
            args.resume = True
            self.assertEqual(cli.run(args)['status'], 'selected')
            resumed = json.loads((report / 'manifest.json').read_text(encoding='utf-8'))
            self.assertTrue(resumed['candidate_status'][0]['cache'])
            self.assertEqual(len(resumed['commands']), len(manifest['commands']))
            self.assertEqual(len(resumed['family_hog_runs']), 1)
            hog_tables = list((report / 'candidates').rglob('N1.tsv'))
            self.assertEqual(len(hog_tables), 1)
            with hog_tables[0].open('a', encoding='utf-8') as handle:
                handle.write('\n')
            self.assertEqual(cli.run(args)['status'], 'failed')
            refused = json.loads((report / 'manifest.json').read_text(encoding='utf-8'))
            self.assertIn('HOG artifact cache integrity mismatch', refused['candidate_status'][0]['reason'])
            self.assertEqual(data.read_tsv(report / 'selected_clusters.tsv'), [])

    def test_incomplete_worker_attempt_is_preserved_and_resume_retries(self):
        fasta = self.root / 'family.fa'
        mapping = self.root / 'map.tsv'
        report = self.root / 'report'
        data.write_fasta(fasta, {gene: 'MPEPTIDE' for gene in self.mapping})
        data.tsv(mapping, [dict(gene_ID=gene, assembly_ID=assembly)
                           for gene, assembly in self.mapping.items()], ['gene_ID', 'assembly_ID'])
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        args = parser.parse_args(['cluster', '-M', 'hog-tree', '-i', str(fasta),
                                 '--gene-map', str(mapping), '--gene-tree', str(self.gene_tree),
                                 '--species-tree', str(self.species_tree), '--outgroup', 'O',
                                 '--hog-level', 'N1', '-o', str(report)])
        interrupted_attempts = []

        def interrupt_worker(command, work, **kwargs):
            request = json.loads(Path(command[-1]).read_text(encoding='utf-8'))
            attempt = Path(request['output_dir'])
            attempt.mkdir(parents=True)
            (attempt / 'partial_diagnostic.txt').write_text('partial attempt retained', encoding='utf-8')
            interrupted_attempts.append(attempt)
            raise TimeoutError('simulated interruption after partial output')

        with patch.object(cli, 'tool_versions', return_value={}):
            with patch.object(cli.Runner, 'run', side_effect=interrupt_worker):
                self.assertEqual(cli.run(args)['status'], 'failed')
            self.assertEqual(len(interrupted_attempts), 1)
            old_attempt = interrupted_attempts[0]
            self.assertFalse((old_attempt / 'inference.json').exists())
            args.resume = True
            self.assertEqual(cli.run(args)['status'], 'selected')
            self.assertEqual((old_attempt / 'partial_diagnostic.txt').read_text(encoding='utf-8'), 'partial attempt retained')
            complete = list((report / 'candidates').rglob('inference.json'))
            self.assertEqual(len(complete), 1)
            self.assertNotEqual(complete[0].parent, old_attempt)
            self.assertTrue(complete[0].parent.parent.samefile(old_attempt.parent))
            manifest = json.loads((report / 'manifest.json').read_text(encoding='utf-8'))
            self.assertEqual(manifest['candidate_status'][0]['status'], 'success')
            self.assertFalse(manifest['candidate_status'][0]['cache'])
            self.assertEqual(len(data.read_tsv(report / 'selected_clusters.tsv')), len(self.mapping))


if __name__ == '__main__':
    unittest.main()
