"""HOG alias, scope and failed-import regressions, without external programs."""
import argparse
import json
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data, hog_import, methods
from orthofinder_process import parse_hogs


class HogImportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.results = self.root/'export'
        (self.results/'Phylogenetic_Hierarchical_Orthogroups').mkdir(parents=True)
        (self.results/'Species_Tree').mkdir()
        self.tree = self.results/'Species_Tree/SpeciesTree_rooted_node_labels.txt'
        self.tree.write_text('(Out:1,(A_AGAT:1,B_AGAT:1)N1:1)N0;')
        self.hog = self.results/'Phylogenetic_Hierarchical_Orthogroups/N1.tsv'
        self.hog.write_text('HOG\tOG\tGene Tree Parent Clade\tA_AGAT\tB_AGAT\tOut\n'
                            'N1.HOG1\tOG1\tn1\ta.1, a.2\tb\t\n'
                            'N1.HOG2\tOG2\tn2\tc\t\t\n')
        self.assemblies = {'a.1': 'A', 'a.2': 'A', 'b': 'B', 'c': 'A', 'missing': 'A'}
        self.fa, self.map = self.root/'family.fa', self.root/'map.tsv'
        data.write_fasta(self.fa, {g: 'MPEPTIDE' for g in self.assemblies})
        data.tsv(self.map, [dict(gene_ID=g, assembly_ID=a) for g, a in self.assemblies.items()],
                 ['gene_ID', 'assembly_ID'])

    def args(self, extra=()):
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        return parser.parse_args(['cluster', '-M', 'orthofinder', '-i', str(self.fa),
            '--gene-map', str(self.map), '--orthofinder-results', str(self.results),
            '--hog-level', 'N1', '-o', str(self.root/'out')] + list(extra))

    def load(self, assemblies=None):
        table = parse_hogs(str(self.results), 'N1', detailed=True)
        scope = hog_import.node_scope(self.results, 'N1', table)
        return hog_import.import_targets(table, assemblies or self.assemblies, scope)

    def test_membership_alias_preserves_isoforms_and_missing_gene(self):
        before = dict(self.assemblies)
        groups, missing, info = self.load()
        self.assertEqual(self.assemblies, before)
        self.assertEqual(missing, ['missing'])
        self.assertIn(['a.1', 'a.2', 'b'], groups)
        self.assertEqual(info['assigned_gene_count'], 4)
        self.assertEqual(info['imported_HOG_count'], 2)
        self.assertEqual(info['import_status'], 'partial')
        self.assertEqual([r['orthofinder_assembly_ID'] for r in info['assembly_mapping']], ['A_AGAT', 'B_AGAT'])
        self.assertTrue(all(r['rule'] == 'unique_target_gene_membership' for r in info['assembly_mapping']))
        row = next(r for r in info['assignments'] if r['gene_ID'] == 'a.2')
        self.assertEqual((row['HOG'], row['OG'], row['gene_tree_parent_clade']), ('N1.HOG1', 'OG1', 'n1'))

    def test_empty_outgroup_column_and_scope(self):
        table = parse_hogs(str(self.results), 'N1', detailed=True)
        self.assertIn('Out', table.attrs['assembly_columns'])
        scope = hog_import.node_scope(self.results, 'N1', table)
        self.assertEqual(scope['leaves'], ['A_AGAT', 'B_AGAT'])
        _, missing, info = self.load(dict(self.assemblies, out_gene='Out'))
        self.assertIn('out_gene', missing)
        self.assertEqual(info['outside_node_input_assemblies'], ['Out'])
        row = next(r for r in info['assignments'] if r['gene_ID'] == 'out_gene')
        self.assertEqual(row['reason'], 'assembly_outside_selected_node')
        with self.hog.open('a') as handle:
            handle.write('N1.HOG3\tOG3\tn3\t\t\tout_gene\n')
        with self.assertRaisesRegex(ValueError, 'outside species-tree node'):
            self.load()

    def test_ambiguous_and_conflicting_aliases_rejected(self):
        self.hog.write_text('HOG\tOG\tGene Tree Parent Clade\tA_AGAT\tB_AGAT\n'
                            'N1.HOG1\tOG1\tn1\tshared, a\tshared, b\n')
        with self.assertRaisesRegex(ValueError, 'ambiguous HOG assembly mapping'):
            self.load({'shared': 'unknown'})
        with self.assertRaisesRegex(ValueError, 'conflicting HOG assembly mapping'):
            self.load({'a': 'unknown', 'b': 'unknown'})
        # Exact assembly names disambiguate a repeated source ID.
        _, missing, info = self.load({'shared': 'A_AGAT'})
        self.assertEqual(missing, [])
        self.assertEqual(info['assembly_mapping'][0]['rule'], 'exact_assembly_name')

    def test_alias_collision_rejected_and_exact_names_not_rewritten(self):
        with self.assertRaisesRegex(ValueError, 'multiple input assemblies'):
            self.load({'a.1': 'first', 'a.2': 'second'})
        _, missing, info = self.load({'a.1': 'B_AGAT', 'b': 'B_AGAT'})
        self.assertEqual(missing, ['a.1'])
        self.assertEqual(info['assignments'][0]['reason'], 'gene_not_in_mapped_assembly')
        self.assertEqual(info['assigned_gene_count'], 1)

    def test_no_prefix_guess_when_alias_has_no_anchors(self):
        _, missing, info = self.load({'a.1': 'A', 'absent': 'B'})
        self.assertEqual(missing, ['absent'])
        self.assertEqual(info['unmapped_input_assemblies'], ['B'])
        self.assertEqual(info['assembly_mapping'][1]['orthofinder_assembly_ID'], '')

    def test_node_must_exist_and_export_requires_tree(self):
        self.tree.write_text('(Out:1,(A_AGAT:1,B_AGAT:1)N2:1)N0;')
        with self.assertRaisesRegex(ValueError, 'exactly once'):
            self.load()
        self.tree.unlink()
        table = parse_hogs(str(self.results), 'N1')
        self.assertEqual(hog_import.node_scope(self.results, 'N1', table)['status'], 'not_verified')
        with self.assertRaisesRegex(ValueError, 'exported HOG results require'):
            hog_import.node_scope(self.results, 'N1', table, required=True)

    def test_nested_export_resolution_does_not_replace_node(self):
        wrapper = self.root/'wrapper'
        wrapper.mkdir()
        self.results.rename(wrapper/'Phylogenetic_Hierarchical_Orthogroups')
        resolved = hog_import.resolve_results(str(wrapper))
        self.assertEqual(resolved, wrapper/'Phylogenetic_Hierarchical_Orthogroups')
        with self.assertRaisesRegex(FileNotFoundError, 'N0.tsv not found'):
            parse_hogs(str(resolved), 'N0')

    def test_export_cli_audit_resume_and_missing_log_guard(self):
        with patch.object(cli, 'tool_versions', return_value={}):
            selection = cli.run(self.args())
        self.assertEqual(selection['status'], 'failed')
        status = data.read_tsv(self.root/'out/method_status.tsv')[0]
        self.assertIn('--orthofinder-export', status['reason'])
        args = self.args(['--orthofinder-export', '-o', str(self.root/'good')])
        with patch.object(cli, 'tool_versions', return_value={}):
            selection = cli.run(args)
        source = json.loads((self.root/'good/orthofinder_import.json').read_text())
        self.assertEqual(source['assigned_gene_count'], 4)
        self.assertEqual(source['node_scope']['status'], 'verified')
        self.assertIn('not independently verified', source['completion_evidence'])
        status = data.read_tsv(self.root/'good/method_status.tsv')[0]
        self.assertEqual(status['status'], 'success')
        self.assertEqual(status['import_status'], 'partial')
        self.assertIn('4/5 genes, 1 unresolved', status['reason'])
        self.assertEqual(selection['orthofinder_import']['assigned_gene_count'], 4)
        table = data.read_tsv(self.root/'good/orthofinder_assignments.tsv')
        self.assertEqual(len(table), 5)
        # Audit outputs are reproducible from a validated resume, too.
        (self.root/'good/orthofinder_assignments.tsv').unlink()
        args.resume = True
        with patch.object(cli, 'tool_versions', return_value={}):
            cli.run(args)
        self.assertEqual(table, data.read_tsv(self.root/'good/orthofinder_assignments.tsv'))
        self.assertFalse((self.results/'Log.txt').exists())
        self.assertEqual(data.mapping(self.map, self.assemblies), self.assemblies)

    def test_zero_import_fails_without_fake_singleton_candidate(self):
        self.hog.write_text('HOG\tOG\tGene Tree Parent Clade\tA_AGAT\tB_AGAT\tOut\n'
                            'N1.HOG1\tOG1\tn1\tother_a\tother_b\t\n')
        with patch.object(cli, 'tool_versions', return_value={}):
            selected = cli.run(self.args(['--orthofinder-export']))
        self.assertEqual(selected['status'], 'failed')
        self.assertEqual(data.read_tsv(self.root/'out/scores.tsv'), [])
        source = json.loads((self.root/'out/orthofinder_import.json').read_text())
        self.assertEqual(source['assigned_gene_count'], 0)
        self.assertEqual(source['unresolved_gene_count'], 5)
        self.assertEqual(source['import_status'], 'failed')
        self.assertEqual(len(data.read_tsv(self.root/'out/orthofinder_assignments.tsv')), 5)
        status = data.read_tsv(self.root/'out/method_status.tsv')[0]
        self.assertEqual(status['status'], 'failed')
        self.assertIn('zero target genes', status['reason'])

    def test_incomplete_log_not_overridden_by_export_flag(self):
        (self.results/'Log.txt').write_text('run started but still processing')
        with patch.object(cli, 'tool_versions', return_value={}):
            selected = cli.run(self.args(['--orthofinder-export']))
        self.assertEqual(selected['status'], 'failed')
        self.assertIn('does not document completed', data.read_tsv(self.root/'out/method_status.tsv')[0]['reason'])

    def test_hybrid_keeps_hog_boundary_and_unresolved_target(self):
        args = self.args(['--orthofinder-export'])
        context = dict(args=args, assemblies=self.assemblies, seqs=data.fasta(self.fa), manifest={},
                       work=self.root, runner=methods.Runner({'commands': []}, self.root, cli.DEFAULTS))
        with patch.object(methods, 'sequence_groups', side_effect=lambda method, genes, params, context, work: ([genes], [])):
            labels, unsupported, factors = methods.execute('orthofinder-mmseqs', {'identity': .8, 'coverage': .8}, context, self.root)
        self.assertEqual(labels['a.1'], labels['b'])
        self.assertNotEqual(labels['a.1'], labels['c'])
        self.assertNotEqual(labels['c'], labels['missing'])
        self.assertEqual(unsupported, ['missing'])

    def test_duplicate_normalized_source_columns_rejected(self):
        self.hog.write_text('HOG\tOG\tGene Tree Parent Clade\tA.fa\tA.pep\n'
                            'N1.HOG1\tOG1\tn1\ta\tb\n')
        with self.assertRaisesRegex(ValueError, 'duplicate assembly columns'):
            parse_hogs(str(self.results), 'N1')


if __name__ == '__main__':
    unittest.main()
