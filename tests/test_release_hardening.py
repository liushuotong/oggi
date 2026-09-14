import importlib.util
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from cluster_engine import data, scoring, report
from run_tests import exit_code


class ReleaseHardeningTests(unittest.TestCase):
    def test_strict_runner_rejects_empty_skipped_and_failed(self):
        result = unittest.TestResult()
        self.assertEqual(exit_code(result), 1)
        result.testsRun = 1
        self.assertEqual(exit_code(result), 0)
        result.skipped.append(('test', 'missing dependency'))
        self.assertEqual(exit_code(result), 1)
        result.skipped.clear()
        result.failures.append(('test', 'failure'))
        self.assertEqual(exit_code(result), 1)

    def test_numeric_schema_preserves_identifiers_and_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp)/'data.tsv'
            data.tsv(path, [{'gene_ID':'NA','Q':None}, {'gene_ID':'0001','Q':0}], ['gene_ID','Q'])
            rows = data.read_tsv(path, numeric_fields=['Q'])
            self.assertEqual(rows, [{'gene_ID':'NA','Q':None}, {'gene_ID':'0001','Q':0.0}])
            with self.assertRaises(ValueError):
                data.tsv(path, [{'gene_ID':None}], ['gene_ID'])
            path.write_text('Q\nnan\n')
            with self.assertRaises(ValueError):
                data.read_tsv(path, numeric_fields=['Q'])

    def test_boolean_is_not_a_score(self):
        with self.assertRaisesRegex(ValueError, 'invalid score metric B'):
            scoring.totals([{'B':True}], scoring.WEIGHTS)

    def test_missing_distance_chart_explains_na(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp)/'chart.svg'
            report.chart(path, 'diagnostics', [('similarity', [])], {0:'No trusted alignment'})
            self.assertIn('No trusted alignment', path.read_text())
            self.assertIn('NA / not evaluable', path.read_text())

    def test_legacy_explicit_mapping_is_required_for_every_gene(self):
        from mmseqs_process import parse_mmseqs_cluster
        from cdhit_process import process_cdhit_result
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            mm = root/'clusters.tsv'
            cd = root/'clusters.clstr'
            mm.write_text('A_g1\tA_g1\nA_g1\tB_g1\n')
            cd.write_text('>Cluster 0\n0 20aa, >A_g1... *\n1 20aa, >B_g1... at 95%\n')
            for parser, path in [(parse_mmseqs_cluster, mm),(process_cdhit_result, cd)]:
                for mapping in [{}, {'A_g1':'sample1'}]:
                    with self.assertRaisesRegex(ValueError, 'missing assembly mapping'):
                        parser(path, mapping)
                result = parser(path, {'A_g1':'sample1','B_g1':'sample2'})
                self.assertEqual(list(result.assembly_ID), ['sample1','sample2'])

    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython required for legacy helper')
    def test_diamond_argument_forwarding_without_shell(self):
        import sub_collinearity_pre_process as pre
        filename = 'family with spaces & symbols.pep'
        for limit in [0,7]:
            with patch.object(pre.subprocess, 'run') as run:
                pre.sub_collinearity_blastp(filename, 1000, max_target_seqs=limit)
                cmd = run.call_args_list[1].args[0]
                self.assertEqual(cmd[cmd.index('--query')+1], filename)
                self.assertEqual(cmd[cmd.index('--max-target-seqs')+1], str(limit))
                self.assertFalse(run.call_args.kwargs.get('shell', False))
        with self.assertRaises(ValueError):
            pre.sub_collinearity_blastp(filename, 1000, max_target_seqs=-1)

