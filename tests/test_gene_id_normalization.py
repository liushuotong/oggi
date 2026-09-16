"""Regression coverage for automatic, assembly-wide ID normalization."""
import argparse
import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gene_id_utils import prepare_gene_inputs


def fasta_ids(path):
    with open(path) as handle:
        return [line[1:].split()[0] for line in handle if line.startswith('>')]


def bed_ids(path):
    with open(path) as handle:
        return [line.rstrip('\n').split('\t')[3] for line in handle
                if line.strip() and not line.startswith('#')]


class GeneInputFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.cache = self.root / 'normalized cache'

    def assembly(self, name, ids, with_bed=True):
        pep = self.root / (name + '.pep')
        pep.write_text(''.join('>%s original description\nMALWMRLLPLL\n' % gene
                               for gene in ids))
        row = {'assembly': name, 'pep': str(pep)}
        if with_bed:
            bed = self.root / (name + '.bed')
            bed.write_text('# coordinates\n' + ''.join(
                'chr1\t%d\t%d\t%s\t0\t+\n' % (i * 100, i * 100 + 80, gene)
                for i, gene in enumerate(ids)))
            row['bed'] = str(bed)
        return row


class GeneInputNormalizationTests(GeneInputFixture):
    def test_unique_inputs_keep_paths_bytes_and_create_no_cache(self):
        rows = [self.assembly('A', ['x', 'y']), self.assembly('B', ['z'])]
        original = {row[key]: Path(row[key]).read_bytes()
                    for row in rows for key in ('pep', 'bed')}
        prepared = prepare_gene_inputs(rows, str(self.cache))
        self.assertEqual(prepared, rows)
        self.assertFalse(self.cache.exists())
        for path, content in original.items():
            self.assertEqual(Path(path).read_bytes(), content)

    def test_any_cross_assembly_duplicate_prefixes_every_id_and_bed(self):
        rows = [self.assembly('A', ['shared', 'onlyA']),
                self.assembly('B', ['shared', 'onlyB']),
                self.assembly('Outgroup', ['outgroup_gene'])]
        original = {row[key]: Path(row[key]).read_bytes()
                    for row in rows for key in ('pep', 'bed')}
        prepared = prepare_gene_inputs(rows, str(self.cache))
        all_ids = []
        for source, normalized in zip(rows, prepared):
            expected = [source['assembly'] + '_' + gene
                        for gene in fasta_ids(source['pep'])]
            self.assertEqual(fasta_ids(normalized['pep']), expected)
            self.assertEqual(bed_ids(normalized['bed']), expected)
            self.assertNotEqual(normalized['pep'], source['pep'])
            self.assertNotEqual(normalized['bed'], source['bed'])
            self.assertTrue(Path(normalized['pep']).is_absolute())
            self.assertTrue(Path(normalized['bed']).is_absolute())
            self.assertIn('original description', Path(normalized['pep']).read_text())
            self.assertEqual(Path(normalized['pep']).read_text().count('MALWMRLLPLL'),
                             len(expected))
            all_ids.extend(expected)
        self.assertEqual(len(all_ids), len(set(all_ids)))
        for path, content in original.items():
            self.assertEqual(Path(path).read_bytes(), content)
        manifests = list(self.cache.rglob('assembly_manifest.tsv'))
        self.assertEqual(len(manifests), 1)
        with manifests[0].open() as handle:
            manifest = list(csv.DictReader(handle, delimiter='\t'))
        for actual, expected in zip(manifest, prepared):
            for key in ('assembly', 'pep', 'bed'):
                self.assertEqual(actual[key], expected[key])
        self.assertEqual(len(manifest), len(prepared))

    def test_repeated_runs_and_already_normalized_input_do_not_stack_prefixes(self):
        rows = [self.assembly('A', ['shared']), self.assembly('B', ['shared'])]
        first = prepare_gene_inputs(rows, str(self.cache))
        second = prepare_gene_inputs(rows, str(self.cache))
        self.assertEqual(first, second)
        other_cache = self.root / 'unused_cache'
        self.assertEqual(prepare_gene_inputs(first, str(other_cache)), first)
        self.assertFalse(other_cache.exists())
        self.assertEqual(fasta_ids(second[0]['pep']), ['A_shared'])
        self.assertEqual(fasta_ids(second[1]['pep']), ['B_shared'])

    def test_cache_does_not_reuse_stale_sequences(self):
        rows = [self.assembly('A', ['shared']), self.assembly('B', ['shared'])]
        first = prepare_gene_inputs(rows, str(self.cache))
        pep = Path(rows[0]['pep'])
        pep.write_text(pep.read_text().replace('MALWMRLLPLL', 'MALWMRLLPLV'))
        second = prepare_gene_inputs(rows, str(self.cache))
        self.assertIn('MALWMRLLPLV', Path(second[0]['pep']).read_text())
        self.assertNotEqual(first[0]['pep'], second[0]['pep'])

    def test_protein_only_inputs_are_supported(self):
        rows = [self.assembly('A', ['shared'], with_bed=False),
                self.assembly('B', ['shared'], with_bed=False)]
        prepared = prepare_gene_inputs(rows, str(self.cache))
        self.assertEqual(fasta_ids(prepared[0]['pep']), ['A_shared'])
        self.assertEqual(fasta_ids(prepared[1]['pep']), ['B_shared'])

    def test_duplicate_inside_one_assembly_is_an_error(self):
        rows = [self.assembly('A', ['duplicate', 'duplicate'])]
        with self.assertRaises(ValueError):
            prepare_gene_inputs(rows, str(self.cache))

    def test_duplicate_assembly_name_is_an_error(self):
        rows = [self.assembly('A', ['a']), self.assembly('B', ['b'])]
        rows[1]['assembly'] = 'A'
        with self.assertRaises(ValueError):
            prepare_gene_inputs(rows, str(self.cache))

    def test_prefix_collision_is_detected(self):
        rows = [self.assembly('A', ['shared', 'B_x']),
                self.assembly('A_B', ['shared', 'x'])]
        with self.assertRaises(ValueError):
            prepare_gene_inputs(rows, str(self.cache))

    def test_protein_missing_from_bed_is_an_error(self):
        rows = [self.assembly('A', ['shared']), self.assembly('B', ['shared'])]
        Path(rows[1]['bed']).write_text('chr1\t0\t80\tsomething_else\n')
        with self.assertRaises(ValueError):
            prepare_gene_inputs(rows, str(self.cache))


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython is required')
class IdentificationNormalizationTests(GeneInputFixture):
    def test_duplicate_outside_family_hits_still_prefixes_all_family_members(self):
        import gene_family_identification as gfi
        rows = [self.assembly('A', ['background', 'centerA']),
                self.assembly('B', ['background', 'centerB'])]

        def identify(hmms, refs, eh, eb, seq_file, output, assembly, cpu):
            self.assertEqual(fasta_ids(seq_file),
                             [assembly + '_background', assembly + '_center' + assembly])
            return {assembly + '_center' + assembly}

        with patch.object(gfi, 'identification_caculation', side_effect=identify) as search:
            mapping, ids, family = gfi.main_identification(
                [[row['assembly'], row['pep']] for row in rows], [], [],
                1e-5, 1e-5, str(self.root / 'family'), 1,
                bed_by_assembly={row['assembly']: row['bed'] for row in rows},
                input_cache_dir=str(self.cache))
        self.assertEqual(search.call_count, 2)
        self.assertEqual(mapping, {'A_centerA': 'A', 'B_centerB': 'B'})
        self.assertEqual(set(ids), set(mapping))
        self.assertEqual(set(fasta_ids(family)), set(mapping))

    def test_invalid_input_fails_before_external_search(self):
        import gene_family_identification as gfi
        row = self.assembly('A', ['duplicate', 'duplicate'])
        with patch.object(gfi, 'identification_caculation') as search:
            with self.assertRaises(ValueError):
                gfi.main_identification([['A', row['pep']]], [], [], 1e-5, 1e-5,
                                        str(self.root / 'family'), 1)
        search.assert_not_called()

    def test_unique_ids_are_unchanged_in_core_identification(self):
        import gene_family_identification as gfi
        rows = [self.assembly('A', ['originalA']), self.assembly('B', ['originalB'])]
        with patch.object(gfi, 'identification_caculation',
                          side_effect=lambda h, r, eh, eb, path, out, asm, cpu:
                          set(fasta_ids(path))):
            mapping, ids, family = gfi.main_identification(
                [[row['assembly'], row['pep']] for row in rows], [], [],
                1e-5, 1e-5, str(self.root / 'family'), 1,
                input_cache_dir=str(self.cache))
        self.assertEqual(mapping, {'originalA': 'A', 'originalB': 'B'})
        self.assertEqual(set(fasta_ids(family)), set(mapping))
        self.assertFalse(self.cache.exists())

    def test_identify_then_subcoli_accepts_original_manifest(self):
        import gene_family_identification as gfi
        import oggi
        import sub_collinearity_pre_process as pre
        rows = [self.assembly(asm, ['g%d' % i for i in range(1, 7)])
                for asm in ('A', 'B')]
        manifest = self.root / 'assembly_manifest.tsv'
        with manifest.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['assembly', 'pep', 'bed'],
                                    delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows(rows)
        family = self.root / 'family'
        identify_args = argparse.Namespace(
            manifest=str(manifest), hmm='mock.hmm', ref='mock.fa',
            output=str(family), evalue_hmm=1e-5, evalue_blastp=1e-5, threads=1)

        def search(h, r, eh, eb, path, out, assembly, cpu):
            self.assertIn(assembly + '_g3', fasta_ids(path))
            return {assembly + '_g3'}

        with patch.object(gfi, 'identification_caculation', side_effect=search):
            oggi.run_identify(identify_args)

        def window_search(**kwargs):
            genes = set()
            for asm, pep in kwargs['assembly_file_dict']:
                expected = [asm + '_g%d' % i for i in range(1, 7)]
                self.assertEqual(fasta_ids(pep), expected)
                self.assertEqual(bed_ids(kwargs['bed_of'][pep]), expected)
                genes.update(expected)
            blast = self.root / 'window.blastp'
            blast.write_text(''.join(
                'A_g%d\tB_g%d\t95\t20\t0\t0\t1\t20\t1\t20\t1e-40\t100\n' % (i, i)
                for i in range(1, 7)))
            return str(blast), genes

        subcoli_args = argparse.Namespace(
            manifest=str(manifest), id_table=str(family / 'gene_to_assembly.tsv'),
            blast=None, output=str(self.root / 'subcoli'), up=6, down=6,
            evalue=1e-5, pvalue=.2, max_pairs=None, threads=1)
        with patch.object(pre, 'seq_BLASTP_for_collinearity',
                          side_effect=window_search) as windows:
            oggi.run_subcoli(subcoli_args)
        windows.assert_called_once()
        pairs = (self.root / 'subcoli.collinear_pairs.tsv').read_text()
        self.assertIn('A_g3\tB_g3', pairs)
        self.assertEqual(fasta_ids(rows[0]['pep']), ['g%d' % i for i in range(1, 7)])
        self.assertEqual(bed_ids(rows[0]['bed']), ['g%d' % i for i in range(1, 7)])


if __name__ == '__main__':
    unittest.main()
