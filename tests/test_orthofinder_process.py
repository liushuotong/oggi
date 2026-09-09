"""Wrapper contract tests; no external OrthoFinder installation required."""
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import orthofinder_process as of


class OrthoFinderWrapperTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.input = self.root / "proteomes"
        self.input.mkdir()
        for i in range(6):
            (self.input / ("sample%d.pep" % i)).write_text(">gene%d\nMPEPTIDE\n" % i)

    def test_full_command_and_current_run_resolution(self):
        output = self.root / "new_run"
        # An unrelated old result must not be selected.
        (self.input / "OrthoFinder" / "Results_old" / "Orthogroups").mkdir(parents=True)

        def fake_run(cmd, check):
            self.assertTrue(check)
            self.assertNotIn("-og", cmd)
            for flag, value in [("-I", "1.2"), ("-M", "msa"),
                                ("-A", "famsa"), ("-T", "fasttree"),
                                ("-S", "diamond")]:
                self.assertEqual(cmd[cmd.index(flag) + 1], value)
            (output / "Results_trial" / "Phylogenetic_Hierarchical_Orthogroups").mkdir(parents=True)

        with patch.object(of.shutil, "which", return_value="orthofinder"), \
                patch.object(of.subprocess, "run", side_effect=fake_run):
            result = of.run_orthofinder(str(self.input), output_dir=str(output), verbose=False)
        self.assertEqual(pathlib.Path(result), output / "Results_trial")

    def test_existing_output_rejected_before_launch(self):
        with patch.object(of.subprocess, "run") as run:
            with self.assertRaises(FileExistsError):
                of.run_orthofinder(str(self.input), output_dir=str(self.input))
            run.assert_not_called()

    def test_legacy_stop_rejected_before_launch(self):
        with patch.object(of.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                of.run_orthofinder(str(self.input), step="og")
            run.assert_not_called()

    def test_ambiguous_results_rejected(self):
        for name in ["Results_one", "Results_two"]:
            (self.root / name / "Orthogroups").mkdir(parents=True)
        with self.assertRaises(ValueError):
            of.find_results_dir(str(self.root))

    def test_string_ids_metadata_and_missing_root(self):
        directory = self.root / "Phylogenetic_Hierarchical_Orthogroups"
        directory.mkdir()
        (directory / "N2.tsv").write_text(
            "HOG\tOG\tGene Tree Parent Clade\ta\tb\n"
            "N2.HOG1\tOG1\tn8\t0001\tNA\n")
        (directory / "N2.long.tsv").write_text("ignore converted files")
        self.assertEqual(of.available_hog_levels(str(self.root)), ["N2"])
        with self.assertRaisesRegex(FileNotFoundError, "available HOG levels: N2"):
            of.parse_hogs(str(self.root))
        output = self.root / "parsed"
        table = of.parse_families(str(self.root), "N2", str(output))
        self.assertEqual(set(table.gene_ID), {"0001", "NA"})
        detail = of.parse_hogs(str(self.root), "N2", detailed=True)
        self.assertEqual(set(detail.Orthogroup), {"OG1"})
        self.assertTrue((output / "N2.metadata.long.tsv").exists())
        with (directory / "N2.tsv").open("a") as handle:
            handle.write("N2.HOG2\tOG1\tn9\t0001\t\n")
        with self.assertRaisesRegex(ValueError, "more than once"):
            of.parse_hogs(str(self.root), "N2")

    def test_parse_selected_hog_without_rerun(self):
        directory = self.root / "Phylogenetic_Hierarchical_Orthogroups"
        directory.mkdir()
        (directory / "N1.tsv").write_text(
            "HOG\tOG\tGene Tree Parent Clade\t01.col_AGAT\t02.tibet_AGAT\n"
            "N1.HOG0001\tOG0001\tn0\tcol_g1, col_g2\ttibet_g1\n"
            "N1.HOG0002\tOG0002\tn1\tcol_g3\t\n")
        with patch.object(of.subprocess, "run") as run:
            table = of.parse_families(str(self.root), level="N1")
            run.assert_not_called()
        self.assertEqual(len(table), 4)
        self.assertEqual(set(table.assembly_ID), {"01.col_AGAT", "02.tibet_AGAT"})
        self.assertEqual(set(table.gene_ID), {"col_g1", "col_g2", "col_g3", "tibet_g1"})
        self.assertTrue((directory / "N1.long.tsv").is_file())
        with self.assertRaises(ValueError):
            of.parse_families(str(self.root), "../N1")


if __name__ == "__main__":
    unittest.main()
