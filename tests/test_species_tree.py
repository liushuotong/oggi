"""BUSCO species-tree contracts; external phylogeny programs are simulated."""
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import species_tree as st


class SpeciesTreeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name).resolve()
        self.busco = self.root / "busco"
        self.output = self.root / "species_tree"
        self.lineage = "viridiplantae_odb12.2"

    def marker(self, assembly, marker, sequence="MPEPTIDE", lineage=None):
        directory = (self.busco / assembly / ("run_" + (lineage or self.lineage)) /
                     "busco_sequences" / "single_copy_busco_sequences")
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / (marker + ".faa")
        path.write_text(
            ">original_{}_{} BUSCO prediction\n{}\n".format(assembly, marker, sequence),
            encoding="utf-8")
        return path

    def args(self, *extra):
        return st.build_parser().parse_args([
            "--busco-dir", str(self.busco), "-o", str(self.output),
            "--stop-after", "extract", *extra])

    def common_markers(self):
        for assembly in ("A", "B", "C", "D"):
            self.marker(assembly, "100at33090")

    def external_tools(self, commands):
        """MAFFT deliberately changes row order; trimAl removes one column."""
        def fake_run(command, **kwargs):
            command = [str(arg) for arg in command]
            commands.append(command)
            program = pathlib.Path(command[0]).stem.lower()
            if "mafft" in program:
                sequences = st.read_fasta(pathlib.Path(command[-1]))
                length = max(map(len, sequences.values()))
                text = "".join(
                    ">{}\n{}\n".format(name, sequence.ljust(length, "-"))
                    for name, sequence in reversed(list(sequences.items())))
                target = kwargs.get("stdout")
                if hasattr(target, "write"):
                    try:
                        target.write(text)
                    except TypeError:
                        target.write(text.encode("utf-8"))
                    target.flush()
                return subprocess.CompletedProcess(command, 0, stdout=text, stderr="")
            if "trimal" in program:
                source = pathlib.Path(command[command.index("-in") + 1])
                target = pathlib.Path(command[command.index("-out") + 1])
                sequences = st.read_fasta(source)
                target.write_text("".join(
                    ">{}\n{}\n".format(name, sequence[:-1])
                    for name, sequence in sequences.items()), encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
            if "iqtree" in program:
                prefix_flag = next(flag for flag in ("--prefix", "-pre", "-prefix")
                                   if flag in command)
                prefix = pathlib.Path(command[command.index(prefix_flag) + 1])
                prefix.parent.mkdir(parents=True, exist_ok=True)
                pathlib.Path(str(prefix) + ".treefile").write_text(
                    "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);\n", encoding="utf-8")
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
            self.fail("Unexpected external command: {}".format(command))
        return fake_run

    def test_extract_uses_shared_busco_ids_and_assembly_headers(self):
        self.common_markers()
        self.marker("A", "200at33090")
        self.marker("B", "300at33090")
        with patch.object(st.subprocess, "run") as run:
            st.run_species_tree(self.args())
        run.assert_not_called()
        loci = self.output / "01_loci"
        self.assertEqual([path.name for path in loci.glob("*.faa")], ["100at33090.faa"])
        self.assertEqual(st.read_fasta(loci / "100at33090.faa"), {
            "A": "MPEPTIDE", "B": "MPEPTIDE", "C": "MPEPTIDE", "D": "MPEPTIDE"})
        self.assertEqual((self.output / "selected_buscos.txt").read_text().splitlines(),
                         ["100at33090"])
        for filename in ("assemblies.tsv", "marker_occupancy.tsv", "sequence_map.tsv",
                         "run_summary.json"):
            self.assertTrue((self.output / filename).is_file(), filename)
        self.assertIn("original_A_100at33090", (self.output / "sequence_map.tsv").read_text())

    def test_no_shared_markers_fails_before_external_programs(self):
        for index, assembly in enumerate(("A", "B", "C", "D")):
            self.marker(assembly, str(index) + "at33090")
        with patch.object(st.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                st.run_species_tree(self.args())
        run.assert_not_called()

    def test_selected_busco_file_with_two_sequences_is_rejected(self):
        self.common_markers()
        path = self.marker("A", "100at33090")
        with path.open("a", encoding="utf-8") as handle:
            handle.write(">second_copy\nMPEPTIDE\n")
        with self.assertRaises(ValueError):
            st.run_species_tree(self.args())

    def test_lineages_cannot_be_mixed(self):
        self.marker("A", "100at33090")
        for assembly in ("B", "C", "D"):
            self.marker(assembly, "100at33090", lineage="embryophyta_odb10")
        with self.assertRaises(ValueError):
            st.run_species_tree(self.args())

    def test_ambiguous_lineage_requires_selection(self):
        self.common_markers()
        for assembly in ("A", "B", "C", "D"):
            self.marker(assembly, "200at33090", lineage="embryophyta_odb10")
        with self.assertRaises(ValueError):
            st.run_species_tree(self.args())

    def test_lineage_can_be_selected_with_or_without_run_prefix(self):
        self.common_markers()
        for assembly in ("A", "B", "C", "D"):
            self.marker(assembly, "200at33090", lineage="embryophyta_odb10")
        for index, lineage in enumerate((self.lineage, "run_" + self.lineage)):
            with self.subTest(lineage=lineage):
                output = self.root / ("chosen_" + str(index))
                st.run_species_tree(self.args("--lineage", lineage, "-o", str(output)))
                self.assertEqual((output / "selected_buscos.txt").read_text().splitlines(),
                                 ["100at33090"])

    def test_manifest_paths_resolve_relative_to_manifest(self):
        self.common_markers()
        manifest = self.root / "assemblies.tsv"
        manifest.write_text("assembly\tbusco_dir\n" + "".join(
            "{}_ref\tbusco/{}\n".format(assembly, assembly)
            for assembly in ("A", "B", "C", "D")), encoding="utf-8")
        args = st.build_parser().parse_args([
            "--manifest", str(manifest), "-o", str(self.output), "--stop-after", "extract"])
        st.run_species_tree(args)
        self.assertEqual(set(st.read_fasta(self.output / "01_loci" / "100at33090.faa")),
                         {"A_ref", "B_ref", "C_ref", "D_ref"})

    def test_missing_manifest_target_is_not_silently_omitted(self):
        self.common_markers()
        manifest = self.root / "assemblies.tsv"
        manifest.write_text("assembly\tbusco_dir\n" + "".join(
            "{}\tbusco/{}\n".format(assembly, assembly)
            for assembly in ("A", "B", "C", "D")) + "missing\tbusco/missing\n",
            encoding="utf-8")
        args = st.build_parser().parse_args([
            "--manifest", str(manifest), "-o", str(self.output), "--stop-after", "extract"])
        with self.assertRaises(FileNotFoundError):
            st.run_species_tree(args)

    def test_nonempty_output_is_preserved_and_rejected(self):
        self.common_markers()
        self.output.mkdir()
        old = self.output / "previous_result.txt"
        old.write_text("keep this result", encoding="utf-8")
        with patch.object(st.subprocess, "run") as run:
            with self.assertRaises(FileExistsError):
                st.run_species_tree(self.args())
        run.assert_not_called()
        self.assertEqual(old.read_text(), "keep this result")

    def test_concatenation_matches_taxon_names_and_pads_absent_loci(self):
        first = self.root / "100at33090.faa"
        second = self.root / "200at33090.faa"
        first.write_text(">B\nCC\n>A\nAA\n>C\nDD\n", encoding="utf-8")
        second.write_text(">C\nEEE\n>A\nFFF\n", encoding="utf-8")
        matrix, partitions = st.concatenate_alignments(
            [first, second], ["A", "B", "C"], self.output)
        self.assertEqual(st.read_fasta(matrix), {"A": "AAFFF", "B": "CC---", "C": "DDEEE"})
        text = partitions.read_text()
        self.assertRegex(text, r"=\s*1-2\s*;")
        self.assertRegex(text, r"=\s*3-5\s*;")

    def test_invalid_alignments_and_duplicate_fasta_ids_are_rejected(self):
        for index, text in enumerate((">A\nAA\n>B\nC\n", ">A\nAA\n>unknown\nCC\n",
                                      ">A\n\n>B\n\n", ">A\nAA\n>A duplicate\nCC\n")):
            with self.subTest(alignment=text):
                path = self.root / ("bad_{}.faa".format(index))
                path.write_text(text, encoding="utf-8")
                with self.assertRaises(ValueError):
                    st.concatenate_alignments([path], ["A", "B"],
                                              self.root / ("bad_output_{}".format(index)))
        duplicate = self.root / "duplicate.faa"
        duplicate.write_text(">gene description\nAA\n>gene other\nCC\n", encoding="utf-8")
        with self.assertRaises(ValueError):
            st.read_fasta(duplicate)

    def test_relaxed_occupancy_adds_gap_block_for_missing_assembly(self):
        self.common_markers()
        for assembly in ("A", "B", "C"):
            self.marker(assembly, "200at33090", "MGK")
        commands = []
        with patch.object(st.shutil, "which", side_effect=lambda program: program), \
                patch.object(st.subprocess, "run", side_effect=self.external_tools(commands)):
            st.run_species_tree(self.args("--min-occupancy", "0.75", "--stop-after", "concat"))
        matrix = st.read_fasta(self.output / "04_supermatrix" / "supermatrix.faa")
        self.assertEqual(matrix, {"A": "MPEPTIDEMGK", "B": "MPEPTIDEMGK",
                                  "C": "MPEPTIDEMGK", "D": "MPEPTIDE---"})
        self.assertEqual(sum("mafft" in pathlib.Path(cmd[0]).stem for cmd in commands), 2)
        self.assertFalse(any("iqtree" in pathlib.Path(cmd[0]).stem for cmd in commands))

    def test_occupancy_boundary_keeps_seven_of_fifty_assemblies_at_point_fourteen(self):
        for index in range(50):
            assembly = "A{:02d}".format(index)
            self.marker(assembly, "100at33090")
            if index < 7:
                self.marker(assembly, "200at33090", "MGK")
        with patch.object(st.subprocess, "run") as run:
            st.run_species_tree(self.args("--min-occupancy", "0.14"))
        run.assert_not_called()
        self.assertEqual((self.output / "selected_buscos.txt").read_text().splitlines(),
                         ["100at33090", "200at33090"])
        self.assertEqual(st.read_fasta(self.output / "01_loci" / "200at33090.faa"),
                         {"A{:02d}".format(index): "MGK" for index in range(7)})

    def test_mafft_failure_records_failed_status_and_never_starts_iqtree(self):
        self.common_markers()
        commands = []

        def fail_mafft(command, **kwargs):
            commands.append(command)
            self.assertIn("mafft", pathlib.Path(command[0]).stem)
            self.assertTrue(kwargs.get("check"))
            kwargs["stderr"].write("simulated alignment failure\n")
            raise subprocess.CalledProcessError(23, command)

        with patch.object(st.shutil, "which", side_effect=lambda program: program), \
                patch.object(st.subprocess, "run", side_effect=fail_mafft):
            with self.assertRaisesRegex(RuntimeError, "exit 23"):
                st.run_species_tree(self.args("--stop-after", "tree"))
        self.assertEqual(len(commands), 1)
        self.assertFalse(any("iqtree" in pathlib.Path(command[0]).stem for command in commands))
        summary = json.loads((self.output / "run_summary.json").read_text(encoding="utf-8"))
        self.assertEqual(summary["status"], "failed")
        self.assertIn("exit 23", summary["error"])
        self.assertIn("simulated alignment failure",
                      (self.output / "logs" / "100at33090.mafft.log").read_text())
        self.assertFalse((self.output / "04_supermatrix" / "supermatrix.faa").exists())
        self.assertFalse(list(self.output.rglob("*.treefile")))

    def test_mocked_pipeline_trims_alignments_and_passes_partitioned_matrix_to_iqtree(self):
        first_sequences = {"A": "MKTA", "B": "MLTA", "C": "MVTA", "D": "MITA"}
        second_sequences = {"A": "GGK", "B": "GAK", "C": "GCK", "D": "GDK"}
        for assembly in ("A", "B", "C", "D"):
            self.marker(assembly, "100at33090", first_sequences[assembly])
            self.marker(assembly, "200at33090", second_sequences[assembly])
        commands = []
        with patch.object(st.shutil, "which", side_effect=lambda program: program), \
                patch.object(st.subprocess, "run", side_effect=self.external_tools(commands)):
            st.run_species_tree(self.args("--stop-after", "tree", "--trim", "automated1",
                                          "--threads", "4", "--jobs", "2"))
        matrix = self.output / "04_supermatrix" / "supermatrix.faa"
        partitions = self.output / "04_supermatrix" / "partitions.nex"
        self.assertEqual(st.read_fasta(matrix), {
            "A": "MKTGG", "B": "MLTGA", "C": "MVTGC", "D": "MITGD"})
        self.assertEqual(len(list((self.output / "03_trimmed").glob("*.faa"))), 2)
        iqtree_commands = [cmd for cmd in commands if "iqtree" in pathlib.Path(cmd[0]).stem]
        self.assertEqual(len(iqtree_commands), 1)
        command = iqtree_commands[0]
        for flag, value in (("-s", str(matrix)), ("-p", str(partitions)), ("-st", "AA"),
                            ("-m", "MFP"), ("-B", "1000"), ("-alrt", "1000")):
            self.assertIn(flag, command)
            self.assertEqual(command[command.index(flag) + 1], value)
        self.assertEqual(len(list(self.output.rglob("*.treefile"))), 1)


if __name__ == "__main__":
    unittest.main()
