"""BUSCO interface contracts; BUSCO itself is simulated in these tests."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import busco_process as bp
import species_tree as st


class BuscoProcessTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.inputs = self.root / "input sequences"
        self.inputs.mkdir()
        self.output = self.root / "BUSCO results"
        self.lineage = "viridiplantae_odb12.2"

    def fasta(self, name="assembly_A.faa", directory=None):
        target = (directory or self.inputs) / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(">protein_1\nMPEPTIDE\n", encoding="utf-8")
        return target

    @staticmethod
    def argument(command, *flags):
        for flag in flags:
            if flag in command:
                return command[command.index(flag) + 1]
        raise AssertionError("Missing argument {} in {}".format(flags, command))

    def fake_busco(self, commands, fail_on=None):
        def fake_run(command, **kwargs):
            self.assertIsInstance(command, (list, tuple))
            self.assertFalse(kwargs.get("shell", False))
            self.assertTrue(all(isinstance(argument, str) for argument in command))
            command = list(command)
            commands.append(command)
            if fail_on == len(commands):
                raise subprocess.CalledProcessError(17, command)
            output = (Path(self.argument(command, "--out_path")) /
                      self.argument(command, "-o", "--out"))
            lineage = Path(self.argument(command, "-l", "--lineage_dataset")).name
            run = output / ("run_" + lineage)
            sequences = run / "busco_sequences" / "single_copy_busco_sequences"
            sequences.mkdir(parents=True)
            (run / "full_table.tsv").write_text(
                "# BUSCO test fixture\n100at33090\tMissing\n", encoding="utf-8")
            (run / "short_summary.json").write_text(
                json.dumps({"results": {"Complete percentage": 0.0}}), encoding="utf-8")
            return subprocess.CompletedProcess(command, 0, stdout="", stderr="")
        return fake_run

    def run_batch(self, **kwargs):
        options = dict(input_path=self.inputs, output_dir=self.output,
                       lineage=self.lineage)
        options.update(kwargs)
        return bp.run_busco_batch(**options)

    def test_parser_requires_one_input_and_explicit_lineage(self):
        parser = bp.build_parser()
        args = parser.parse_args(["-i", str(self.inputs), "-o", str(self.output),
                                  "-l", self.lineage])
        self.assertEqual(args.mode, "proteins")
        self.assertEqual(args.threads, 8)
        for extra in (["-i", str(self.inputs)],
                      ["-l", self.lineage],
                      ["-i", str(self.inputs), "--manifest", "assemblies.tsv",
                       "-l", self.lineage]):
            with self.subTest(arguments=extra):
                with self.assertRaises(SystemExit), patch("sys.stderr"):
                    parser.parse_args(["-o", str(self.output), *extra])

    def test_single_run_preserves_paths_modes_and_optional_arguments(self):
        source = self.fasta()
        download = self.root / "download cache"
        commands = []
        with patch.object(bp.shutil, "which", return_value="/tools/busco"), \
                patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands)):
            result = bp.run_busco(
                source, self.output / "assembly_A", self.lineage,
                mode="genome", threads=3, busco="custom-busco",
                download_path=download, offline=True)
        self.assertEqual(len(commands), 1)
        command = commands[0]
        self.assertEqual(self.argument(command, "-i", "--in"), str(source))
        self.assertEqual(self.argument(command, "-m", "--mode"), "genome")
        self.assertEqual(self.argument(command, "-c", "--cpu"), "3")
        self.assertEqual(self.argument(command, "--download_path"), str(download))
        self.assertIn("--offline", command)
        self.assertEqual(result["assembly"], "assembly_A")
        self.assertEqual(Path(result["input_file"]), source)
        self.assertEqual(Path(result["busco_dir"]), self.output / "assembly_A")
        self.assertEqual(Path(result["run_dir"]), self.output / "assembly_A" /
                         ("run_" + self.lineage))
        self.assertEqual(result["command"], command)

    def test_local_lineage_directory_uses_dataset_basename_for_run(self):
        dataset = self.root / "local lineages" / self.lineage
        dataset.mkdir(parents=True)
        (dataset / "dataset.cfg").write_text("name=viridiplantae\n", encoding="utf-8")
        commands = []
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands)):
            result = bp.run_busco(self.fasta(), self.output / "assembly_A", dataset)
        self.assertEqual(self.argument(commands[0], "-l", "--lineage_dataset"),
                         str(dataset))
        self.assertEqual(Path(result["run_dir"]).name, "run_" + self.lineage)

    def test_unversioned_lineage_is_rejected_before_single_or_batch_run(self):
        source = self.fasta()
        original = source.read_bytes()
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                bp.run_busco(source, self.output / "assembly_A", "viridiplantae")
            self.assertFalse(self.output.exists())
            with self.assertRaises(ValueError):
                self.run_batch(lineage="viridiplantae")
        run.assert_not_called()
        self.assertFalse(self.output.exists())
        self.assertEqual(source.read_bytes(), original)

    def test_download_cache_cannot_replace_batch_metadata_or_its_children(self):
        source = self.fasta()
        original = source.read_bytes()
        for metadata in ("busco_manifest.tsv", "busco_run_summary.json"):
            for suffix in (Path(), Path("nested") / "cache"):
                with self.subTest(metadata=metadata, suffix=str(suffix)):
                    download = self.output / metadata / suffix
                    with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                            patch.object(bp.subprocess, "run") as run:
                        with self.assertRaises(ValueError):
                            self.run_batch(download_path=download)
                    run.assert_not_called()
                    self.assertFalse(self.output.exists())
                    self.assertEqual(source.read_bytes(), original)

    def test_batch_writes_absolute_manifest_accepted_by_species_tree(self):
        for name in ("assembly_B.faa", "assembly_A.fasta"):
            self.fasta(name)
        commands = []
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands)):
            summary = self.run_batch()
        self.assertEqual(len(commands), 2)
        manifest = self.output / "busco_manifest.tsv"
        with manifest.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["assembly"] for row in rows], ["assembly_A", "assembly_B"])
        self.assertTrue(all(Path(row["busco_dir"]).is_absolute() for row in rows))
        inputs, skipped = st.discover_inputs(manifest=manifest, lineage=self.lineage)
        self.assertEqual(set(inputs), {"assembly_A", "assembly_B"})
        self.assertEqual(skipped, [])
        saved = json.loads((self.output / "busco_run_summary.json").read_text())
        self.assertEqual(saved, summary)

    def test_duplicate_stems_ids_and_existing_output_fail_before_busco(self):
        self.fasta("assembly_A.fa")
        self.fasta("assembly_A.faa")
        manifest = self.root / "duplicates.tsv"
        manifest.write_text("assembly\tfasta\nA\tinput sequences/assembly_A.fa\n"
                            "A\tinput sequences/assembly_A.faa\n", encoding="utf-8")
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                self.run_batch()
            with self.assertRaises(ValueError):
                self.run_batch(input_path=None, manifest=manifest)
            existing = self.root / "existing_output"
            existing.mkdir()
            retained = existing / "previous_result.txt"
            retained.write_text("preserve this result", encoding="utf-8")
            with self.assertRaises(FileExistsError):
                bp.run_busco(self.inputs / "assembly_A.fa", existing, self.lineage)
            self.assertEqual(retained.read_text(), "preserve this result")
        run.assert_not_called()

    def test_missing_manifest_source_fails_before_any_batch_run(self):
        source = self.fasta()
        manifest = self.root / "missing.tsv"
        manifest.write_text("assembly\tpep\nA\t{}\nB\tmissing.faa\n".format(source),
                            encoding="utf-8")
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run") as run:
            with self.assertRaises(FileNotFoundError):
                self.run_batch(input_path=None, manifest=manifest)
        run.assert_not_called()

    def test_manifest_supports_mode_columns_and_relative_reduce_paths(self):
        source = self.fasta()
        commands = []
        for index, (mode, column) in enumerate((("proteins", "pep"), ("genome", "genome"),
                                                ("transcriptome", "transcriptome"),
                                                ("proteins", "fasta"))):
            with self.subTest(mode=mode, column=column):
                manifest = self.root / ("manifest_{}.tsv".format(index))
                manifest.write_text("assembly\t{}\nrenamed_A\tinput sequences/{}\n".format(
                    column, source.name), encoding="utf-8")
                with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                        patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands)):
                    self.run_batch(input_path=None, manifest=manifest, mode=mode,
                                   output_dir=self.root / ("out_{}".format(index)))
                self.assertEqual(self.argument(commands[-1], "-i", "--in"), str(source))
                self.assertEqual(self.argument(commands[-1], "-o", "--out"), "renamed_A")
                self.assertEqual(self.argument(commands[-1], "-m", "--mode"), mode)

    def test_manifest_cwd_relative_paths_work_and_ambiguity_is_rejected(self):
        source = self.fasta()
        manifest_dir = self.root / "manifests"
        manifest_dir.mkdir()
        manifest = manifest_dir / "reduced.tsv"
        manifest.write_text("assembly\tpep\nA\tinput sequences/assembly_A.faa\n",
                            encoding="utf-8")
        original_cwd = Path.cwd()
        self.addCleanup(os.chdir, original_cwd)
        os.chdir(self.root)
        commands = []
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands)):
            self.run_batch(input_path=None, manifest=manifest)
        self.assertEqual(self.argument(commands[0], "-i", "--in"), str(source))
        self.fasta(directory=manifest_dir / "input sequences")
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run") as run:
            with self.assertRaises(ValueError):
                self.run_batch(input_path=None, manifest=manifest,
                               output_dir=self.root / "ambiguous_output")
        run.assert_not_called()

    def test_failure_stops_batch_and_does_not_publish_partial_manifest(self):
        for assembly in ("A", "B", "C"):
            self.fasta(assembly + ".faa")
        commands = []
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run", side_effect=self.fake_busco(commands, fail_on=2)):
            with self.assertRaises((RuntimeError, subprocess.CalledProcessError)):
                self.run_batch()
        self.assertEqual(len(commands), 2)
        self.assertFalse((self.output / "busco_manifest.tsv").exists())
        summary = json.loads((self.output / "busco_run_summary.json").read_text())
        self.assertEqual(summary["status"], "failed")
        self.assertFalse((self.output / "C").exists())

    def test_missing_executable_and_incomplete_success_are_errors(self):
        source = self.fasta()
        with patch.object(bp.shutil, "which", return_value=None), \
                patch.object(bp.subprocess, "run") as run:
            with self.assertRaises((FileNotFoundError, RuntimeError)):
                bp.run_busco(source, self.output / "unavailable", self.lineage)
        run.assert_not_called()
        with patch.object(bp.shutil, "which", side_effect=lambda name: name), \
                patch.object(bp.subprocess, "run", return_value=subprocess.CompletedProcess(
                    ["busco"], 0, stdout="", stderr="")):
            with self.assertRaises((FileNotFoundError, RuntimeError, ValueError)):
                bp.run_busco(source, self.output / "incomplete", self.lineage)


if __name__ == "__main__":
    unittest.main()
