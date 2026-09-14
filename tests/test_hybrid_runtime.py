"""Within-HOG resource limits, checkpoint integrity and interruption recovery."""
import argparse
import contextlib
import io
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import time
import unittest
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data, methods
from cluster_engine.hybrid import refine_hogs


class HybridRuntimeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        self.groups = [['a', 'b'], ['c', 'd'], ['e', 'f']]
        self.seqs = {g: 'MPEPTIDE' * 3 for group in self.groups for g in group}
        self.args = SimpleNamespace(threads=32, hybrid_threads=1, hybrid_timeout=300)
        self.context = dict(args=self.args, seqs=self.seqs, manifest={}, fingerprint='same-run',
                            candidate_directory=self.root/'candidate')
        self.params = dict(identity=.8, coverage=.8)

    def run_hybrid(self, callback):
        with contextlib.redirect_stdout(io.StringIO()):
            return refine_hogs('orthofinder-mmseqs', self.groups, [], self.params,
                               self.context, self.root/'work', callback)

    def test_per_hog_thread_env_limit_leaves_standalone_unchanged(self):
        calls = []
        class Recorder:
            def run(inner, cmd, work, **options):
                calls.append((cmd, options))
                genes = list(data.fasta(work/'target.fa'))
                (work/'cluster_cluster.tsv').write_text(''.join(g+'\t'+g+'\n' for g in genes))
        self.context['runner'] = Recorder()
        result, unsupported = self.run_hybrid(methods.sequence_groups)
        self.assertEqual(len(result), 6)
        self.assertEqual(unsupported, [])
        for cmd, options in calls:
            self.assertEqual(cmd[cmd.index('--threads')+1], '1')
            self.assertEqual(options, dict(timeout=300, env={'MMSEQS_NUM_THREADS': '1'}))
        self.assertEqual(self.args.threads, 32)
        methods.sequence_groups('mmseqs', list(self.seqs), self.params, self.context, self.root/'standalone')
        cmd, options = calls[-1]
        self.assertEqual(cmd[cmd.index('--threads')+1], '32')
        self.assertEqual(options, {})

    def test_interruption_resumes_only_remaining_hogs(self):
        completed = []
        def interrupt(method, genes, params, context, work):
            if genes == ['e', 'f']:
                raise KeyboardInterrupt
            completed.append(genes)
            return [genes], []
        with self.assertRaises(KeyboardInterrupt):
            self.run_hybrid(interrupt)
        progress = json.loads((self.root/'candidate/hybrid_progress.json').read_text())
        self.assertEqual(progress['state'], 'interrupted')
        self.assertEqual(progress['completed_groups'], 2)
        self.assertEqual(len(list((self.root/'candidate/hog_checkpoints').glob('*.json'))), 2)
        called = []
        def finish(method, genes, params, context, work):
            called.append(genes)
            return [genes], []
        result, _ = self.run_hybrid(finish)
        self.assertEqual(called, [['e', 'f']])
        self.assertEqual(result, self.groups)
        progress = json.loads((self.root/'candidate/hybrid_progress.json').read_text())
        self.assertEqual((progress['state'], progress['completed_groups'], progress['cached_groups']), ('complete', 3, 2))

    def test_corrupt_checkpoint_rejected_and_failed_hog_not_cached(self):
        self.run_hybrid(lambda m, g, p, c, w: ([g], []))
        checkpoint = next((self.root/'candidate/hog_checkpoints').glob('*.json'))
        record = json.loads(checkpoint.read_text())
        record['result_sha256'] = 'invalid'
        checkpoint.write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            self.run_hybrid(lambda *a: self.fail('must not call a tool for a valid cached HOG'))
        other = dict(self.context, candidate_directory=self.root/'bad')
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, 'invalid partition'):
            refine_hogs('orthofinder-mmseqs', self.groups, [], self.params, other, self.root/'bad_work',
                        lambda m, g, p, c, w: ([['foreign_gene']], []))
        self.assertFalse(list((self.root/'bad/hog_checkpoints').glob('*.json')))

    def test_changed_sequence_changes_checkpoint_identity(self):
        self.run_hybrid(lambda m, g, p, c, w: ([g], []))
        self.context['seqs']['a'] += 'K'
        calls = []
        def record(m, g, p, c, w):
            calls.append(g)
            return [g], []
        self.run_hybrid(record)
        self.assertEqual(calls, [['a', 'b']])

    def test_runner_real_timeout_and_child_env_override(self):
        manifest = {'commands': []}
        runner = methods.Runner(manifest, self.root, dict(cli.DEFAULTS))
        runner.heartbeat_seconds = .1
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(TimeoutError, 'timed out'):
            runner.run([sys.executable, '-c', 'import time; time.sleep(10)'], self.root, timeout=.4)
        self.assertIsNotNone(manifest['commands'][0]['returncode'])
        self.assertEqual(manifest['commands'][0]['status'], 'failed')
        with patch.dict(os.environ, {'MMSEQS_NUM_THREADS': '32'}):
            runner.run([sys.executable, '-c', 'import os; print(os.environ["MMSEQS_NUM_THREADS"])'],
                       self.root, env={'MMSEQS_NUM_THREADS': '1'})
            self.assertEqual(os.environ['MMSEQS_NUM_THREADS'], '32')
        entry = manifest['commands'][-1]
        self.assertEqual(pathlib.Path(entry['log']).read_text().strip(), '1')
        self.assertEqual(entry['environment_overrides'], {'MMSEQS_NUM_THREADS': '1'})

    def test_runner_candidate_budget_spans_multiple_commands(self):
        manifest = {'commands': []}
        runner = methods.Runner(manifest, self.root, dict(cli.DEFAULTS, candidate_timeout=.8))
        runner.begin_candidate()
        runner.run([sys.executable, '-c', 'import time; time.sleep(.3)'], self.root)
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(TimeoutError):
            runner.run([sys.executable, '-c', 'import time; time.sleep(10)'], self.root)
        self.assertLess(manifest['commands'][1]['timeout_seconds'], manifest['commands'][0]['timeout_seconds']-.2)
        runner.end_candidate()
        self.assertIsNone(runner.candidate_deadline)

    def test_posix_cleanup_signals_only_own_process_group(self):
        process = MagicMock(pid=987654)
        process.wait.side_effect = [subprocess.TimeoutExpired('child', 2), 0]
        with patch.object(methods.os, 'name', 'posix'), patch.object(methods.os, 'killpg', create=True) as kill, \
             patch.object(methods.signal, 'SIGKILL', 9, create=True):
            methods.Runner._stop(process)
        self.assertEqual(kill.call_args_list[0].args, (987654, methods.signal.SIGTERM))
        self.assertEqual(kill.call_args_list[1].args, (987654, 9))

    def test_scheduler_records_interrupt_and_releases_lock(self):
        fa, mapping = self.root/'family.fa', self.root/'map.tsv'
        data.write_fasta(fa, self.seqs)
        data.tsv(mapping, [dict(gene_ID=g, assembly_ID='A') for g in self.seqs], ['gene_ID', 'assembly_ID'])
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        args = parser.parse_args(['cluster', '-M', 'mmseqs', '-i', str(fa), '--gene-map', str(mapping),
                                  '-o', str(self.root/'out')])
        with patch.object(cli, 'tool_versions', return_value={}), patch.object(cli, 'applicable', return_value=None), \
             patch.object(cli, 'execute', side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
            cli.run(args)
        self.assertFalse((self.root/'out/.cluster.lock').exists())
        manifest = json.loads((self.root/'out/manifest.json').read_text())
        self.assertEqual(manifest['state'], 'interrupted')
        self.assertEqual(manifest['candidate_status'][0]['status'], 'interrupted')
        self.assertEqual(data.read_tsv(self.root/'out/method_status.tsv')[0]['status'], 'interrupted')
        args.resume = True
        with patch.object(cli, 'tool_versions', return_value={}), patch.object(cli, 'applicable', return_value=None), \
             patch.object(cli, 'execute', return_value=(data.partition(self.groups, self.seqs), [], [])):
            cli.run(args)
        self.assertEqual(json.loads((self.root/'out/manifest.json').read_text())['state'], 'complete')

    def test_hybrid_resource_options_validated(self):
        parser = argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        for flag, value in [('--hybrid-threads', '0'), ('--hybrid-timeout', '0'), ('--hybrid-timeout', 'nan')]:
            args = parser.parse_args(['cluster', '-i', 'family.fa', '--gene-map', 'map.tsv', '-o', 'out', flag, value])
            with self.assertRaises(ValueError):
                cli.configure(args)


if __name__ == '__main__':
    unittest.main()
