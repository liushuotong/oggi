"""Differential tests: legacy implementation is the independent reference."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import collinearity as c
import sub_collinearity as sub


class NativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if c.resolve_backend() != "rust":
            raise unittest.SkipTest("Build subcoli_rs first to run native differential tests")

    def compare(self, points, options=()):
        expected = c.collinearity(options, points.copy()).run()
        actual = c.run_blocks(options, points.copy(), "rust")
        self.assertEqual(len(actual), len(expected))
        for (a, ap, asc), (b, bp, bsc) in zip(actual, expected):
            pd.testing.assert_frame_equal(a, b)
            self.assertEqual(ap, bp)
            self.assertEqual(asc, bsc)

    def test_randomized_order_ties_overlap_and_parameters(self):
        rng = np.random.default_rng(210921)
        for seed in range(160):
            n = int(rng.integers(1, 95))
            points = pd.DataFrame(dict(loc1=rng.integers(0, 35, n),
                loc2=rng.integers(0, 35, n), grading=rng.choice([25, 40, 50], n)))
            if seed % 2:
                points = points.sort_values(["loc1", "loc2"])
            if seed % 4 == 1:
                # pandas 3 rejects fractional updates to the reference's int column.
                points["grading"] = points["grading"].astype(float)
            points.index = ["anchor_%d" % i for i in range(n)]
            options = [("gap_penalty", [-1, -.5, 0, -3][seed % 4]),
                ("mg", ["40,40", "5,9", "1,1", "15,7"][seed % 4]),
                ("over_gap", 1 + seed % 5), ("coverage_ratio", [0, .5, .8, 1][seed % 4]),
                ("pvalue", [.05, .2, 1][seed % 3])]
            with self.subTest(seed=seed):
                self.compare(points, options)

    def test_empty_and_long_plus_minus(self):
        self.compare(pd.DataFrame(columns=["loc1", "loc2", "grading"]))
        for reverse in (False, True):
            n = 100
            self.compare(pd.DataFrame(dict(loc1=np.arange(n),
                loc2=np.arange(n)[::-1] if reverse else np.arange(n), grading=[50]*n)))

    def test_batch_outputs_and_pair_limit(self):
        rng = np.random.default_rng(9)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows, members, hits = [], {}, []
            for asm in ("A", "B", "C"):
                bed = root / (asm + ".bed")
                bed.write_text("".join(f"chr{i//20}\t{i*10}\t{i*10+9}\t{asm}{i}\n" for i in range(40)))
                rows.append(dict(assembly=asm, bed=str(bed)))
                members.update({asm+str(i): asm for i in range(0, 40, 5)})
            for a, b in (("A", "B"), ("A", "C"), ("B", "C")):
                for i in range(40):
                    targets = [i, 39-i, int(rng.integers(40))]
                    for j in targets:
                        hits.append((a+str(i), b+str(j), 1e-30, 100))
                        hits.append((b+str(j), a+str(i), 1e-30, 80))
            blast = pd.DataFrame(hits, columns=["qseqid", "sseqid", "evalue", "bitscore"])
            for limit in (None, 0, 3):
                results = []
                for backend in ("python", "rust"):
                    pairs, blocks = root/(backend+".pairs"), root/(backend+".blocks")
                    df = sub.batch_member_pair_collinearity(rows, members, blast,
                        backend=backend, verbose=False, max_pairs=limit,
                        pairs_out=pairs, blocks_out=blocks)
                    results.append((df, pairs.read_bytes(), blocks.read_bytes()))
                pd.testing.assert_frame_equal(results[0][0], results[1][0])
                self.assertEqual(results[0][1:], results[1][1:])

    def test_pairwise_public_api(self):
        w1 = pd.DataFrame(dict(gene_id=['a'+str(i) for i in range(12)],loc=range(12)))
        w2 = pd.DataFrame(dict(gene_id=['b'+str(i) for i in range(12)],loc=range(12)))
        hits = pd.DataFrame(dict(qseqid=w1.gene_id,sseqid=w2.gene_id[::-1].to_numpy(),
                                 evalue=[1e-30]*12,bitscore=[100]*12))
        expected = sub.pairwise_comparison(w1,w2,hits,backend='python')
        actual = sub.pairwise_comparison(w1,w2,hits,backend='rust')
        for a,b in zip(expected,actual):
            pd.testing.assert_frame_equal(a,b)


class SelectionTests(unittest.TestCase):
    def test_empty_blast_is_valid_no_hits(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'empty.blastp'
            path.touch()
            result = sub.batch_member_pair_collinearity([],{},str(path),backend='python',verbose=False)
            self.assertEqual(len(result),0)

    def test_no_implicit_build_and_explicit_failure(self):
        with patch.object(c, "_native_library", return_value=None):
            self.assertEqual(c.resolve_backend(), "python")
            with self.assertRaises(RuntimeError):
                c.resolve_backend("rust")
        with patch.dict(os.environ, {"OGGI_SUBCOLI_LIB": "nonexistent-subcoli-library"}):
            self.assertEqual(c.resolve_backend("python"), "python")
            with self.assertRaises(FileNotFoundError):
                c.resolve_backend("rust")


if __name__ == "__main__":
    unittest.main()
