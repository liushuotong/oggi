import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_engine.metric_reports import percent, independent_metrics
from cluster_engine.scoring import quality


class MetricReportsTests(unittest.TestCase):
    def test_fixed_transforms_and_missing(self):
        self.assertEqual(percent('silhouette', -1), 0)
        self.assertEqual(percent('silhouette', 0), 50)
        self.assertEqual(percent('silhouette', 1), 100)
        self.assertEqual(percent('dunn', 0), 0)
        self.assertEqual(percent('dunn', float('inf')), 100)
        self.assertEqual(percent('davies_bouldin', 0), 100)
        self.assertIsNone(percent('dbcv', None))
        self.assertLess(percent('calinski_harabasz', 10), percent('calinski_harabasz', 20))

    def test_real_libraries_and_label_invariance(self):
        from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
        import networkx as nx
        from hdbscan.validity import validity_index
        genes = list('abcdefgh')
        x = np.array([[0,0],[0,.2],[.2,0],[.2,.2],[10,10],[10,10.2],[10.2,10],[10.2,10.2]])
        labels = {g: i//4 for i,g in enumerate(genes)}
        d = np.linalg.norm(x[:,None]-x[None,:], axis=2)
        q, _ = quality(labels, genes, d)
        edges = [(a,b,1.) for i,a in enumerate(genes) for b in genes[i+1:] if labels[a] == labels[b]]
        inputs = dict(genes=genes, features=x, edges=edges)
        rows = {r['metric']: r for r in independent_metrics(labels, genes, q, inputs)}
        self.assertTrue(all(r['status'] == 'evaluated' for r in rows.values()), rows)
        y = np.array(list(labels.values()))
        self.assertAlmostEqual(rows['silhouette']['raw_score'], silhouette_score(d, y, metric='precomputed'))
        self.assertAlmostEqual(rows['davies_bouldin']['raw_score'], davies_bouldin_score(x,y))
        self.assertAlmostEqual(rows['calinski_harabasz']['raw_score'], calinski_harabasz_score(x,y))
        self.assertAlmostEqual(rows['dbcv']['raw_score'], validity_index(x,y))
        self.assertAlmostEqual(rows['modularity']['raw_score'], .5)
        renamed = {g: 'cluster'+str(v) for g,v in labels.items()}
        self.assertEqual(list(rows.values()), independent_metrics(renamed, genes, q, inputs))
        unavailable = independent_metrics(labels, genes, q)
        self.assertTrue(all(r['score_100'] is None for r in unavailable[2:]))

    def test_cross_cluster_zero_distance(self):
        genes = list('abcd')
        labels = dict(zip(genes, [0,0,1,1]))
        d = np.array([[0,1,0,1],[1,0,1,1],[0,1,0,1],[1,1,1,0]], dtype=float)
        q,_ = quality(labels, genes, d)
        rows = independent_metrics(labels, genes, q)
        self.assertEqual(rows[1]['raw_score'], 0)
        self.assertEqual(rows[1]['score_100'], 0)


if __name__ == '__main__':
    unittest.main()
