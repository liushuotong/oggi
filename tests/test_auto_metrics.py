"""Auto evaluates all metrics without manually supplied features or graphs."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data
from cluster_engine.auto_metrics import pcoa_features, knn_graph, sequence_distances, prepare_inputs


class AutoMetricTests(unittest.TestCase):
    def test_pcoa_preserves_euclidean_distances_and_duplicates(self):
        from scipy.spatial.distance import pdist, squareform
        x = np.array([[0.,0.],[0,0],[.1,.2],[2,2],[2.2,2]])
        d = squareform(pdist(x))
        features, info = pcoa_features(d)
        np.testing.assert_allclose(squareform(pdist(features))*info['distance_scale'], d, atol=1e-10)
        np.testing.assert_array_equal(features[0], features[1])
        self.assertLess(info['distance_stress'], 1e-10)
        # A non-Euclidean distance needs a declared, measurable approximation.
        _, distorted = pcoa_features(np.array([[0.,1,1],[1,0,3],[1,3,0]]))
        self.assertGreater(distorted['negative_inertia_fraction'], 0)
        self.assertGreater(distorted['distance_stress'], 0)

    def test_graph_and_sampling_are_order_independent(self):
        genes = list('abcd')
        d = np.array([[0.,1,1,4],[1,0,2,3],[1,2,0,3],[4,3,3,0]])
        edges,_ = knn_graph(genes,d,1)
        self.assertIn(('a','b',np.exp(-1/2.5)), edges)
        self.assertIn(('a','c',np.exp(-1/2.5)), edges)
        reverse,_ = knn_graph(genes[::-1],d[::-1,::-1],1)
        canonical = lambda es: {tuple(sorted((a,b))):w for a,b,w in es}
        self.assertEqual(canonical(edges),canonical(reverse))
        seqs = {str(i):'MACDEFGHIK'*3+'AC'*i for i in range(12)}
        a,da,info=sequence_distances(seqs,6,42)
        b,db,_=sequence_distances(dict(reversed(list(seqs.items()))),6,42)
        self.assertEqual(a,b)
        np.testing.assert_array_equal(da,db)
        self.assertTrue(info['sampled'])

    def test_auto_cli_no_extra_inputs_six_metrics_and_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            rng=np.random.default_rng(42)
            genes=list('abcdefgh')
            seqs={g:''.join(rng.choice(list('ACDEFGHIKL' if i<4 else 'MNPQRSTVWY'),100)) for i,g in enumerate(genes)}
            data.write_fasta(root/'input.fa',seqs)
            data.tsv(root/'map.tsv',[dict(gene_ID=g,assembly_ID=g) for g in genes],['gene_ID','assembly_ID'])
            p=argparse.ArgumentParser()
            cli.add_parser(p.add_subparsers())
            args=p.parse_args(['cluster','-M','auto','-i',str(root/'input.fa'),'--gene-map',str(root/'map.tsv'),'-o',str(root/'out')])
            labels={g:str(i//4) for i,g in enumerate(genes)}
            with patch.object(cli,'tool_versions',return_value={}), \
                 patch.object(cli,'applicable',side_effect=lambda m,c: None if m in ('mmseqs','cdhit') else 'test adapter unavailable'), \
                 patch.object(cli,'execute',return_value=(labels,[],[])):
                selection=cli.run(args)
                rows=data.read_tsv(root/'out'/'metric_scores.tsv')
                self.assertEqual(len(rows),12)
                self.assertTrue(all(r['status']=='evaluated' for r in rows), rows)
                self.assertTrue(all(0 <= float(r['score_100']) <= 100 for r in rows))
                self.assertEqual(selection['status'],'equivalent_best')
                self.assertTrue((root/'out'/'evaluation'/'features.tsv').exists())
                self.assertTrue((root/'out'/'evaluation'/'graph.tsv').exists())
                manifest=json.loads((root/'out'/'manifest.json').read_text())
                self.assertEqual(manifest['distance']['source'],'auto-dipeptide')
                self.assertEqual(len(manifest['scoring']['diagnostic_metrics']),6)
                args.resume=True
                cli.run(args)
                self.assertEqual(rows,data.read_tsv(root/'out'/'metric_scores.tsv'))
                # A supplied unusable gene tree must not silently switch to composition.
                (root/'bad.tree').write_text('(outsider:1,other:1);')
                args.gene_tree=str(root/'bad.tree')
                args.output=str(root/'bad_out')
                args.resume=False
                invalid=cli.run(args)
                self.assertEqual(invalid['score_status'],'not_evaluable')
                bad_rows=data.read_tsv(root/'bad_out'/'metric_scores.tsv')
                self.assertTrue(all(r['score_100']=='NA' for r in bad_rows))

    def test_auto_respects_explicit_inputs_on_one_shared_scope(self):
        with tempfile.TemporaryDirectory() as temp:
            genes=list('abcd')
            x=np.array([[0.,0],[0,1],[2,1],[2,0]])
            explicit=dict(genes=genes,features=x,edges=[('a','b',3.),('b','c',2.),('c','d',4.)])
            inputs,info=prepare_inputs(['a','b','c'],np.ones((3,3))-np.eye(3),explicit,cli.DEFAULTS,temp,'test')
            np.testing.assert_array_equal(inputs['features'],x[:3])
            self.assertEqual(inputs['edges'],explicit['edges'][:2])
            self.assertEqual(inputs['genes'],['a','b','c'])


if __name__=='__main__':
    unittest.main()
