"""Real Louvain and subprocess integration on synthetic graphs."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import networkx as nx

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from cluster_engine import cli, data, methods
from cluster_engine.louvain_worker import cluster


class LouvainTests(unittest.TestCase):
    def test_real_weighted_communities_isolates_and_order(self):
        genes=list('abcdefg')
        edges=[(a,b,1.) for group in ('abc','def') for i,a in enumerate(group) for b in group[i+1:]]
        edges.append(('c','d',.001))
        result=cluster(genes,edges,1,42)
        self.assertEqual(result['groups'],[['a','b','c'],['d','e','f'],['g']])
        self.assertEqual(result,cluster(genes[::-1],edges[::-1],1,42))
        self.assertEqual(cluster(genes,[],1,42)['groups'],[[g] for g in genes])
        with self.assertRaisesRegex(ValueError,'self-loop'):
            cluster(genes,[('a','a',1.)])

    def test_resolution_parameter_changes_granularity(self):
        genes=list('abcdef')
        edges=[(a,b,1.) for i,a in enumerate(genes) for b in genes[i+1:]]
        self.assertEqual(len(cluster(genes,edges,.5)['groups']),1)
        self.assertEqual(len(cluster(genes,edges,3.)['groups']),6)
        for invalid in (0,-1,float('inf'),float('nan')):
            with self.assertRaises(ValueError):
                cluster(genes,edges,invalid)

    def test_shared_mcl_weights_and_no_evidence_fallback(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            (root/'synteny.tsv').write_text('a\tb\n')
            args=argparse.Namespace(tree=None,collinear_pairs=str(root/'synteny.tsv'),similarity='unused')
            context=dict(args=args,seqs={g:'M'*20 for g in 'abc'},config=cli.DEFAULTS,
                         edges={('a','b'):(.8,.8,1),('a','c'):(.8,.8,1)},
                         constraints=[dict(gene_a='a',gene_b='c',role='construction',relation='different',weight=1)])
            params=dict(identity=.5,coverage=.8)
            a,_=methods.build_graph_edges('weighted-mcl',params,context)
            b,_=methods.build_graph_edges('weighted-louvain',params,context)
            self.assertEqual(a,b)
            self.assertEqual(a,{('a','b'):1.6,('a','c'):.4})
            args.collinear_pairs=None
            context['constraints']=[]
            self.assertEqual(methods.build_graph_edges('weighted-louvain',params,context)[0],
                             {('a','b'):.8,('a','c'):.8})
            self.assertEqual(context['edges'][('a','b')][0],.8)

    def test_auto_real_subprocess_six_scores_alias_and_resume(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            genes=list('abcdefg')
            data.write_fasta(root/'input.fa',{g:'MACDEFGHIKLMNPQRSTVWY' for g in genes})
            data.tsv(root/'map.tsv',[dict(gene_ID=g,assembly_ID=g) for g in genes],['gene_ID','assembly_ID'])
            edges=[(a,b) for group in ('abc','def') for i,a in enumerate(group) for b in group[i+1:]]
            (root/'hits.tsv').write_text(''.join('%s\t%s\t90\t20\t0\t0\t1\t20\t1\t20\t0\t100\n'%pair for pair in edges))
            p=argparse.ArgumentParser()
            cli.add_parser(p.add_subparsers())
            base=['cluster','-i',str(root/'input.fa'),'--gene-map',str(root/'map.tsv'),
                  '--similarity',str(root/'hits.tsv'),'-o',str(root/'out')]
            args=p.parse_args(base+['-M','auto'])
            # Test-only target installs must also be visible in the real worker process.
            pythonpath=os.pathsep.join(filter(None,[str(Path(nx.__file__).resolve().parents[1]),os.environ.get('PYTHONPATH','')]))
            with patch.dict(os.environ,{'PYTHONPATH':pythonpath}), \
                 patch.object(cli,'tool_versions',return_value={}), \
                 patch.object(methods.shutil,'which',return_value=None):
                selection=cli.run(args)
                scores=data.read_tsv(root/'out'/'scores.tsv')
                self.assertEqual(len(scores),1)
                self.assertEqual(scores[0]['method'],'weighted-louvain')
                self.assertEqual(scores[0]['OGG_count'],'3')
                metric=data.read_tsv(root/'out'/'metric_scores.tsv')
                self.assertEqual(len(metric),6)
                self.assertEqual({r['method'] for r in metric},{'weighted-louvain'})
                manifest=json.loads((root/'out'/'manifest.json').read_text())
                self.assertEqual(len(manifest['louvain_runs']),1)
                graph=Path(manifest['louvain_runs'][0]['graph_file']).read_text()
                self.assertTrue(all(line.split()[0]!=line.split()[1] for line in graph.splitlines()))
                self.assertEqual(manifest['commands'][0]['returncode'],0)
                first=data.read_tsv(root/'out'/'selected_clusters.tsv')
                self.assertEqual(len(first),7)
                args.resume=True
                cli.run(args)
                self.assertEqual(first,data.read_tsv(root/'out'/'selected_clusters.tsv'))
                self.assertEqual(metric,data.read_tsv(root/'out'/'metric_scores.tsv'))
                self.assertEqual(len(json.loads((root/'out'/'manifest.json').read_text())['louvain_runs']),1)
                alias=p.parse_args(base+['-M','gephi'])
                alias.output=str(root/'alias')
                cli.run(alias)
                self.assertEqual(data.read_tsv(root/'alias'/'scores.tsv')[0]['method'],'weighted-louvain')
                self.assertEqual(first,data.read_tsv(next((root/'alias'/'candidates').glob('*/clusters.tsv'))))

    def test_candidate_grid_resolution_validation(self):
        parser=argparse.ArgumentParser()
        cli.add_parser(parser.add_subparsers())
        args=parser.parse_args(['cluster','-i','x','--gene-map','m','-o','out','--louvain-resolution','2'])
        plans=cli.candidates_for(['weighted-louvain'],args,cli.DEFAULTS)
        self.assertEqual(plans[0][1]['resolution'],2)
        config=dict(cli.DEFAULTS,grid={'weighted-louvain':[{'resolution':0}]})
        with self.assertRaises(ValueError):
            cli.candidates_for(['weighted-louvain'],args,config)


if __name__=='__main__':
    unittest.main()
