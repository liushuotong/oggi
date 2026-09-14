"""Louvain subprocess so the existing candidate/global timeouts apply.

Gephi's built-in Modularity uses the same algorithm family. This is NetworkX's
implementation, with its gamma convention, not a bit-for-bit Gephi port.
"""
import argparse
import json
import math
from pathlib import Path
import networkx as nx


def cluster(genes, edges, resolution=1., seed=20260914):
    if not math.isfinite(resolution) or resolution <= 0:
        raise ValueError('Louvain resolution must be finite and positive')
    if len(genes) != len(set(genes)):
        raise ValueError('duplicate target node IDs')
    graph = nx.Graph()
    graph.add_nodes_from(sorted(genes))
    for a,b,w in sorted(edges):
        if a not in graph or b not in graph or a == b or graph.has_edge(a,b):
            raise ValueError('unknown node, self-loop or duplicate undirected edge')
        if not math.isfinite(w) or w <= 0:
            raise ValueError('edge weights must be finite and positive')
        graph.add_edge(a,b,weight=w)
    if not math.isfinite(graph.size(weight='weight')):
        raise ValueError('total edge weight overflow')
    if graph.number_of_edges():
        groups = nx.community.louvain_communities(graph, weight='weight', resolution=resolution,
                                                  threshold=1e-7, seed=seed)
        objective = nx.community.modularity(graph, groups, weight='weight', resolution=resolution)
        standard = nx.community.modularity(graph, groups, weight='weight', resolution=1)
    else:
        groups = [{g} for g in sorted(genes)]
        objective = standard = None
    reciprocal = 1/resolution
    return dict(groups=sorted(sorted(g) for g in groups), metadata=dict(
        implementation='networkx.community.louvain_communities', networkx_version=nx.__version__,
        node_count=len(genes), edge_count=graph.number_of_edges(), seed=seed,
        resolution_networkx=resolution, gephi_resolution_objective_equivalent=reciprocal if math.isfinite(reciprocal) else None,
        modularity_gain_threshold=1e-7, artificial_self_loops=False,
        construction_objective=objective, construction_modularity_gamma1=standard,
        note='Construction objective only; auto evaluation uses its separate common inputs. '
             'The reciprocal matches objective conventions, not guaranteed identical Gephi partitions.'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',required=True)
    parser.add_argument('--graph',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--resolution',type=float,default=1.)
    parser.add_argument('--seed',type=int,default=20260914)
    args=parser.parse_args()
    genes=json.loads(Path(args.nodes).read_text(encoding='utf-8'))
    edges=[]
    with Path(args.graph).open(encoding='utf-8') as handle:
        for line in handle:
            if line.strip():
                a,b,w=line.rstrip('\n').split('\t')
                edges.append((a,b,float(w)))
    result=cluster(genes,edges,args.resolution,args.seed)
    Path(args.output).write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print('Louvain: %d nodes, %d edges, %d groups; gamma=%g; seed=%d' %
          (len(genes),len(edges),len(result['groups']),args.resolution,args.seed))


if __name__=='__main__':
    main()
