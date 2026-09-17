"""Regression checks for independently verified audit findings; tools are mocked."""
import argparse
import importlib.util
import json
import pathlib
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import collinearity_matrix as cm
import wgdi_all_vs_all as wg
import sub_collinearity as sc
from cluster_engine import scoring, evidence


class ParserRegressionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=pathlib.Path(self.temp.name)

    def test_scientific_notation_and_mcscanx(self):
        file=self.root/'wgdi.txt'
        file.write_text('# Alignment 1: score=2.9e2 pvalue=1.16e-08 N=6 A&B minus\na 0 b 5\n')
        table=wg.parse_wgdi_collinearity(file)
        self.assertEqual(len(table),1)
        self.assertEqual(table.iloc[0].pvalue,1.16e-8)
        file.write_text('## Alignment 0: score=290 e_value=1e-10 N=6 A&B minus\n  0-0: a b 1e-50\n')
        self.assertEqual(cm.parse_collinearity_pairs(file).values.tolist(),[['a','b']])
        file.write_text('## Alignment 0: score=290 e_value=1e-10 N=6 A&B minus\n  0-  0: a b 1e-50\n')
        self.assertEqual(cm.parse_collinearity_pairs(file).values.tolist(),[['a','b']])

    def test_duplicate_bed_is_explicit_error_not_guessed(self):
        file=self.root/'a.bed'
        file.write_text('chr1\t1\t5\tx\nchr2\t2\t6\tx\n')
        for reader in (sc.load_bed,wg.read_agat_bed):
            with self.assertRaisesRegex(ValueError,'ambiguous'):
                reader(file)

    def test_nonfinite_and_empty_scoring_rejected(self):
        with self.assertRaises(ValueError):
            scoring.totals([dict(B=float('nan'),R=1,A=1,Q=1)],scoring.WEIGHTS)
        with self.assertRaises(ValueError):
            evidence.boundary({},[])

    def test_unknown_constraint_valueerror(self):
        file=self.root/'e.tsv'
        file.write_text('gene_a\tgene_b\trelation\tweight\tsource\tblock_id\trole\na\tunknown\tsame\t1\tx\tb\tevaluation\n')
        with self.assertRaisesRegex(ValueError,'non-target'):
            evidence.read_constraints(file,['a'],'hog','hog')

    def test_minus_blocks_and_nonprefix_ids(self):
        rows=[]
        hits=[]
        membership={}
        for asm,offset in [('A',0),('B',6)]:
            file=self.root/(asm+'.bed')
            file.write_text(''.join('chr1\t%d\t%d\tx%d\n'%(i*100,i*100+80,i+offset) for i in range(1,7)))
            rows.append(dict(assembly=asm,pep='unused',bed=str(file)))
        # IDs deliberately share the same prefix; only centers are in target map.
        membership={'x3':'A','x10':'B'}
        for i in range(1,7):
            hits.append('x%d\tx%d\t95\t20\t0\t0\t1\t20\t1\t20\t1e-40\t100\n'%(i,13-i))
        file=self.root/'hits'
        file.write_text(''.join(hits))
        out=self.root/'blocks'
        result=sc.batch_member_pair_collinearity(rows,membership,str(file),up=6,down=6,
                      pvalue_accept=.2,blocks_out=str(out),verbose=False)
        self.assertEqual(len(result),1)
        self.assertTrue(result.iloc[0].in_collinear_block)
        self.assertIn('minus',out.read_text())
        self.assertEqual(set(map(tuple,cm.parse_collinearity_pairs(out).values)),
                         {tuple(sorted(('x%d'%i,'x%d'%(13-i)))) for i in range(1,7)})


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Biopython absent: upstream integration tests require Bio')
class UpstreamIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=pathlib.Path(self.temp.name)

    def test_mcscanx_pair_adapter_no_name_collision(self):
        import mcscan_all_vs_all as m
        with patch.object(m,'_match_stems',return_value=[('A','A.gff','A.fa'),('B','B.gff','B.fa')]), \
             patch.object(m.pre,'sub_collinearity_gff_process',side_effect=lambda g,s:(g,'cds',s,s+'.bed')), \
             patch.object(m,'bed_to_gff_for_mcscanx'), patch.object(m,'mcscanx_BLASTP'), \
             patch.object(m,'run_mcscanx_pair') as pair:
            names,prefixes=m.run_mcscanx('gff','fa',verbose=False)
        self.assertEqual(names,['A','B'])
        pair.assert_called_once()

    def test_reduce_identify_subcoli_chain_with_mocked_external_tools(self):
        import oggi
        import sub_collinearity_pre_process as pre
        import gene_family_identification as gfi
        gffs,genomes,processed=[self.root/n for n in ('gff','genomes','processed')]
        for d in (gffs,genomes,processed): d.mkdir()
        for asm in ['A','B']:
            (gffs/(asm+'.gff')).write_text('##gff-version 3\n')
            (genomes/(asm+'.fa')).write_text('>chr1\nATGATGATG\n')
        def fake_agat(cmd):
            out=pathlib.Path(cmd[cmd.index('-o')+1])
            asm=out.stem
            if out.suffix=='.bed':
                out.write_text(''.join('chr1\t%d\t%d\t%s%d\n'%(i*100,i*100+80,asm,i) for i in range(1,7)))
            elif out.suffix=='.pep':
                out.write_text(''.join('>%s%d\n%s\n'%(asm,i,'M'*20) for i in range(1,7)))
            else: out.write_text('mock external AGAT output\n')
        reduce=argparse.Namespace(gff_dir=str(gffs),genome_dir=str(genomes),output=str(processed),skip_existing=False)
        with patch.object(pre,'run_agat',side_effect=fake_agat):
            oggi.run_reduce(reduce)
            (processed/'A.bed').unlink()
            reduce.skip_existing=True
            oggi.run_reduce(reduce)
        self.assertTrue((processed/'A.bed').exists())
        family=self.root/'family'
        with patch.object(gfi,'identification_caculation',side_effect=lambda h,r,eh,eb,fa,out,asm,cpu:{asm+'3'}):
            mapping,ids,fa=gfi.main_identification([['A',str(processed/'A.pep')],['B',str(processed/'B.pep')]],[],[],1e-5,1e-5,str(family),1)
        blast=self.root/'hits'
        blast.write_text(''.join('A%d\tB%d\t95\t20\t0\t0\t1\t20\t1\t20\t1e-40\t100\n'%(i,i) for i in range(1,7)))
        args=argparse.Namespace(manifest=str(processed/'assembly_manifest.tsv'),id_table=str(family/'gene_to_assembly.tsv'),
                   blast=str(blast),output=str(self.root/'subcoli'),up=6,down=6,evalue=1e-5,pvalue=.2,max_pairs=None,threads=1)
        oggi.run_subcoli(args)
        self.assertTrue((self.root/'subcoli.collinear_pairs.tsv').stat().st_size>0)
        with patch.object(gfi,'identification_caculation',return_value=set()):
            with self.assertRaisesRegex(ValueError,'no family genes'):
                gfi.main_identification([['A',str(processed/'A.pep')]],[],[],1e-5,1e-5,str(self.root/'empty'),1)


if __name__=='__main__': unittest.main()
