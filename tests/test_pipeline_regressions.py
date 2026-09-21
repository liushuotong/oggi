"""Regression checks for independently verified audit findings; tools are mocked."""
import argparse
import importlib.util
import json
import os
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
        def fake_agat(cmd,**kwargs):
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


class ReduceFastModeTests(unittest.TestCase):
    """`oggi reduce --fast` must drive the Rust engine and never the AGAT path."""

    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=pathlib.Path(self.temp.name)
        self.gffs=self.root/'gff'; self.genomes=self.root/'genomes'; self.out=self.root/'out'
        for d in (self.gffs,self.genomes): d.mkdir()
        (self.gffs/'A.gff').write_text('##gff-version 3\n')
        (self.genomes/'A.fa').write_text('>chr1\n'+'ATG'*40+'\n')

    def _args(self,**over):
        ns=dict(gff_dir=str(self.gffs),genome_dir=str(self.genomes),output=str(self.out),
                skip_existing=False,fast_mode=True,reduce_rs=None)
        ns.update(over)
        return argparse.Namespace(**ns)

    def test_only_fast_enables_rust_even_with_engine_environment(self):
        import oggi
        import sub_collinearity_pre_process as pre
        parser = argparse.ArgumentParser()
        oggi.add_reduce_parser(parser.add_subparsers())
        argv = ['reduce', '--gff-dir', str(self.gffs), '--genome-dir', str(self.genomes), '-o', str(self.out)]
        with patch.dict(os.environ, {oggi.REDUCE_RS_ENV: 'nonexistent'}), \
             patch.object(oggi, 'require_reduce_rs_engine') as rust, \
             patch.object(pre, 'run_agat') as perl:
            args = parser.parse_args(argv)
            self.assertFalse(args.fast_mode)
            oggi.run_reduce(args)
            rust.assert_not_called()
            self.assertEqual(perl.call_count, 3)
        self.assertTrue(parser.parse_args(argv + ['--fast']).fast_mode)
        self.assertFalse(parser.parse_args(argv + ['--reduce-rs', 'nonexistent']).fast_mode)
        for old in ['--no-fast', '--fast-mode', '--no-fast-mode']:
            with self.assertRaises(SystemExit):
                parser.parse_args(argv + [old])

    def test_windows_perl_is_launched_explicitly(self):
        import sub_collinearity_pre_process as pre
        paths = {'perl': r'C:\perl\perl.exe', 'agat_sp_keep_longest_isoform.pl': r'C:\AGAT tools\agat_sp_keep_longest_isoform.pl'}
        with patch.object(pre.os, 'name', 'nt'), \
             patch.object(pre.shutil, 'which', side_effect=paths.get), \
             patch.object(pre.subprocess, 'run') as run, \
             patch.object(pre, '_cleanup_agat_logs'):
            pre.run_agat(['agat_sp_keep_longest_isoform.pl', '--gff', 'file with spaces.gff'], quiet=True)
        self.assertEqual(run.call_args.args[0], [paths['perl'], paths['agat_sp_keep_longest_isoform.pl'], '--gff', 'file with spaces.gff'])
        self.assertFalse(run.call_args.kwargs['shell'])

    def test_rust_build_cache_changes_with_source(self):
        import oggi
        workspace = self.root / 'sources'
        workspace.mkdir()
        source = workspace / 'Cargo.toml'
        source.write_text('first')
        with patch.dict(os.environ, {'OGGI_REDUCE_CACHE': str(self.root / 'cache')}):
            first = oggi._reduce_rs_cache(str(workspace))
            self.assertEqual(first, oggi._reduce_rs_cache(str(workspace)))
            source.write_text('second')
            self.assertNotEqual(first, oggi._reduce_rs_cache(str(workspace)))

    def test_fast_mode_uses_rust_engine_and_skips_fasta_wrap(self):
        import oggi
        import sub_collinearity_pre_process as pre
        calls=[]
        def fake(cmd,shell=False,quiet=False): calls.append(cmd)
        def fake_fast(engine,keep,extract,bed,progress=None):
            calls.append([engine or 'PATH',keep,extract,bed])
            pathlib.Path(keep).write_text('rust gff\n')
            pathlib.Path(extract[2]).write_text('>A1\nMMM\n')
            pathlib.Path(bed).write_text('chr1\t0\t9\tA1\n')
        with patch.object(oggi,'find_reduce_rs_engine',return_value=r'X:\engine') as finder, \
             patch.object(oggi,'_run_reduce_fast',side_effect=fake_fast), \
             patch.object(pre,'run_agat',side_effect=fake):
            oggi.run_reduce(self._args(reduce_rs=r'X:\explicit'))
        # --reduce-rs must reach the discovery function
        self.assertEqual(finder.call_args[0][0],r'X:\explicit')
        # the perl engine must not run at all, and the Bio::DB::Fasta
        # workaround must not rewrite the genome in fast mode
        self.assertEqual(calls[0][0],r'X:\engine')
        self.assertEqual(len(calls),1)
        for ext in ('.gff','.pep','.bed'):
            self.assertTrue((self.out/('A'+ext)).is_file())

    def test_missing_engine_exits_with_build_instructions(self):
        import oggi
        with patch.object(oggi,'find_reduce_rs_engine',return_value=None):
            with self.assertRaises(SystemExit) as ctx:
                oggi.run_reduce(self._args())
        msg=str(ctx.exception)
        self.assertIn('--fast',msg)
        self.assertIn('cargo build --release',msg)
        self.assertFalse((self.out/'assembly_manifest.tsv').exists())

    def test_explicit_engine_dir_is_authoritative(self):
        """A bad --reduce-rs / $OGGI_REDUCE_RS must error, not silently use
        whatever engine happens to be installed."""
        import oggi
        empty=self.root/'empty-engine'; empty.mkdir()
        with self.assertRaises(SystemExit) as ctx:
            oggi.find_reduce_rs_engine(explicit=str(empty))
        self.assertIn('--reduce-rs',str(ctx.exception))
        with patch.dict(os.environ,{oggi.REDUCE_RS_ENV:str(empty)}):
            with self.assertRaises(SystemExit) as ctx:
                oggi.find_reduce_rs_engine()
        self.assertIn(oggi.REDUCE_RS_ENV,str(ctx.exception))
        # a directory holding all three binaries is accepted
        for n in oggi.REDUCE_RS_BINARIES:
            (empty/oggi._exe(n)).write_text('')
        self.assertEqual(oggi.find_reduce_rs_engine(explicit=str(empty),
                                                    build=False),str(empty))

    def test_both_engines_drive_equivalent_steps(self):
        """--fast must issue the same three steps with the same arguments as the
        AGAT path; the Rust engine additionally needs --force on keep_longest
        (it refuses to overwrite by default, like AGAT does)."""
        import oggi
        import sub_collinearity_pre_process as pre
        def norm(cmds):
            out=[]
            for c in cmds:
                c=list(c)
                c[0]=os.path.basename(c[0]).replace('.pl','').replace('.exe','')
                out.append(c)
            return out
        fast=[]
        def fake_run(argv,**kw):
            fast.append(list(argv))
            pathlib.Path(argv[argv.index('-o')+1]).write_text('x\n')
            class R: returncode=0
            return R()
        with patch.object(oggi,'find_reduce_rs_engine',return_value=r'X:\eng'), \
             patch.object(oggi.subprocess,'run',side_effect=fake_run):
            oggi.run_reduce(self._args())
        agat=[]
        def fake_agat(cmd,shell=False,quiet=False):
            agat.append(list(cmd))
            pathlib.Path(cmd[cmd.index('-o')+1]).write_text('x\n')
        with patch.object(pre,'run_agat',side_effect=fake_agat):
            oggi.run_reduce(self._args(fast_mode=False))
        self.assertTrue(any('--force' in c for c in fast))
        strip=lambda cs:[ [a for a in c if a!='--force'] for c in cs]
        self.assertEqual(strip(norm(fast)),norm(agat))


if __name__=='__main__': unittest.main()
