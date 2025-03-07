import unittest
import os
import subprocess
import shutil


KAMRAT = f"./apps/kamrat"

DATADIR = "./toyroom/data/"
COUNTTAB = f"{DATADIR}kmer-counts.subset4toy.tsv.gz"
NFFILE = f"{DATADIR}nf_values.txt"
FASTA = f"{DATADIR}sequence.toy.fa"
SMPFLT = f"{DATADIR}sample-indications.toy.tsv"
DSGN = f"{DATADIR}sample-states.toy.tsv"

OUTDIR = "./tests/tmp_test/"
if not os.path.exists(OUTDIR):
    os.mkdir(OUTDIR)


class TestCommands(unittest.TestCase):

    def test_basic_command(self):
        output = subprocess.check_output(
            KAMRAT, stderr=subprocess.STDOUT, shell=True, text=True
        )
        self.assertTrue(
            "================== k-mer Matrix, Really Tremendous! =================="
            in output
        )

    def test_index(self):
        # Test folder
        idxdir = f"{OUTDIR}index/"
        if os.path.exists(idxdir):
            shutil.rmtree(idxdir)
        # Arguments, outputs, and their md5sums
        arg_out_md5 = {
            "": {  # no normalisation
                "idx-meta.bin": "c6ac53928b1cd9ce946de29496684cd0",
                "idx-pos.bin": "86eca32eb46f0f701a8ba2c56e6393aa",
                "idx-mat.bin": "654fdaf169953e8aa6c191b2367690dd",
            },
            "-nfbase 1000000": {  # normalisation with given base value
                "idx-meta.bin": "2767605d73fac371b5c027d0a9af2644",
                "idx-pos.bin": "86eca32eb46f0f701a8ba2c56e6393aa",
                "idx-mat.bin": "728e1d2b099a1e2b92ac782c6929e6a9",
            },
            f"-nffile {NFFILE}": {  # normalisation with given normalisation factor
                "idx-meta.bin": "2767605d73fac371b5c027d0a9af2644",
                "idx-pos.bin": "86eca32eb46f0f701a8ba2c56e6393aa",
                "idx-mat.bin": "728e1d2b099a1e2b92ac782c6929e6a9",
            },
        }
        for i, arg in enumerate(arg_out_md5.keys()):
            os.mkdir(idxdir)
            cmd = f"{KAMRAT} index -intab {COUNTTAB} -outdir {idxdir} -klen 31 -unstrand {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}index{i + 1}.log", "w") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            self.assertEqual(3, len(os.listdir(idxdir)))  # 3 output files
            ## Test if the outputs are correct
            for out, md5 in arg_out_md5[arg].items():
                out_path = f"{idxdir}{out}"
                self.assertTrue(os.path.exists(out_path))
                stream = os.popen(f"md5sum {out_path}")
                self.assertTrue(stream.read().split(" ")[0] == md5)
                stream.close()
            ## Cleaning test folder
            shutil.rmtree(idxdir)

    def test_mask(self):
        # Test folder
        mskdir = f"{OUTDIR}mask/"
        if os.path.exists(mskdir):
            shutil.rmtree(mskdir)
        # Arguments, outputs, and their md5sums
        arg_md5 = {
            "": "083e48a1b7e599da135d5ce6bb68e81a",  # no reverse, tab, int
            "-counts float": "36c864675446d253d1128afe2f668f15",  # no reverse, tab, float
            "-reverse -counts float": "f08906abb46f8ef9081fff2552985074",  # reverse, tab, float
        }
        # Run index
        os.mkdir(mskdir)
        idx = f"{KAMRAT} index -intab {COUNTTAB} -outdir {mskdir} -klen 31 -unstrand -nffile {NFFILE}"
        process = None
        outerr = subprocess.DEVNULL
        process = subprocess.run(idx.split(" "), stdout=outerr, stderr=outerr)
        # Test mask
        for i, arg in enumerate(arg_md5.keys()):
            if os.path.exists(f"{mskdir}out"):
                os.remove(f"{mskdir}out")
            cmd = f"{KAMRAT} mask -idxdir {mskdir} -fasta {FASTA} -outpath {mskdir}out {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}mask{i + 1}.log", "w") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            ## Test if the outputs are correct
            md5 = arg_md5[arg]
            self.assertTrue(os.path.exists(f"{mskdir}out"))
            stream = os.popen(f"md5sum {mskdir}out")
            self.assertTrue(stream.read().split(" ")[0] == md5)
            stream.close()
        ## Cleaning test folder
        shutil.rmtree(mskdir)

    def test_filter(self):
        # Test folder
        fltdir = f"{OUTDIR}filter/"
        if os.path.exists(fltdir):
            shutil.rmtree(fltdir)
        # Arguments, outputs, and their md5sums
        arg_md5 = {
            "-upmin 5:5 -downmax 0:5": "6f41031fe656f9d24c13a2abca2fc83e",  # tab, int
            "-upmin 5:5 -downmax 0:5 -counts float": "07dced55477e8f41357bb1ec78a8f38b",  # tab, float
            "-upmin 5:5 -downmax 0:5 -reverse": "2e1d65280faca08b9df60b9cc3a0ca54",  # reverse, tab, int
            "-upmin 5:5 -downmax 0:5 -outfmt fa": "3082d4af223c0f2581f2620c2881ef40",  # fa
            "-upmin 5:5 -downmax 0:5 -outfmt bin": "1a32d48128b2cca75382b7cd9886f755",  # bin
        }
        # Run index
        os.mkdir(fltdir)
        idx = f"{KAMRAT} index -intab {COUNTTAB} -outdir {fltdir} -klen 31 -unstrand -nffile {NFFILE}"
        process = None
        outerr = subprocess.DEVNULL
        process = subprocess.run(idx.split(" "), stdout=outerr, stderr=outerr)
        # Test filter
        for i, arg in enumerate(arg_md5.keys()):
            if os.path.exists(f"{fltdir}out"):
                os.remove(f"{fltdir}out")
            cmd = f"{KAMRAT} filter -idxdir {fltdir} -design {SMPFLT} -outpath {fltdir}out {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}filter{i + 1}.log", "w") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            ## Test if the outputs are correct
            md5 = arg_md5[arg]
            self.assertTrue(os.path.exists(f"{fltdir}out"))
            stream = os.popen(f"md5sum {fltdir}out")
            self.assertTrue(stream.read().split(" ")[0] == md5)
            stream.close()
        ## Cleaning test folder
        shutil.rmtree(fltdir)

    def test_merge(self):
        # Test folder
        mgdir = f"{OUTDIR}merge/"
        if os.path.exists(mgdir):
            shutil.rmtree(mgdir)
        # Arguments, outputs, and their md5sums
        arg_md5 = {
            "": "d4078093f865ee9615f2bcb7e78143e7",  # tab, rep, int
            "-counts mean": "0fb8cf5abf86ee0c6480e2238c972f17",  # tab, mean, int
            "-counts median": "7993a9cc071d83f5c8d13927a432c3bd",  # tab, median, int
            "-counts mean:float": "7e02083c2562e77418b16c7da88e8dd9",  # tab, mean, float
            "-outfmt fa": "8e74e45ae8ab3434a265ab9d845c9498",  # fasta
            f"-with {mgdir}flt:min -counts mean:float": "7e16baf588cb3cbc172f4b86d5315260",  # with, mean, float
        }
        # Run index
        os.mkdir(mgdir)
        idx = f"{KAMRAT} index -intab {COUNTTAB} -outdir {mgdir} -klen 31 -unstrand -nffile {NFFILE}"
        process = None
        outerr = subprocess.DEVNULL
        process = subprocess.run(idx.split(" "), stdout=outerr, stderr=outerr)
        # Test merge
        for i, arg in enumerate(arg_md5.keys()):
            if os.path.exists(f"{mgdir}out"):
                os.remove(f"{mgdir}out")
            if arg.startswith("-with"):
                ## Run filter
                flt = f"{KAMRAT} filter -idxdir {mgdir} -design {SMPFLT} -upmin 5:5 -downmax 0:5 -outfmt bin -outpath {mgdir}flt"
                process = None
                outerr = subprocess.DEVNULL
                process = subprocess.run(flt.split(" "), stdout=outerr, stderr=outerr)
            cmd = f"{KAMRAT} merge -idxdir {mgdir} -overlap 30-15 -outpath {mgdir}out {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}merge{i + 1}.log", "a") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            ## Test if the outputs are correct
            md5 = arg_md5[arg]
            self.assertTrue(os.path.exists(f"{mgdir}out"))
            stream = os.popen(f"md5sum {mgdir}out")
            self.assertTrue(stream.read().split(" ")[0] == md5)
            stream.close()
        ## Cleaning test folder
        shutil.rmtree(mgdir)

    def test_score(self):
        # Test folder
        scdir = f"{OUTDIR}score/"
        if os.path.exists(scdir):
            shutil.rmtree(scdir)
        # Arguments, outputs, and their md5sums
        arg_md5 = {
            "": "6c4b831d47741122fd9854b07baea329",  # tab, int
            "-counts float": "8861fcecd8288213551768dc3723b13f",  # tab, float
            "-outfmt fa": "d92384e661a41708e44ac8d07ee0a2a0",  # fasta
            "-outfmt bin": "9ec2639ff71a21d164ffb9fb93c35f54",  # bin
            f"-with {scdir}mg:mean -counts int": "d7bf39494b7a9526cba8cafb8297ff3e",  # with, mean, int
        }
        # Run index
        os.mkdir(scdir)
        idx = f"{KAMRAT} index -intab {COUNTTAB} -outdir {scdir} -klen 31 -unstrand -nffile {NFFILE}"
        process = None
        outerr = subprocess.DEVNULL
        process = subprocess.run(idx.split(" "), stdout=outerr, stderr=outerr)
        # Test score
        for i, arg in enumerate(arg_md5.keys()):
            if os.path.exists(f"{scdir}out"):
                os.remove(f"{scdir}out")
            if arg.startswith("-with"):
                ## Run filter
                flt = f"{KAMRAT} merge -idxdir {scdir} -overlap 30-15 -outpath {scdir} -outfmt bin -outpath {scdir}mg"
                process = None
                outerr = subprocess.DEVNULL
                process = subprocess.run(flt.split(" "), stdout=outerr, stderr=outerr)
            cmd = f"{KAMRAT} score -idxdir {scdir} -scoreby ttest.padj -design {DSGN} -seltop 0.1 -outpath {scdir}out {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}score{i + 1}.log", "a") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            ## Test if the outputs are correct
            md5 = arg_md5[arg]
            self.assertTrue(os.path.exists(f"{scdir}out"))
            stream = os.popen(f"md5sum {scdir}out")
            self.assertTrue(stream.read().split(" ")[0] == md5)
            stream.close()
        ## Cleaning test folder
        shutil.rmtree(scdir)

    def test_query(self):
        # Test folder
        qydir = f"{OUTDIR}query/"
        if os.path.exists(qydir):
            shutil.rmtree(qydir)
        # Arguments, outputs, and their md5sums
        arg_md5 = {
            "-toquery mean -withabsent": "157b4014960f6a977f0d303524ce5fac",  # mean, withabsent, int
            "-toquery mean -counts float": "64e067d3e600deef8f7acfb30eaca495",  # mean, float
            "-toquery median -counts float": "05858c862941ca53f57d670e0106a2a0",  # median, float
        }
        # Run index
        os.mkdir(qydir)
        idx = f"{KAMRAT} index -intab {COUNTTAB} -outdir {qydir} -klen 31 -unstrand -nffile {NFFILE}"
        process = None
        outerr = subprocess.DEVNULL
        process = subprocess.run(idx.split(" "), stdout=outerr, stderr=outerr)
        # Test merge
        for i, arg in enumerate(arg_md5.keys()):
            if os.path.exists(f"{qydir}out"):
                os.remove(f"{qydir}out")
            cmd = f"{KAMRAT} query -idxdir {qydir} -fasta {FASTA} -outpath {qydir}out {arg}".rstrip()
            ## Test if command runs
            process = None
            with open(f"{OUTDIR}query{i + 1}.log", "a") as outerr:
                process = subprocess.run(cmd.split(" "), stdout=outerr, stderr=outerr)
            self.assertEqual(0, process.returncode)  # normal exit code
            ## Test if the outputs are correct
            md5 = arg_md5[arg]
            self.assertTrue(os.path.exists(f"{qydir}out"))
            stream = os.popen(f"md5sum {qydir}out")
            self.assertTrue(stream.read().split(" ")[0] == md5)
            stream.close()
        ## Cleaning test folder
        shutil.rmtree(qydir)


if __name__ == "__main__":
    unittest.main()
