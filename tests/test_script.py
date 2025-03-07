import unittest
import os
import subprocess
import shutil


KAMRAT = f"./apps/kamrat"

DATADIR = "./toyroom/data/"
COUNTTAB = f"{DATADIR}kmer-counts.subset4toy.tsv.gz"
NFFILE = f"{DATADIR}nf_values2.txt"
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
            "-nfbase 1000000000": {
                "idx-meta.bin": "8036df91ef0fd678b431b65a83863bdb",
                "idx-pos.bin": "86eca32eb46f0f701a8ba2c56e6393aa",
                "idx-mat.bin": "e5904e8c5e9cce3be07e3d96a9a91bc0",
            },
            f"-nffile {NFFILE}": {
                "idx-meta.bin": "8036df91ef0fd678b431b65a83863bdb",
                "idx-pos.bin": "86eca32eb46f0f701a8ba2c56e6393aa",
                "idx-mat.bin": "e5904e8c5e9cce3be07e3d96a9a91bc0",
            },
        }
        for i, arg in enumerate(arg_out_md5.keys()):
            os.mkdir(idxdir)
            cmd = f"{KAMRAT} index -intab {COUNTTAB} -outdir {idxdir} -klen 31 -unstrand {arg}"
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
            "": "61f72f7489cab2f3053919bdfa224f59",
            "-counts float": "fffe5be701f4002cd72c870bb059daa3",
            "-reverse": "0ecea4551994dc67d8bbba52b4102c93",
            "-reverse -counts float": "4d9f7ffa5b2100252d554f3fc7f7282f",
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
            "-upmin 5:5 -downmax 0:10": "8f841aa7bd3e9560c0b62955e19cfdcb",  # tab, int
            "-upmin 5:5 -downmax 0:10 -counts float": "f609ef4f8f4799c09e37174f49fe83af",  # tab, float
            "-upmin 5:5 -downmax 0:10 -reverse": "a58f74274658d08f6e0be75d93f34a45",  # reverse, tab, int
            "-upmin 5:5 -downmax 0:10 -outfmt fa": "28300e219e9e9591cb8d019a4398bb22",  # fa
            "-upmin 5:5 -downmax 0:10 -outfmt bin": "6e43e7e58362415a1e252c98f5b8dd7b",  # bin
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
            "": "90978714b0f6b397f4b9ee71905d5099",  # tab, rep, int
            "-counts mean": "279d58144fe77a90e5ee2e04e9e9a5b3",  # tab, mean, int
            "-counts median": "331a9af85afc4e2d14685945b16e5bea",  # tab, median, int
            "-counts mean:float": "a6f576c240b127e1d4caf12a03c9565d",  # tab, mean, float
            "-outfmt fa": "8e74e45ae8ab3434a265ab9d845c9498",  # fasta
            f"-with {mgdir}flt:min -counts mean:float": "10e910446f777f30922f97c148a0404a",  # with, mean, float
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
                flt = f"{KAMRAT} filter -idxdir {mgdir} -design {SMPFLT} -upmin 5:5 -downmax 0:10 -outfmt bin -outpath {mgdir}flt"
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
            "": "dc15382ae02524a260c8b8fa6ab26946",  # tab, int
            "-counts float": "f07ec39e3d69d8bd453fd22f82c214ef",  # tab, float
            "-outfmt fa": "3e5e88ed24936a33af47ebef87729902",  # fasta
            "-outfmt bin": "faf89801a833751a64c02cb9db88770f",  # bin
            f"-with {scdir}mg -counts int": "24e7d45f19b63638e196c2fe79f30bc8",  # with, rep, int
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
            "-toquery mean -withabsent": "33f375e98e717122f352b25a6afbffb8",  # mean, withabsent, int
            "-toquery mean -counts float": "a034b8b18dab962938ec8c9a1562f55a",  # mean, float
            "-toquery median -counts float": "ff72412b200a4c9780fc9133c18273b6",  # median, float
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
