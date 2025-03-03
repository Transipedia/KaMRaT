import unittest
from os import path, mkdir, listdir
import os
from shutil import rmtree
import subprocess


kamrat = path.join(".", "apps", "kamrat")


class TestCommands(unittest.TestCase):

    def test_basic_command(self):
        output = subprocess.check_output(kamrat, stderr=subprocess.STDOUT, shell=True, text=True)
        self.assertTrue("================== k-mer Matrix, Really Tremendous! ==================" in output)

    def test_index(self):
        test_dir = "index_tmp_test"
        data = path.join("toyroom", "data")

        # Remove previous test remainings
        if path.exists(test_dir):
            rmtree(test_dir)
        mkdir(test_dir)

        # Define inputs and outputs
        intab = path.join(data, "kmer-counts.subset4toy.tsv.gz")
        outdir = path.join(test_dir, "kamrat.idx")
        mkdir(outdir)
        index_stdout = path.join(test_dir, "index.stdout")

        cmd = f"{kamrat} index -intab {intab} -outdir {outdir} -klen 31 -unstrand -nfbase 1000000000"
        # Command ok ?
        process = None
        with open(index_stdout, "w") as idx_out:
            process = subprocess.run(cmd.split(" "), stdout=idx_out, stderr=idx_out)
        self.assertEqual(0, process.returncode)
        # Number of outputs ok ?
        self.assertEqual(3, len(listdir(outdir)))

        # correct outputs ?
        idx_mat = path.join(outdir, "idx-mat.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("e5904e8c5e9cce3be07e3d96a9a91bc0"))
        stream.close()

        idx_mat = path.join(outdir, "idx-meta.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("8036df91ef0fd678b431b65a83863bdb"))
        stream.close()

        idx_mat = path.join(outdir, "idx-pos.bin")
        self.assertTrue(path.exists(idx_mat))
        stream = os.popen(f"md5sum {idx_mat}")
        self.assertTrue(stream.read().startswith("86eca32eb46f0f701a8ba2c56e6393aa"))
        stream.close()

        # Cleaning
        rmtree(test_dir)


    def test_filter(self):
        test_dir = "filter_tmp_test"
        data = path.join("toyroom", "data")

        # Remove previous test remainings
        if path.exists(test_dir):
            rmtree(test_dir)
        mkdir(test_dir)

        # Define inputs and outputs for index
        intab = path.join(data, "kmer-counts.subset4toy.tsv.gz")
        idx_dir = path.join(test_dir, "kamrat.idx")
        mkdir(idx_dir)
        index_stdout = path.join(test_dir, "index.stdout")

        # Index
        cmd = f"{kamrat} index -intab {intab} -outdir {idx_dir} -klen 31 -unstrand -nfbase 1000000000"
        process = None
        with open(index_stdout, "w") as idx_out:
            process = subprocess.run(cmd.split(" "), stdout=idx_out, stderr=idx_out)
        self.assertEqual(0, process.returncode)

        # Define filter i/o
        toy_conditions = path.join(data, "sample-states.toy.tsv")
        conditions = path.join(test_dir, "conditions.tsv")
        filtered = path.join(test_dir, "top-kmers.bin")
        filter_stdout = path.join(test_dir, "filter.stdout")

        # Prepare design file
        with open(conditions, "w") as cdt_out:
            process = subprocess.run(f"sed 's/normal/DOWN/g' {toy_conditions} | sed 's/tumor/UP/g'", shell=True, stdout=cdt_out)
        self.assertEqual(0, process.returncode)
        # Filter
        cmd = f"{kamrat} filter -idxdir {idx_dir} -design {conditions} -upmin 5:5 -downmax 0:10 -outpath {filtered}"
        with open(filter_stdout, "w") as flt_out:
            process = subprocess.run(cmd.split(" "), stdout=flt_out, stderr=flt_out)
        self.assertEqual(0, process.returncode)

        self.assertTrue(path.exists(filtered))
        stream = os.popen(f"md5sum {filtered}")
        self.assertTrue(stream.read().startswith("6e43e7e58362415a1e252c98f5b8dd7b"))
        stream.close()

        rmtree(test_dir)

    def test_merge(self):
        test_dir = "merge_tmp_test"
        data = path.join("toyroom", "data")

        # Remove previous test remainings
        if path.exists(test_dir):
            rmtree(test_dir)
        mkdir(test_dir)

        # Define inputs and outputs for index
        intab = path.join(data, "kmer-counts.subset4toy.tsv.gz")
        idx_dir = path.join(test_dir, "kamrat.idx")
        mkdir(idx_dir)
        index_stdout = path.join(test_dir, "index.stdout")

        # Index
        cmd = f"{kamrat} index -intab {intab} -outdir {idx_dir} -klen 31 -unstrand -nfbase 1000000000"
        process = None
        with open(index_stdout, "w") as idx_out:
            process = subprocess.run(cmd.split(" "), stdout=idx_out, stderr=idx_out)
        self.assertEqual(0, process.returncode)

        # Define filter i/o
        toy_conditions = path.join(data, "sample-states.toy.tsv")
        conditions = path.join(test_dir, "conditions.tsv")
        filtered = path.join(test_dir, "top-kmers.bin")
        filter_stdout = path.join(test_dir, "filter.stdout")

        # Prepare design file
        with open(conditions, "w") as cdt_out:
            process = subprocess.run(f"sed 's/normal/DOWN/g' {toy_conditions} | sed 's/tumor/UP/g'", shell=True, stdout=cdt_out)
        self.assertEqual(0, process.returncode)
        # Filter
        cmd = f"{kamrat} filter -idxdir {idx_dir} -design {conditions} -upmin 5:5 -downmax 0:10 -outpath {filtered}"
        with open(filter_stdout, "w") as flt_out:
            process = subprocess.run(cmd.split(" "), stdout=flt_out, stderr=flt_out)
        self.assertEqual(0, process.returncode)

        self.assertTrue(path.exists(filtered))
        stream = os.popen(f"md5sum {filtered}")
        self.assertTrue(stream.read().startswith("6e43e7e58362415a1e252c98f5b8dd7b"))
        stream.close()

        # Define merge i/o
        merged = path.join(test_dir, "top-ctg-counts.tsv")

        # Merge
        cmd = f"{kamrat} merge -idxdir {idx_dir} -overlap 30-15 -with {filtered}:min -outpath {merged} -withcounts mean"
        with open(filter_stdout, "w") as flt_out:
            process = subprocess.run(cmd.split(" "), stdout=flt_out, stderr=flt_out)
        self.assertEqual(0, process.returncode)

        self.assertTrue(path.exists(merged))
        stream = os.popen(f"md5sum {merged}")
        self.assertTrue(stream.read().startswith("10e910446f777f30922f97c148a0404a"))
        stream.close()

        rmtree(test_dir)


    def test_rank(self):
        test_dir = "rank_tmp_test"
        data = path.join("toyroom", "data")

        # Remove previous test remainings
        if path.exists(test_dir):
            rmtree(test_dir)
        mkdir(test_dir)

        # Define inputs and outputs for index
        intab = path.join(data, "kmer-counts.subset4toy.tsv.gz")
        idx_dir = path.join(test_dir, "kamrat.idx")
        mkdir(idx_dir)
        index_stdout = path.join(test_dir, "index.stdout")

        # Index
        cmd = f"{kamrat} index -intab {intab} -outdir {idx_dir} -klen 31 -unstrand -nfbase 1000000000"
        process = None
        with open(index_stdout, "w") as idx_out:
            process = subprocess.run(cmd.split(" "), stdout=idx_out, stderr=idx_out)
        self.assertEqual(0, process.returncode)

        # Define rank i/o
        condition = path.join(data, "sample-states.toy.tsv")
        ranked = path.join(test_dir, "top-kmers.bin")
        rank_out = path.join(test_dir, "rank.stdout")
        
        # Rank
        cmd = f"{kamrat} score -idxdir {idx_dir} -scoreby ttest.padj -design {condition} -seltop 0.1 -outpath {ranked}"
        with open(rank_out, "w") as rk_out:
            process = subprocess.run(cmd.split(" "), stdout=rk_out, stderr=rk_out)
        self.assertEqual(0, process.returncode)

        self.assertTrue(path.exists(ranked))
        stream = os.popen(f"md5sum {ranked}")
        self.assertTrue(stream.read().startswith("faf89801a833751a64c02cb9db88770f"))
        stream.close()

        rmtree(test_dir)

if __name__ == '__main__':
  unittest.main()