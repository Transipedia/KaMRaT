import unittest
from os import path, mkdir, listdir
import os
from shutil import rmtree
import subprocess


kamrat = path.join(".", "bin", "kamrat")


class TestCommands(unittest.TestCase):

    def test_basic_command(self):
        output = subprocess.check_output(f"./bin/kamrat", stderr=subprocess.STDOUT, shell=True, text=True)
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
        self.assertTrue(stream.read().startswith("c6ac53928b1cd9ce946de29496684cd0"))
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
        toy_conditions = path.join(data, "sample-condition.toy.tsv")
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
        toy_conditions = path.join(data, "sample-condition.toy.tsv")
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
        self.assertTrue(stream.read().startswith("54fba78dade59136999f8f5e1dfb6b08"))
        stream.close()

        rmtree(test_dir)

if __name__ == '__main__':
  unittest.main()