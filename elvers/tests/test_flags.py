#!/usr/bin/env python
from unittest import TestCase
from .const import (here, elvers_cmd, test_config_yaml)
from .utils import capture_stdouterr


class TestFlags(TestCase):
    """
    This class runs tests of elvers with basic flags:
    -h --help
    -n --dry-run
    """
    @classmethod
    def setUpClass(self):
        """Set up a flags test"""
        pass

    def test_help_flag(self):
        """Test the -h --help flag"""
        command = [elvers_cmd,'-h']
        p_out, p_err = capture_stdouterr(command,here)
        self.assertIn('elvers',p_out)
        self.assertIn('snakemake',p_out)

    def test_dry_run_flag(self):
        """Test the -n --dry-run flag"""
        which_workflow = 'default'
        flags = '-n'
        command = [elvers_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)
        verify_present = '''
Job counts:
	count	jobs
	1	dammit_annotate
	1	elvers
	20	fastqc_pretrim
	20	fastqc_trimmed
	10	http_get_fq1
	10	http_get_fq2
	10	khmer_pe_diginorm
	10	khmer_split_paired
	1	multiqc
	1	rename_trinity_fasta
	1	rename_trinity_gene_trans_map
	1	salmon_index
	2	salmon_quant_combine_units
	1	sourmash_compute_assembly
	10	sourmash_compute_pe_interleaved
	10	trimmomatic_pe
	1	trinity
	110'''
        # Skip the first element b/c empty line
        for line in verify_present.split("\n")[1:]:
            self.assertIn(line,p_out)

        # Check that job descriptions are printed as they are
        # added to the list of dry run tasks
        verify_jobs = '''
--- Quality trimming PE read data with Trimmomatic. ---
--- khmer trimming of low-abundance kmers and digital normalization ---
--- Computing a MinHash signature of the kmer-trimmed reads with Sourmash ---
--- Assembling read data with Trinity ---
--- Indexing the transcriptome with Salmon ---'''
        # Skip the first element b/c empty line
        for line in verify_jobs.split("\n")[1:]:
            self.assertIn(line,p_out)

    @classmethod
    def tearDownClass(self):
        """Tear down after a test"""
        pass

