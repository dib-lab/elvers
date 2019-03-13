#!/usr/bin/env python
from unittest import TestCase
from .const import (here, elvers_cmd, test_config_yaml)
from .utils import capture_stdouterr


class TestPrint(TestCase):
    """
    This class runs tests of elvers print functionality:
    -w --print_workflows
    -r --print_rules
    -p --print_params
    """
    def test_print_workflows(self):
        """Test the --print_workflows flag"""
        command = [elvers_cmd, '--print-workflows']
        p_out, p_err = capture_stdouterr(command,here)

        gold_output = '''
  ####################  Available Eelpond Workflows  ####################


  default:
	get_data
	trimmomatic
	khmer
	trinity
	fastqc
	dammit
	salmon
	sourmash

  protein_assembly:
	get_data
	trimmomatic
	khmer
	plass
	pear
	fastqc
	paladin
	sourmash

  input_data:
	get_data

  preprocess:
	get_data
	fastqc
	trimmomatic

  kmer_trim:
	get_data
	trimmomatic
	khmer

  assemble:
	get_data
	trimmomatic
	khmer
	trinity

  assemblyinput:
	assemblyinput

  annotate:
	dammit

  quantify:
	get_data
	trimmomatic
	salmon

  diffexp:
	get_data
	trimmomatic
	salmon
	deseq2

  sourmash_compute:
	get_data
	trimmomatic
	khmer
	sourmash

  plass_assemble:
	get_data
	trimmomatic
	khmer
	plass

  paladin_map:
	get_data
	trimmomatic
	khmer
	plass
	pear
	paladin

  correct_reads:
	get_data
	trimmomatic
	rcorrector'''

        #self.assertIn(gold_output,p_out)
        #self.assertIn(gold_output,p_err)
        pass

    def test_print_rules(self):
        """Test the --print_rules flag"""
        command = [elvers_cmd, '--print-rules']
        p_out, p_err = capture_stdouterr(command,here)
        pass

    def test_print_params(self):
        """Test the --print_params flag"""
        command = [elvers_cmd, '--print-params']
        p_out, p_err = capture_stdouterr(command,here)
        pass

