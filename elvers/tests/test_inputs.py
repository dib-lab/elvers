#!/usr/bin/env python
import os
from unittest import TestCase
from .const import (here, elvers_cmd, test_config_yaml)
from .utils import capture_stdouterr


class TestInputs(TestCase):
    """
    This class runs tests with improperly formatted inputs
    """
    def test_bad_samples_entry(self):
        """Test a configfile with no samples tsv specified"""
        configfile = os.path.realpath(os.path.join(here,'test_files/bad-samples-entry.yaml'))
        command = [elvers_cmd, configfile, '-n']
        p_out, p_err = capture_stdouterr(command,here)
        assert "Error: trying to run `get_data` workflow, but the samples tsv file is not specified in your configfile. Please fix." in p_err

    def test_bad_samples_tsv(self):
        """Test an improperly formatted samples tsv file (spaces instead of tabs)"""
        configfile = os.path.realpath(os.path.join(here,'test_files/bad-samples-tsv.yaml'))
        samples_tsv = os.path.realpath(os.path.join(here,'test_files/bad-samples-tsv.tsv'))
        command = [elvers_cmd, configfile, '-n']
        p_out, p_err = capture_stdouterr(command,here)
        assert f"{samples_tsv} file is not properly formatted. Please fix." in p_err

    def test_bad_inputs_entry(self):
        """Test a configfile with a typo in the inputs section for a program"""
        configfile = os.path.realpath(os.path.join(here,'test_files/bad-inputs-entry.yaml'))
        command = [elvers_cmd, configfile, '-n']
        p_out, p_err = capture_stdouterr(command,here)
        assert f"cannot find corresponding outputs for inputs" in p_err

