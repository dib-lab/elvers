#!/usr/bin/env python
from unittest import TestCase
from .const import (here, elvers_cmd, test_config_yaml)
from .utils import capture_stdouterr


class TestConfig(TestCase):
    """
    This class runs tests of elvers config utilities:
    --build_config
    --extra-config
    """
    def test_build_config(self):
        """Test the --build_config flag"""
        command = [elvers_cmd, '--build_config']
        p_out, p_err = capture_stdouterr(command,here)
        pass

    def test_extra_config(self):
        """Test the --build_config flag"""
        command = [elvers_cmd, '-n', '--extra_config']
        p_out, p_err = capture_stdouterr(command,here)
        pass


