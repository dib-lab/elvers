from unittest import TestCase
from .const import (here, run_eelpond_cmd, test_config_yaml)
from .utils import capture_stdouterr
import subprocess
import os


class TestDag(TestCase):
    """
    This class runs tests of eelpond dag flags and functionality:
    --dag
    --dagfile
    --dagpng
    """
    def test_dag_flag(self):
        which_workflow = 'default'
        flags = '--dag'
        command = [run_eelpond_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)
        self.assertIn('Building DAG of jobs...',p_err)

        # Verify some Graphviz notation showed up
        self.assertIn('digraph',p_out)
        self.assertIn('node',p_out)
        self.assertIn('edge',p_out)

    def test_dagfile_flag(self):
        """Test the --dagfile=<dotfile> flag
        """
        dotfile_name = 'dag.dot'
        dotfile_fullpath = os.path.join(here,dotfile_name)
        which_workflow = 'default'
        flags = '--dagfile=%s'%(dotfile_fullpath)
        command = [run_eelpond_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)
        self.assertIn('Building DAG of jobs...',p_err)

        # Output should print where dotfile went
        self.assertIn('Printed workflow dag to dot file',p_out)
        self.assertIn(dotfile_fullpath, p_out)

        # The dotfile should now be a file on disk
        self.assertTrue(os.path.exists(dotfile_fullpath))
        self.assertTrue(os.path.isfile(dotfile_fullpath))

        # Clean up
        subprocess.call(['rm','-f',dotfile_fullpath])

    def test_dagpng_flag(self):
        """Test the --dagpng=<pngfile> flag
        """
        pass

