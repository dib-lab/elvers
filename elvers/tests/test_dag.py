from unittest import TestCase
from .const import (here, elvers_cmd, test_config_yaml)
from .utils import capture_stdouterr
import subprocess
import os


class TestDag(TestCase):
    """
    This class runs tests of elvers dag flags and functionality:
    --dag
    --dagfile
    --dagpng
    """
    def test_dag_flag(self):
        which_workflow = 'default'
        flags = '--dag'
        command = [elvers_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)

        # Look for this message in snakemake output
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
        command = [elvers_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)

        # Look for this message in snakemake output
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
        pngfile_name = 'dag.png'
        pngfile_fullpath = os.path.join(here,pngfile_name)
        which_workflow = 'default'
        flags = '--dagpng=%s'%(pngfile_fullpath)
        command = [elvers_cmd, test_config_yaml, which_workflow, flags]
        p_out, p_err = capture_stdouterr(command,here)

        # Look for this message in snakemake output
        self.assertIn('Building DAG of jobs...',p_err)

        # Output should print where pngfile went
        self.assertIn('Printed workflow dag to png file',p_out)
        self.assertIn(pngfile_fullpath, p_out)

        # The pngfile should now be a file on disk
        self.assertTrue(os.path.exists(pngfile_fullpath))
        self.assertTrue(os.path.isfile(pngfile_fullpath))

        # Clean up
        subprocess.call(['rm','-f',pngfile_fullpath])

