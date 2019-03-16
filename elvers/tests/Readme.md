# elvers tests

This directory contains tests for elvers.

The test strategy is as follows:

* Use the `unittest` library to write unit tests
* Use the `subprocess` library to run elvers workflows from the command line
* Capture the output of each `subprocess` command and make assertions about its contents

To run all tests, use the `pytest` command to automatically search for and run all
tests. This command can be run either from this `tests/` directory or from the top level
of the repository: 

```
pytest -v
```

If all goes well, you should see output like this:

```
$ pytest -v
=============================================== test session starts ================================================
platform darwin -- Python 3.6.3, pytest-4.2.0, py-1.7.0, pluggy-0.8.1 -- /Users/charles/.pyenv/versions/miniconda3-4.3.30/bin/python
cachedir: .pytest_cache
rootdir: /temp/elvers/tests, inifile:
collected 10 items

test_config.py::TestConfig::test_build_config PASSED                                                         [ 10%]
test_config.py::TestConfig::test_extra_config PASSED                                                         [ 20%]
test_dag.py::TestDag::test_dag_flag PASSED                                                                   [ 30%]
test_dag.py::TestDag::test_dagfile_flag PASSED                                                               [ 40%]
test_dag.py::TestDag::test_dagpng_flag PASSED                                                                [ 50%]
test_flags.py::TestFlags::test_dry_run_flag PASSED                                                           [ 60%]
test_flags.py::TestFlags::test_help_flag PASSED                                                              [ 70%]
test_print.py::TestPrint::test_print_params PASSED                                                           [ 80%]
test_print.py::TestPrint::test_print_rules PASSED                                                            [ 90%]
test_print.py::TestPrint::test_print_workflows PASSED                                                        [100%]

============================================ 10 passed in 12.61 seconds ============================================
```

