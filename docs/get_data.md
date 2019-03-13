# Get Data Utility Rule

To run any workflow, we need to get reads (and any assemblies already generated) into the right format for `elvers`. We can either (1) link the data from another location on your machine, or (2) download the data via http or ftp. This is done through a utility rule called `get_data`.

**At the moment, get_data only works with gzipped files. Checks and improved functionality coming soon** 


## Modifying Params for get_data:

Be sure to set up your sample info and build a configfile first (see [Understanding and Configuring Workflows](configure.md)).

To see the available parameters for the `get_data` utility rule, run
```
./run_elvers config get_data --print_params
```

In here, you'll see a section for "get_data" parameters that looks like this:

```
  ####################  get_data  ####################
get_data:
  download_data: false
  use_ftp: false
  #####################################################
```

If you want to download data, you'll want to change the `download_data` and `use_ftp` parameters appropriately. To download via http, set `download_data: True`. To download via ftp, set `download_data: True` and `use_ftp: True`.


## Output Files

The output of the `get_data` step is all your input data files in a subdirectory (`input_data`) within your output directory (basename_out). These will either be links or downloaded files, depending on the options you specified in the config file.
