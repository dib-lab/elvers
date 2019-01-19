# Protein Assembly Workflow

*in progress*

The protein assembly workflow relies on the PLASS assembler and some downstream protein mapping tools.

It consists of:  

  - [PLASS](plass.md)
  - [pear](pear.md)
  - [paladin](paladin.md)



## Running Test Data

```
./run_eelpond nema-test.yaml plass_assemble paladin_map
```
This will run a small set of Nematostella vectensis test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16)).


## Running Your Own Data

To run your own data, you'll need to create two files, a tsv file containing your sample info, and a yaml file containing basic configuration info.

IMPORTANT: The sample info must be in a **properly formatted** tsv file. The easiest way to do this is to copy the test data tsv and modify:
```
cp nema_samles.tsv my-samples.tsv
```
Now modify  `my-samples.tsv` with your sample information.

Next, build a configfile to edit:
```
./run_eelpond my-config.yaml plass_assemble paladin_map--build_config
```
This configfile will contain all the default paramters for each step of the targeted workflows. 

Please see the documentation file for each individual program (linked above) for what paramters to modify.
