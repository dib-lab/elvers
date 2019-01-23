# EelPond


```
                           ___
                        .-'   `'.
                       /         \
                      |           ;
                      |           |           ___.--,
             _.._     |O)  ~  (O) |    _.---'`__.-( (_.       
      __.--'`_.. '.__.\      '--. \_.-' ,.--'`     `""`
     ( ,.--'`   ',__ /./;     ;, '.__.'`    __
     _`) )  .---.__.' / |     |\   \__..--""  """--.,_
    `---' .'.''-._.-'`_./    /\ '.  \_.-~~~````~~~-.__`-.__.'
          | |  .' _.-' |    |  \  \  '.
           \ \/ .'     \    \   '. '-._)
            \/ /        \    \    `=.__`-~-.
            / /\         `)   )     / / `"".`\
      , _.-'.'\ \        /   /     (  (   /  /
       `--~`  )  )    .-'  .'       '.'. |  (
             (/`     (   (`           ) ) `-;
              `       '--;            (' 

```
eelpond started as a snakemake update of the Eel Pond Protocal for *de novo* RNAseq analysis. It has evolved slightly to enable a number of workflows for (mostly) RNA data, which can all be run via the `eelpond` workflow wrapper. `eelpond` uses [snakemake](https://snakemake.readthedocs.io) for workflow management and [conda](https://conda.io/docs/) for software installation. The code can be found [here](https://github.com/dib-lab/eelpond). 


## Getting Started

Linux is the recommended OS. Nearly everything also works on MacOSX, but some programs (fastqc, Trinity) are troublesome.

If you don't have conda yet, install [miniconda](https://conda.io/miniconda.html) (for Ubuntu 16.04 [Jetstream image](https://use.jetstream-cloud.org/application/images/107)):
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
echo export PATH="$HOME/miniconda3/bin:$PATH" >> ~/.bash_profile
source ~/.bash_profile
```

Now, get the eelpond code
```
git clone https://github.com/dib-lab/eelpond.git
cd eelpond
```

Create a conda environment with all the dependencies for eelpond
```
conda env create --file environment.yml -n eelpond
```

Activate that environment. You'll need to do this anytime you want to run eelpond
```
conda activate eelpond
```
Now, grab the test data, and untar:
```
curl -L https://osf.io/chb7z/download -o nema_testdata.tar.gz
tar xvf nema_testdata.tar.gz
```

Now you can start running workflows!


## Default workflow: Eel Pond Protocol for *de novo* RNAseq analysis

The Eel Pond protocol (which inspired the `eelpond` name) included line-by-line commands that the user could follow along with using a test dataset provided in the instructions. We have re-implemented the protocol here to enable automated *de novo* transcriptome assembly, annotation, and quick differential expression analysis on a set of short-read Illumina data using a single command. See more about this protocol [here](Eel_Pond_workflow.md).

To test the default workflow:
```
./run_eelpond nema-test eel_pond
```
This will run a small set of _Nematostella vectensis_ test data (from [Tulin et al., 2013](https://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16))

## Running Your Own Data

To run your own data, you'll need to create two files, a `tsv` file containing 
your sample info, and a `yaml` file containing basic configuration info. To start,
copy the test data files so you can modify them.
```
cp nema_samples.tsv <my-tsv-samples.tsv>
```
Modify this tab-separated file with your sample names and locations for each file. 

Notes:  
    - `eelpond` is a bit finicky here: all columns must be separated by tabs, not spaces. If you get an immediate error about your samples file, this is likely the cause.   
    - the `unit` column allows you to enter multiple files or file pairs for each samples, such as files that are from different lanes. If you don't have a `unit` bit of information, just add something short as placeholder (e.g. `a`). All units for a single sample are combined prior to the differential expression steps. At some point, we may enable a `no_units` version of the `eelpond`, but it's not in our immediate plans :). 

Next, build a configfile to edit:
```
./run_eelpond config_name --build_config
```

This configfile will contain all the default paramters for each step of the workflow you target.
If you don't specify any targets, it will default to the full "eel_pond" pipeline, which executes read
preprocessing, assembly, annotation, and quantification. 

Then, modify this configfile as necessary. 

**The configfile must contain at least:**
```
samples: path/to/my-tsv-samples.tsv
```
**which directs `eelpond` to your sample files.**


## Additional Info

Each independent step is split into a smaller workflow that can be run independently, if desired, e.g. `./run_eelpond nema-test preprocess`. Individual tools can also be run independently, see [Advanced Usage](advanced_usage.md).

See the help, here:
```
./run_eelpond -h
```
**available workflows:**  

  - preprocess: Read Quality Trimming and Filtering (fastqc, trimmomatic)
  - kmer_trim: Kmer Trimming and/or Digital Normalization (khmer)
  - assemble: Transcriptome Assembly (trinity)
  - assemblyinput: Specify assembly for downstream steps
  - annotate : Annotate the transcriptome (dammit, sourmash)
  - quantify: Quantify transcripts (salmon) 
  - eel_pond: preprocess, kmer_trim, assemble, annotate, quantify 

You can see the available workflows (and which programs they run) by using the `--print_workflows` flag:
```
./run_eelpond nema-test --print_workflows
```



**References:**  

  * [original eel-pond protocol docs, last updated 2015](https://khmer-protocols.readthedocs.io/en/ctb/mrnaseq/)
  * [eel-pond protocol docs, last updated 2016](http://eel-pond.readthedocs.io/en/latest/)
  * [DIBSI, nonmodel RNAseq workshop, July 2017](http://dibsi-rnaseq.readthedocs.io/en/latest/)
  * [SIO-BUG, nonmodel RNAseq workshop, October 2017](http://rnaseq-workshop-2017.readthedocs.io/en/latest/index.html)


