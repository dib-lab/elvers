# Annotating de novo transcriptomes with dammit

## Quickstart: Running dammit via eelpond

Test data:
```
./run_eelpond nema-test annotate
```
Note that you either need to 1) have already run an assembly, such that a `fasta` file is sitting in the `eelpond/assembly` directory, 2) Run an assembly at the same time, or 3) pass an assembly in via `assemblyinput`

If you have not already run `./run_eelpond nema-test assemble`:

   2) Run trinity assembly at the same time:
   ```
   ./run_eelpond nema-test assemble annotate
   ```
   3) OR, Pass an assembly in via `assemblyinput`
   ```
   ./run_eelpond assemblyinput annotate
   ```
   with an assembly in your `yaml` configfile, e.g.:
   ```
   assemblyinput:
     assembly: rna_testdata/nema.fasta
     gene_trans_map:  rna_testdata/nema.fasta.gene_trans_map
     assembly_extension: '_input'
     ```
    This is commented out in the test data yaml, but go ahead and uncomment (remove leading `#`) in order to use this option

## dammit!

[dammit](http://www.camillescott.org/dammit/index.html) is an annotation pipeline written by [Camille Scott](http://www.camillescott.org/). dammit runs a relatively standard annotation protocol for transcriptomes: it begins by building gene models with [Transdecoder](http://transdecoder.github.io/), and then uses the following protein databases as evidence for annotation:
  -  [Swiss-Prot](https://www.ebi.ac.uk/uniprot) (manually reviewed and curated)
  -  [Pfam-A](http://pfam.xfam.org/)
  -  [Rfam](http://rfam.xfam.org/)
  -  [OrthoDB](http://www.orthodb.org/)
   - [uniref90](http://www.uniprot.org/help/uniref) (uniref is optional with `--full`).
   - nr (nr is optional with `--nr`).

If a protein dataset is available, this can also be supplied to the
`dammit` pipeline with `--user-databases` as optional evidence for
annotation.

In addition, [BUSCO](http://busco.ezlab.org/) v3 is run, which will compare the gene content in your transcriptome with a lineage-specific data set. The output is a proportion of your transcriptome that matches with the data set, which can be used as an estimate of the completeness of your transcriptome based on evolutionary expectation ([Simho et al.2015](http://bioinformatics.oxfordjournals.org/content/31/19/3210.full)). There are several lineage-specific datasets available from the authors of BUSCO. We will use the `metazoa` dataset for this transcriptome.


## Computational Requirements
For the standard pipeline, dammit needs ~18GB of storage space to store its prepared databases, plus a few hundred MB per BUSCO database. For the standard annotation pipeline, we recommend at least 16GB of RAM. This can be reduced by editing LAST parameters via a custom configuration file.

The full pipeline, which uses uniref90, needs several hundred GB of space and considerable RAM to prepare the databases. You'll also want either a fat internet connection or a big cup of patience to download uniref.

For some species, we have found that the amount of RAM required can be proportional to the size of the transcriptome being annotated.

While dammit runs, it will print out which tasks its running to the terminal. dammit is written with a library called [pydoit](http://www.pydoit.org), which is a python workflow library similar to GNU Make and Snakemake. This not only helps organize the underlying workflow, but also means that if we interrupt it, it should properly resume! Caveat: if your job dies, without properly shutting down, snakemake will leave your directory "locked" (a safety feature to prevent two runs/programs from editing the same file simultaneously). If this happens, you'll need to run `eelpond` with the `--unlock` flag.

## Dammit Commands

dammit has two major subcommands: `dammit databases` and `dammit annotate`. `databases` checks that the databases are installed and prepared, and if run with the `--install` flag, will perform that installation and preparation. `annotate` then performs the annotation using installed tools and databases. These commands are automated and integrated into a single rule in `eelpond`. At some point, these may be split into `databases` and `annotate` if that functionality is desired, but it's not in
our immediate plans. 

By default, databases are placed at `eelpond/databases`.

Both `databases` and `annotate` have a `--quick` option, that only installs or runs a "quick" version of the pipeline: transdecoder

On the command line, the commands we would run are as follows:

### Databases
```
#install databases
dammit databases --install --busco-group metazoa --busco-group eukaryota
```
Note: if you have limited space on your instance, you can also install these databases in a different location (e.g. on an external volume) by adding `--database-dir /path/to/databases`. 

You can also use a custom protein database for your species. If your critter is a non-model organism, you will likely need to create your own with proteins from closely-related species. This will rely on your knowledge of your system!

### Annotation

```
dammit annotate trinity.nema.fasta --busco-group metazoa --user-databases <your-database-here> --n_threads 6
```
If you want to run a quick version of the pipeline, add a parameter, `--quick`, to omit OrthoDB, Uniref, Pfam, and Rfam. A "full" run will take longer to install and run, but you'll have access to the full annotation pipeline.

## Customizing Dammit parameters

To modify any program params, you need to add a couple lines to the config file you provide to `eelpond`.

To get a dammit configfile you can modify, run:
```
./run_eelpond dammit.yaml dammit --build_config
```
The output should be a small `yaml` configfile that will look something like this:
```
  ####################  dammit  ####################
dammit:
  busco_group:     # specify all busco groups below here
  - metazoa
  - eukaryota
  db_dir: databases   # specify location for databases (or previously installed databases)
  db_install_only: False   # just install databases, don't run annotation
  db_extra:  
  annot_extra: ' --quick '
```

In addition to changing parameters we've specifically enabled, you can modify the `extra` param to pass any extra parameters.

In dammit, both `databases` and `annotation` take an `extra` param:
```
  db_extra: '--someflag someparam --someotherflag thatotherparam'
  annot_extra: '--someflag someparam --someotherflag thatotherparam'
```
Be sure the modified lines go into the config file you're using to run `eelpond`. Here, we just generated params for `dammit`, but if you're running a larger workflow, we recommend that you generate all params for your workflow in a single file, e.g. `./run_eelpond my-workflow.yaml full --build_config` and edit parameters there.

## Dammit Output

After a successful run, you'll have a new directory called `BASENAME.fasta.dammit`. If you look inside, you'll see a lot of files. For example, for a transcriptome with basename `trinity.nema`, the folder `trinity.nema.fasta.dammit` should contain the following files after a standard (not `--quick`) run:

```
ls trinity.nema.fasta.dammit/
```    
```    
    annotate.doit.db                              trinity.nema.fasta.dammit.namemap.csv  trinity.nema.fasta.transdecoder.pep
    dammit.log                                    trinity.nema.fasta.dammit.stats.json   trinity.nema.fasta.x.nema.reference.prot.faa.crbl.csv
    run_trinity.nema.fasta.metazoa.busco.results  trinity.nema.fasta.transdecoder.bed    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.gff3
    tmp                                           trinity.nema.fasta.transdecoder.cds    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.csv
    trinity.nema.fasta                            trinity.nema.fasta.transdecoder_dir    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.plot.pdf
    trinity.nema.fasta.dammit.fasta               trinity.nema.fasta.transdecoder.gff3
    trinity.nema.fasta.dammit.gff3                trinity.nema.fasta.transdecoder.mRNA
```

As part of eelpond, we copy the two most important files, `trinity.nema.fasta.dammit.fasta` and `trinity.nema.fasta.dammit.gff3` into the main `annotation` directory. `trinity.nema.fasta.dammit.stats.json` also gives summary stats that are quite useful.

If the above `dammit` command is run again, there will be a message:
`**Pipeline is already completed!**`

If you'd like to rerun the dammit pipeline, you'll need to use the `--forceall` flag, like so:
```
./run_eelpond nema-test annotation --forceall
```

## Additional Notes (non-eelpond): Parsing Dammit GFF3 files

Camille wrote dammit in Python, which includes a library to parse gff3 dammit output. If you want to work with this gff3 downstream, use htis parsing library: 

First, enter in a dammit environment. you can find the one eelpond uses, but it might also be easier to just make a new one using the dammit environment file:
```
# from within the eelpond directory
conda env create -n dammit --file rules/dammit/dammit-env.yaml
source activate dammit
```
Remember you can exit your conda environments with `source deactivate`

Then:
```
cd trinity.nema.fasta.dammit
python
```

Use gff3 librarys to output a list of gene IDs:

```
import pandas as pd
from dammit.fileio.gff3 import GFF3Parser
gff_file = "trinity.nema.fasta.dammit.gff3"
annotations = GFF3Parser(filename=gff_file).read()
names = annotations.sort_values(by=['seqid', 'score'], ascending=True).query('score < 1e-05').drop_duplicates(subset='seqid')[['seqid', 'Name']]
new_file = names.dropna(axis=0,how='all')
new_file.head()
new_file.to_csv("nema_gene_name_id.csv")
exit()
```
This will output a table of genes with 'seqid' and 'Name' in a .csv file: `nema_gene_name_id.csv`. Let's take a look at that file:

```
less nema_gene_name_id.csv
```

Notice there are multiple transcripts per gene model prediction. This `.csv` file can be used in `tximport` in downstream DE analysis.

