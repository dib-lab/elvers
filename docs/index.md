# elvers


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


**`elvers`** is an automated workflow system, designed to facilitate running a number of tools on a set of data without the need to babysit each individual program. **`elvers`** is free and open source, and relies solely on programs that are also free, open source, and reasonably installable (which we consider part of being open).

## Why use elvers?

Some parts of bioinformatic analysis are standardized and just need to be run. For others, we may want to compare results between different input parameters or different programs at individual steps. Both of these cases are facilitated by using a workflow system which 1) handles all installations, 2) executes steps in a standard and repeatable manner, 3) keeps track of all steps that have been run, and 4) provides a way to pick up from where you left off, should any issues arise during
execution. **`elvers`** is designed to be easy to begin executing default workflows, but flexible and extensible, allowing full configuration of program parameters and programs within each workflow. 

## What do I need?

To run **`elvers`**, you need to input some data, either **reads, a reference transcriptome, or both**. To run, you'll need access to a linux system via a command-line interface. Most programs also work on Mac, but a few are troublesome. If you'll be running *de novo* assembly, we recomment a high-memory machine.


## What programs are run?

**`elvers`** is a **workflow system** or **workflow playbook**, meaning that it provides a number of different workflows that can be run. Each workflow consists of one or more programs to run an analysis on your data. The workflow(s) you choose will depend on your data and analysis needs. A list of available workflows and tools can be found on the [workflows](workflows.md) page.


## Name disambiguation

**`elvers`**, formerly `eelpond`, is an automated workflow system. It evolved from a snakemake update of the "Eel Pond" protocol for *de novo* RNAseq analysis, previously developed by CT Brown and members of the dib-lab. An automated version of this protocol is the `default` workflow in `elvers`, and the older,  step-by-step Eel Pond protocol can be found [here](https://eel-pond.readthedocs.io/en/latest/).

## Authors

**`elvers`** was developed by N Tessa Pierce with invaluable support and feedback from CT Brown, Charles Reid, Lisa Johnson, Taylor Reiter, Luiz Irber, and the entire dib-lab. 

## Citation

_coming soon_, please contact Tessa at ntpierce@gmail.com if you need to cite and this hasn't been filled in yet!
