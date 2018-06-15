# Notes for Advanced Users

Some useful conda, snakemake, workflow hints:

# optional: to make conda installs simpler, set up conda configuration
    conda config --set always_yes yes --set changeps1 no
    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels bioconda

# If you need to modify a conda package:

 you'll need to work with a local install of that package. 
 Here's how to use conda to install the dependenciesfrom the conda recipe.

 install conda-build  

    ```
    conda install conda-build
    ```

 clone the repo of interest and cd into it  
     
     ```
     git clone dammit-repo
     cd dammit-repo
     ```

 There should be a folder called recipe. 
 Use conda-build to build it.  
    
    ``` 
    conda build recipe
    ```

 Install the local code  
    
    ```
    conda install dammit --use-local
    # or, you can use pip:
    # pip install -e . —no-deps
    ```
