# Get Data Utility Rule

To run any workflow, we need to get reads (and any assemblies already generated) into the right format for `eelpond`. We can either (1) link the data from another location on your machine, or (2) download the data via http or ftp. This is done through a utility rule called `get_data`.

First, you need to specify your data (file path if on your machine, or download link if remote) via a *properly formatted* tsv file. For example, our test data looks like this:


If we had our test data within the `nema_testdata` folder in the main `eelpond` directory:
    ```
    sample  unit    fq1 fq2 condition
    0Hour   001 nema_testdata/0Hour_ATCACG_L002_R1_001.extract.fastq.gz nema_testdata/0Hour_ATCACG_L002_R2_001.extract.fastq.gz time0
    0Hour   002 nema_testdata/0Hour_ATCACG_L002_R1_002.extract.fastq.gz nema_testdata/0Hour_ATCACG_L002_R2_002.extract.fastq.gz time0
    6Hour   001 nema_testdata/6Hour_CGATGT_L002_R1_001.extract.fastq.gz nema_testdata/6Hour_CGATGT_L002_R2_001.extract.fastq.gz time6
    6Hour   002 nema_testdata/6Hour_CGATGT_L002_R1_002.extract.fastq.gz nema_testdata/6Hour_CGATGT_L002_R2_002.extract.fastq.gz time6
    ```

If we want to download the test data instead:
    ```
    sample  unit    fq1 fq2 condition
    0Hour   001 https://osf.io/vw4dt/download   https://osf.io/b47s2/download   time0
    0Hour   002 https://osf.io/92jr6/download   https://osf.io/qzea8/download   time0
    6Hour   001 https://osf.io/ft3am/download   https://osf.io/jqmsx/download   time6
    6Hour   002 https://osf.io/rnkq3/download   https://osf.io/bt9vh/download   time6
    ```

**Note that proper formatting means all columns must be separated by tabs, not spaces!**

