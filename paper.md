---
title: 'elvers: an automated, extensible bioinformatics workflow system'
tags:
  - snakemake 
  - RNA-seq
  - Python
authors:
 - name: N Tessa Pierce
   orcid: 0000-0002-2942-5331 
   affiliation: University of California, Davis
 - name: Charles Reid
   orcid:
   affiliation: University of California, Davis
   
 - name: C. Titus Brown
   orcid: 0000-0001-6001-2677
   affiliation: University of California, Davis
date: 22 Mar 2019
bibliography: paper.bib
---

# Summary

elvers is a bioinformatics workflow system that leverages snakemake 
[@10.1093/bioinformatics/bts480] to facilitate automated and reproducible analyses.

elvers provides an extensible playbook of snakemake rules and worklows that 
can be mixed and matched for completely customizable automated analyses.
The elvers command-line interface enables full user customization of all
included programs, while handling specification of any intermediate files.

elvers requires users to input either data files (via spreadsheet), reference
transcriptomes or genomes (via YAML configuration file), or both.

# References
