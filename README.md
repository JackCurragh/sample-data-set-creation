# Sample Dataset Creation
Contains methods used to create test data (primarily for nf-core modules)

In my experience test datasets that are made available via projects such as [nf-core](https://github.com/nf-core/test-datasets) are not suitable for use with Ribo-Seq applications. Although the correct file types are available there are still a number of issues. For example, short read RNA-Seq files do not show 3nt periodicity that man Ribo-Seq tools look for. 

In this repo each subdirectory is it's own project for the creation of a specific type of test data file. 

## Required 
- conda
- docker
- jupyter notebooks
