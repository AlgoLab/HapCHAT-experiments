# Phasing comparison experiments

This repository contains the instructions reproducing the experiments described
in our study “WhatsHap: fast and accurate read-based phasing”, for which a
pre-print is available at <http://dx.doi.org/10.1101/085050>.

The workflow is described as a “Snakefile” for use with
[snakemake](http://snakemake.bitbucket.org). We also provide a Dockerfile from
which it should be possible to reproduce all results with a single command.

## Overview

The workflow is designed to be run from the directory in which this README file
is located. Downloaded data and the results of computations are stored in
subdirectories.

Further things to know:

* About 600 GB of hard disk space is needed (better 1 TB to be on the safe side).
* No data is included in this repository since all required files will be
  downloaded.
* About 30 GB of data will be downloaded. Most of this are BAM files with
  PacBio and Oxford Nanopore reads, but also the GRCh38 human genome reference.
* The pipeline requires proprietary software which we cannot download and
  install automatically, so this software will need to be installed manually
  (instructions below).
* 8 cores should work, but more is better. Presumably 8 GB RAM works, but more
  is better. The smallest configuration we tested was a 16-core machine with
  16 GB RAM.


## Running the pipeline with Docker

The easiest way to run the pipeline is by using Docker. These are the steps:

1. Install docker (not described here).
2. Install GATK. See below.
3. Build the Docker image:

        sudo docker build -t whatshap-experiments docker/

4. Run the pipeline:

        sudo docker run -it -v $PWD:/io/ whatshap-experiments snakemake -p -j


## Installing GATK

The workflow uses the ReadBackedPhasing part of the
[Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/),
which we cannot re-distribute due to its license. Before starting the workflow,
you need to install it manually.

To exactly reproduce the results of our pipeline, you need to use GATK 3.5.
That said, the ReadBackedPhasing component of GATK has not changed in a while,
so probably any version that is not too old will give the same results.

- Go to
  [the download page for “archived versions” of the GATK](https://software.broadinstitute.org/gatk/download/archive)
- Download GATK 3.5-0-g36282e4
- Unpack the downloaded file and place the `GenomeAnalysisTK.jar` file into
  `restricted-software/`.
- Optionally, run `( cd restricted-software && md5sum -c MD5SUM )` to check that
  the binary is the correct one. (You want to see an output that says `OK`.)


## Running the pipeline without Docker

Without Docker, you will need to install all dependencies yourself. Please
inspect the file `docker/Dockerfile` for this, which gives the commands that are
needed: Everything after a `RUN` directive is a normal shell command.

With the dependencies in place, run snakemake to start the workflow:

    snakemake -p -j 16

Adjust the 16, which is the number of cores you want to use.


## Results

These are the result files:

* `eval/summary.eval` -- Table with evaluation results
* `eval/summary.components` -- Table with sizes of connected components
* `eval/*.pdf` -- Plots
