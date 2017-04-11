# Phasing comparison experiments

This repository contains the instructions reproducing the experiments described
in our study “WhatsHap: fast and accurate read-based phasing”, for which a
pre-print is available at <http://dx.doi.org/10.1101/085050>.

The workflow is described as a “Snakefile”, which needs to be run with
[snakemake](http://snakemake.bitbucket.org). We also provide a Dockerfile from
which it should be possible to reproduce all results with a single command.
This does not quite work, yet, so for now it merely serves as a compact way of
describing how to install all dependencies.

## Overview

The workflow is designed to be run from the directory in which this README file
is located. In principle, you only need to run `snakemake` here and then the
rest is done automatically, but all dependencies need to be installed first.
Further things to know:

* About 600 GB of hard disk space is needed (better 1 TB to be on the safe side).
* No data is included in this repository since all required files will be
  downloaded.
* About 100 GB of data will be downloaded. Most of this is BAM files from the
  Genome in a Bottle dataset, but also the human genome reference (hg37).
* The pipeline requires proprietary software which we cannot download and install
  automatically, so those tools needs to be installed manually.
* 8 cores should work, but more is better. Presumably 8 GB RAM works, but more
  is better. We tested on a 48-core machine with 256 GB RAM.

## Install non-free dependencies

The workflow uses the statistical-phasing tool [SHAPEIT](http://shapeit.fr) and
the [Genome Analysis Toolkit (GATK)](https://software.broadinstitute.org/gatk/),
both of which cannot be re-distributed due to the way they are licensed. Before
starting the workflow, you need to install them manually.

### SHAPEIT

- Download version v2.r837 of SHAPEIT from the
  [download page](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download).
  Get the file named `shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz`.
- Unpack the tgz file and place the `shapeit` binary that it contains into
  `restricted-software/`.

### GATK

To exactly reproduce the results of our pipeline, you need to use GATK 3.5.
That said, the ReadBackedPhasing component of GATK has not changed in a while,
so probably any version that is not too old will give the same results.

- Go to
  [the download page for “archived versions” of the GATK](https://software.broadinstitute.org/gatk/download/archive)
- Download GATK 3.5-0-g36282e4
- Unpack the downloaded file and place the `GenomeAnalysisTK.jar` file into
  `restricted-software/`.

### Checksums

This step is optional, but if you want to ensure that the files you downloaded
are actualy the ones that we used, go into the `restricted-software` directory
and run

    md5sum -c MD5SUM

For both files, you should see an `OK` message. This check is not run as part
of the workflow because we want to give you the option of using different
versions of the tools.

## Install all other dependencies

1. Install necessary system packages. We assume you use Debian or Ubuntu. If you
    do not have root access on your system, you can try to ignore this step.
    These packages are typically already installed.

        sudo apt-get install bzip2 wget unzip build-essential zlib1g-dev ncurses-dev pigz

2. Install hapCUT

        wget -O hapcut.zip https://github.com/vibansal/hapcut/archive/844af08c.zip
        unzip hapcut.zip
        cd hapcut-*/
        make

    Then copy the `extractHAIRS` and `HAPCUT` binaries to some location that is on
    your `$PATH`.

3. Install hapCUT2

    TODO (very similar to above)

4. Install Miniconda. Most of the remaining dependencies are available as
  packages in the [bioconda](http://bioconda.github.io/) Conda channel, for which
  we need conda.

        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

    Then follow the instructions on screen. After conda is installed, add the
    required channels to the configuration.

        conda config --add channels bioconda --add channels r --add channels conda-forge

5. Install all dependencies available as Conda packages.

        conda install -y python=3.5.2 snakemake=3.7.1 samtools=1.2 picard=1.126 \
            bedtools=2.23.0 vcftools=0.1.14 bwa=0.7.12 pbsim=1.0.3 whatshap=0.13 \
            biopython=1.68 htslib=1.4

6. Install phASER into a separate conda environment.

        conda create -y -n phaser python=2.7 intervaltree pyvcf scipy numpy samtools bedtools
        wget -O phaser-05832c1.zip https://github.com/secastel/phaser/archive/05832c1.zip
        unzip phaser-05832c1.zip

    Then copy `phaser-*/phaser/*.py` to a directory of your choosing.
    TO DO: The description of this step is incomplete. You need to somehow make
    sure that running `phaser` on the command-line works. See also the wrapper
    `docker/phaser`.

## Run the workflow

Start the workflow with

    nice snakemake -p -j 16

Adjust the 16, which is the number of cores you want to use.

## Results

These are the result files:

* `eval/summary.eval` -- Table with evaluation results
* `eval/summary.components` -- Table with sizes of connected components
* `eval/*.pdf` -- Plots
