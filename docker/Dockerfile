FROM debian:jessie

ENV LC_ALL C

RUN apt-get -y update && \
    apt-get -y install bzip2 wget unzip
    #&& rm -rf /var/lib/apt/lists/*

# Install hapCUT
RUN apt-get install -y build-essential zlib1g-dev ncurses-dev
RUN wget -O hapcut.zip https://github.com/vibansal/hapcut/archive/844af08c.zip && \
    unzip hapcut.zip && \
    cd hapcut-*/ && \
    make && \
    mv extractHAIRS HAPCUT /usr/bin/ && \
    rm -rf hapcut*

# Install Miniconda
RUN echo 'export PATH=/opt/miniconda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda/bin:$PATH
RUN conda config --add channels bioconda --add channels r --add channels conda-forge

# Install all dependencies available as conda packages
RUN conda install -y python=3.5.2 snakemake=3.7.1 samtools=1.2 picard=1.141 \
        bedtools=2.23.0 vcftools=0.1.14 bwa=0.7.12 pbsim=1.0.3 whatshap=0.13 \
        gatk=3.5
# desired picard 1.130 is not available
# desired vcftools r940 is not available
# pbsim version unclear, chose most probable one

# Install phASER into a separate env
RUN conda create -y -n phaser python=2.7 intervaltree pyvcf scipy numpy samtools bedtools
RUN wget -O phaser.zip https://github.com/secastel/phaser/archive/6b527927.zip && \
    unzip phaser.zip && mkdir /opt/phaser/ && mv phaser-*/phaser/*.py /opt/phaser/ && rm -r phaser-*/
COPY phaser /usr/bin/

RUN apt-get install -y pigz

WORKDIR /io/
