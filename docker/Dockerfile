FROM debian:jessie

ENV LC_ALL C

RUN apt-get -y update && \
    apt-get -y install bzip2 wget unzip git pigz
    #&& rm -rf /var/lib/apt/lists/*

# Install hapCUT
RUN apt-get install -y build-essential zlib1g-dev ncurses-dev
RUN wget -O hapcut.zip https://github.com/vibansal/hapcut/archive/844af08c.zip && \
    unzip hapcut.zip && \
    cd hapcut-*/ && \
    make && \
    mv extractHAIRS HAPCUT /usr/bin/ && \
    rm -rf hapcut*

# Install hapCUT2
RUN git clone https://github.com/pjedge/hapcut2.git && \
    cd hapcut2 && \
    git reset --hard d96542fa839e11dc319e85d228dcec444e563b61 && \
    make && \
    mv build/extractHAIRS /usr/bin/extractHAIRS2 && \
    mv build/HAPCUT2 /usr/bin && \
    rm -rf hapcut2

# Install Miniconda
RUN echo 'export PATH=/opt/miniconda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH /opt/miniconda/bin:$PATH
RUN conda config --add channels bioconda --add channels r --add channels conda-forge

# Install all dependencies available as conda packages
RUN conda install -y python=3.5.2 snakemake=3.7.1 samtools=1.2 picard=1.126 \
        bedtools=2.23.0 vcftools=0.1.14 bwa=0.7.12 pbsim=1.0.3 whatshap=0.13 \
        biopython=1.68 htslib=1.4 seaborn=0.7.1 pandas=0.19.2 matplotlib=2.0.0

# Install phASER into a separate env
RUN conda create -y -n phaser python=2.7 intervaltree pyvcf scipy numpy samtools bedtools Cython
RUN wget -O phaser.zip https://github.com/secastel/phaser/archive/05832c1.zip && \
    unzip phaser.zip && \
    cd phaser-*/phaser && \
    mkdir /opt/phaser/ && \
    /opt/miniconda/envs/phaser/bin/python setup.py build_ext -i && \
    mv *.py *.so /opt/phaser/ && \
    cd ../.. && \
    rm -r phaser-*/

COPY phaser /usr/bin/

RUN conda install -y Cython && conda remove whatshap && pip install https://bitbucket.org/whatshap/whatshap/get/459899e25e93.zip

WORKDIR /io/
