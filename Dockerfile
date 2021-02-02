FROM continuumio/miniconda:latest

RUN apt-get update && \
    apt-get -y install gcc mono-mcs libz-dev && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/usr/bin:$PATH"
ENV PATH="/opt/conda/bin:$PATH"
RUN conda config --add channels defaults
RUN conda config --add channels anaconda
RUN conda config --add channels bioconda
RUN conda config --set always_yes yes --set changeps1 no
RUN conda install -c bioconda python=3.6 cython numpy networkx seaborn pyBigwig six pysam \
    ujson pytest scipy matplotlib samtools future pytest-cov codecov

COPY . /opt/sequencing_tools

RUN cd /opt/sequencing_tools; pip install .
ENV PATH="/opt/conda/bin:${PATH}"
CMD ["/opt/conda/bin/seqtools"]
