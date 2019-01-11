FROM continuumio/miniconda:4.5.11
LABEL maintainer="Fadhil Abubaker"

WORKDIR /home

RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge

COPY requirements.txt /home

RUN conda install --yes --file requirements.txt

COPY bin/* /usr/bin/

COPY . /home
