FROM continuumio/miniconda:4.5.11

WORKDIR /home

RUN conda config --add channels bioconda

COPY requirements.txt /home

RUN conda install --yes --file requirements.txt

CMD ["/bin/bash"]
