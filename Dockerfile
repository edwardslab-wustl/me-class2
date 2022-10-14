FROM continuumio/miniconda3:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && apt-get clean

ADD https://api.github.com/repos/edwardslab-wustl/me-class2/git/refs/heads/main/version.json
RUN /opt/conda/bin/python3 -m pip install git+https://github.com/edwardslab-wustl/me-class2.git

RUN mkdir -p /usr/local/meclass2 
ADD example_data.tgz /usr/local/meclass2
ADD  utils-0.2.0.tar.gz /usr/local/meclass2

ENV PATH=/opt/conda/bin:/usr/local/utils/$PATH

