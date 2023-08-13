FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get install -y wget git python3 python3-pip default-jre
RUN apt-get install -y libbz2-dev liblzma-dev
RUN pip install ipython
RUN pip install pandas scipy vcfpy seaborn

RUN mkdir /pipeline
WORKDIR /pipeline
