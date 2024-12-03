FROM continuumio/anaconda3:5.2.0

ENV MPLBACKEND=Agg

RUN pip install --no-cache-dir progressbar33 pxl

RUN mkdir /work
WORKDIR /work
