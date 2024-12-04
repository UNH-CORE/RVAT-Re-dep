FROM continuumio/anaconda3:5.2.0

ENV MPLBACKEND=Agg

RUN pip install --no-cache-dir progressbar33 pxl>=0.0.10

RUN mkdir /work
WORKDIR /work
