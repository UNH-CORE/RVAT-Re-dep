FROM continuumio/anaconda3:4.2.0

ENV MPLBACKEND=Agg

RUN pip install \
    --trusted-host pypi.org \
    --trusted-host pypi.python.org \
    --trusted-host files.pythonhosted.org \
    --no-cache-dir \
    progressbar33 \
    seaborn==0.7.0 \
    pxl==0.0.9 \
    --no-deps

RUN mkdir /work
WORKDIR /work
