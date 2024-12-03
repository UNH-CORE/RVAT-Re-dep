FROM continuumio/anaconda3:5.2.0

# Add debian archive for packages
# See https://unix.stackexchange.com/a/743863
RUN echo "deb http://archive.debian.org/debian stretch main contrib non-free" > /etc/apt/sources.list

RUN apt-get update && apt-get install -y texlive-full

ENV MPLBACKEND=Agg

RUN pip install --no-cache-dir progressbar33 pxl

RUN mkdir /work
WORKDIR /work
