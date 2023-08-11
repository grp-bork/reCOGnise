FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="0.1"
LABEL description="This is a Docker Image for the reCOGnise tool."


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
RUN apt clean

RUN mkdir -p /opt/software && \
	cd /opt/software && \
	git clone https://github.com/cschu/reCOGnise.git && \
	cd reCOGnise && \
	pip install .
  
RUN cd /opt/software && \
	wget -q https://github.com/jfmrod/MAPseq/releases/download/v1.2.6/mapseq-1.2.6-linux.tar.gz && \
	tar xzf mapseq-1.2.6-linux.tar.gz && \
	rm mapseq-1.2.6-linux.tar.gz && \
	mv mapseq-1.2.6-linux mapseq && \
	ln -s /opt/software/mapseq/mapseq /usr/bin/ && \
	ln -s /opt/software/mapseq/share /usr/bin/

RUN cd /opt/software && \
	git clone https://github.com/motu-tool/fetchMGs.git

CMD ["recognise"]
