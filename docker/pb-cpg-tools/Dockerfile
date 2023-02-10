FROM python:3.9-buster

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update \
	&& apt-get -qq install \
		wget \
		git \
	&& rm -rf /var/lib/apt/lists/*

COPY requirements.txt /opt/
RUN python3 -m pip install -r /opt/requirements.txt

ARG PB_CPG_TOOLS_GIT_HASH
WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/pb-CpG-tools.git \
	&& cd pb-CpG-tools \
	&& git checkout "${PB_CPG_TOOLS_GIT_HASH}"

RUN chmod +x /opt/pb-CpG-tools/aligned_bam_to_cpg_scores.py \
	&& ln -s /opt/pb-CpG-tools/aligned_bam_to_cpg_scores.py /usr/local/bin/

ENV PILEUP_MODEL_DIR /opt/pb-CpG-tools/pileup_calling_model
