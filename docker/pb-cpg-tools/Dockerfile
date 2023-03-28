FROM ubuntu:focal

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update && \
	apt-get -qq install \
	wget \
	build-essential \
	&& rm -rf /var/lib/apt/lists/*

ARG PB_CPG_TOOLS_VERSION
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/${PB_CPG_TOOLS_VERSION}/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar -zxvf pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /usr/local/bin/

ENV PILEUP_MODEL_DIR /opt/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/models/
