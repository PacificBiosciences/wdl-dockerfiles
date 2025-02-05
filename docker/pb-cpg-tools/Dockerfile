FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PB_CPG_TOOLS_VERSION
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v${PB_CPG_TOOLS_VERSION}/pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/pb-CpG-tools-v${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /usr/local/bin/
