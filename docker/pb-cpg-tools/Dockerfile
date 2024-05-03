FROM pb_wdl_base:build1

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PB_CPG_TOOLS_VERSION
RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/${PB_CPG_TOOLS_VERSION}/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz && \
	tar --no-same-owner -zxvf pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz --directory /opt && \
	rm pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu.tar.gz
RUN ln -s /opt/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores /usr/local/bin/

ENV PILEUP_MODEL_DIR /opt/pb-CpG-tools-${PB_CPG_TOOLS_VERSION}-x86_64-unknown-linux-gnu/models/
