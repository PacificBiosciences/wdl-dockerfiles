FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG PBTK_VERSION
RUN wget https://github.com/PacificBiosciences/pbtk/releases/download/v${PBTK_VERSION}/pbtk.tar.gz \
	&& tar --no-same-owner -zxvf pbtk.tar.gz --directory /usr/local/bin \
	&& rm pbtk.tar.gz \
	&& chmod +x /usr/local/bin/bam2fasta \
		/usr/local/bin/bam2fastq \
		/usr/local/bin/ccs-kinetics-bystrandify \
		/usr/local/bin/extracthifi \
		/usr/local/bin/pbindex \
		/usr/local/bin/pbindexdump \
		/usr/local/bin/pbmerge \
		/usr/local/bin/zmwfilter \
		/usr/local/bin/stitchread
