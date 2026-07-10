FROM quay.io/pacbio/pb_wdl_base@sha256:fa70d211e884a7348b53ff0120eed1f8070f053dd4a1c4b0bfeec74e1a42c5b0

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG IMAGE_DESCRIPTION
ARG IMAGE_VENDOR
ARG IMAGE_URL
ARG IMAGE_MAINTAINER

LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
	org.opencontainers.image.description="${IMAGE_DESCRIPTION}" \
	org.opencontainers.image.authors="${IMAGE_MAINTAINER}" \
	org.opencontainers.image.vendor="${IMAGE_VENDOR}" \
	org.opencontainers.image.url="${IMAGE_URL}" \
	org.opencontainers.image.version="${IMAGE_TAG}"

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

RUN bam2fasta --version \
	&& bam2fastq --version \
	&& ccs-kinetics-bystrandify --version \
	&& extracthifi --version \
	&& pbindex --version \
	&& pbindexdump --version \
	&& pbmerge --version \
	&& zmwfilter --version \
	&& stitchread --version
