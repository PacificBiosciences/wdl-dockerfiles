FROM mambaorg/micromamba:2.5.0
LABEL maintainer="Zev Kronenberg <zkronenberg@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG SEVERUS_VERSION
ARG HTSLIB_VERSION
ARG BCFTOOLS_VERSION
ARG SAMTOOLS_VERSION
RUN micromamba create --yes \
	--name severus \
	--channel defaults \
	--channel bioconda \
	--channel conda-forge \
	severus=${SEVERUS_VERSION} \
	wget bcftools=${BCFTOOLS_VERSION} \
	samtools=${SAMTOOLS_VERSION} \
	htslib=${HTSLIB_VERSION} \
	&& micromamba env export --name severus --explicit \
	&& micromamba clean --all --yes
