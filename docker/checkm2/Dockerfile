FROM ubuntu:focal AS builder

RUN apt-get update \
	&& apt-get -y install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG CONDA_ENVIRONMENT_REPO_HASH
WORKDIR /opt
RUN git clone \
	--recurse-submodules \
	https://github.com/PacificBiosciences/pb-metagenomics-tools.git
RUN cd /opt/pb-metagenomics-tools \
	&& git checkout "${CONDA_ENVIRONMENT_REPO_HASH}"


FROM continuumio/miniconda3:4.12.0

MAINTAINER Karen Fang <karen@dnastack.com>

RUN apt-get -qq update \
	&& apt-get -qq install \
		gawk \
		libarchive13 \
	&& rm -rf /var/lib/apt/lists/*

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CONDA_ENVIRONMENT
ENV CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"

COPY --from=builder /opt/pb-metagenomics-tools/HiFi-MAG-Pipeline/scripts /opt/scripts
# Additional or overwritten scripts
# Must copy a file that definitely exists (build.env) to allow optional files to be copied in (scripts/*)
COPY build.env scripts/* /opt/scripts/
RUN rm /opt/scripts/build.env && \
	chmod +x /opt/scripts/*

COPY --from=builder "/opt/pb-metagenomics-tools/HiFi-MAG-Pipeline/envs/${CONDA_ENVIRONMENT}.yml" /opt/envs/
RUN conda install mamba -c conda-forge && \
	mamba env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:/opt/scripts:$PATH

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
