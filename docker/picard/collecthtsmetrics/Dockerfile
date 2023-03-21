FROM continuumio/miniconda3:4.12.0

MAINTAINER Heather Ward <heather@dnastack.com>

RUN apt-get -qq update \
	&& apt-get -qq install \
		gawk \
	&& rm -rf /var/lib/apt/lists/*

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG BASE_TOOL
ARG CONDA_ENVIRONMENT
ARG SNAKEMAKE_WRAPPER

COPY "${CONDA_ENVIRONMENT}.yaml" /opt/envs/${CONDA_ENVIRONMENT}.yaml
RUN conda env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yaml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:$PATH

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
