FROM ubuntu:focal AS builder

RUN apt-get update \
	&& apt-get -y install \
		git \
	&& rm -rf /var/lib/apt/lists/*

ARG CONDA_ENVIRONMENT_REPO_HASH
WORKDIR /opt
RUN git clone \
	--recurse-submodules \
	https://github.com/PacificBiosciences/HiFiTargetEnrichment.git
RUN cd /opt/HiFiTargetEnrichment \
	&& git checkout "${CONDA_ENVIRONMENT_REPO_HASH}"


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
ARG CONDA_ENVIRONMENT
ENV CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"

COPY --from=builder "/opt/HiFiTargetEnrichment/rules/envs/${CONDA_ENVIRONMENT}.yaml" /opt/envs/
RUN conda env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yaml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:/opt/scripts:$PATH

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc

# workaround for no haplotypes found
RUN sed -i \
	's~^        self.haplotypes = None~        self.haplotypes = []~' \
	/opt/conda/envs/${CONDA_ENVIRONMENT}/lib/python3.10/site-packages/${CONDA_ENVIRONMENT}/cyp2d6_typer.py