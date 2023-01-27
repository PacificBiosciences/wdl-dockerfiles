FROM continuumio/miniconda3:4.12.0

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update && \
	apt-get -qq install \
	git

ARG SVPACK_GIT_HASH
ENV SVPACK_GIT_HASH "${SVPACK_GIT_HASH}"

WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/svpack.git
RUN cd svpack && \
	git checkout "${SVPACK_GIT_HASH}"

ENV CONDA_ENVIRONMENT svpack
COPY "${CONDA_ENVIRONMENT}.yaml" /opt/envs/

RUN conda env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yaml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:/opt/svpack:$PATH

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
