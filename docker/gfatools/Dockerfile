FROM continuumio/miniconda3:4.12.0

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"
ARG CONDA_ENVIRONMENT
ENV CONDA_ENVIRONMENT "${CONDA_ENVIRONMENT}"
ARG CALN50_GIT_HASH
ENV CALN50_GIT_HASH "${CALN50_GIT_HASH}"

RUN apt-get -qq update && \
	apt-get -qq install \
	git

WORKDIR /opt
RUN git clone https://github.com/lh3/calN50.git
RUN cd calN50 && \
	git checkout "${CALN50_GIT_HASH}"

COPY "${CONDA_ENVIRONMENT}.populated.yaml" "/opt/envs/${CONDA_ENVIRONMENT}.yaml"

RUN conda env create --file "/opt/envs/${CONDA_ENVIRONMENT}.yaml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:$PATH

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
