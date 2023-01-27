FROM continuumio/miniconda3:4.12.0

MAINTAINER Heather Ward <heather@dnastack.com>

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

RUN apt-get -qq update && \
	apt-get -qq install \
	git

ARG PB_CPG_TOOLS_GIT_HASH
ENV PB_CPG_TOOLS_GIT_HASH "${PB_CPG_TOOLS_GIT_HASH}"

WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/pb-CpG-tools.git
RUN cd pb-CpG-tools && \
	git checkout "${PB_CPG_TOOLS_GIT_HASH}"

ENV CONDA_ENVIRONMENT pb-cpg-tools

RUN conda env create --file "/opt/pb-CpG-tools/conda_env_cpg.yaml" --name "${CONDA_ENVIRONMENT}"
ENV PATH /opt/conda/envs/${CONDA_ENVIRONMENT}/bin:/opt/pb-CpG-tools:$PATH

RUN chmod +x /opt/pb-CpG-tools/aligned_bam_to_cpg_scores.py

ENV PILEUP_MODEL_DIR /opt/pb-CpG-tools/pileup_calling_model

RUN echo "source activate ${CONDA_ENVIRONMENT}" >> ~/.bashrc
