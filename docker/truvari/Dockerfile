FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

LABEL maintainer="Zev Kronenberg <zev@pacbio.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME="${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG="${IMAGE_TAG}"

ARG MAFFT_VERSION=7.505-1
ARG TRUVARI_VERSION=5.4.0
ARG PYABPOA_VERSION=1.4.3
ARG PYWFA_VERSION=0.5.1
ARG SETUPTOOLS_VERSION=58.5.3
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_${MAFFT_VERSION}_amd64.deb \
	&& dpkg -i mafft_${MAFFT_VERSION}_amd64.deb \
	&& rm mafft_${MAFFT_VERSION}_amd64.deb
RUN python3 -m pip install --upgrade pip \
	&& python3 -m pip install setuptools==${SETUPTOOLS_VERSION} wheel \
	&& python3 -m pip install --no-build-isolation pyabpoa==${PYABPOA_VERSION} \
	&& python3 -m pip install pywfa==${PYWFA_VERSION} \
	&& python3 -m pip install "truvari[bwa]==${TRUVARI_VERSION}"
