FROM quay.io/pacbio/pb_wdl_base@sha256:4b889a1f21a6a7fecf18820613cf610103966a93218de772caba126ab70a8e87

COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

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

ARG WAKHAN_VERSION
RUN uv pip install --system git+https://github.com/KolmogorovLab/Wakhan@${WAKHAN_VERSION} \
	pysam pyfaidx "numpy<1.26.0" pandas plotly \
  "scikit-learn==1.2.0" "scipy==1.9.2" ruptures vcf_parser \
  kaleido intervaltree matplotlib
