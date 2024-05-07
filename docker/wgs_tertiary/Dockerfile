FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG ZENODO_RECORD
ARG WORKFLOW_VERSION
# # TODO: For the next zenodo release, package these files separately from the reference files and update this.
# RUN curl -s https://zenodo.org/records/${ZENODO_RECORD}/files/wdl-humanwgs.v${WORKFLOW_VERSION}.resource.tgz \
# 	| tar zxv static_resources/genes static_resources/hpo \
# 	&& mv static_resources/* /opt/

COPY scripts/* /opt/scripts/
ENV PATH "${PATH}:/opt/scripts/"

COPY data/hpo/* /opt/data/hpo/
ENV HPO_TERMS_TSV "/opt/data/hpo/hpoTerms.txt"
ENV HPO_DAG_TSV "/opt/data/hpo/hpoDag.txt"
ENV ENSEMBL_TO_HPO_TSV "/opt/data/hpo/ensembl.hpoPhenotype.tsv"
COPY data/genes/* /opt/data/genes/
ENV ENSEMBL_TO_HGNC "/opt/data/genes/ensembl.hgncSymbol.tsv"
