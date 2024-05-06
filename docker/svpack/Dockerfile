FROM quay.io/pacbio/pb_wdl_base@sha256:963c8ef4cd011cd044d5efcff2759aa37b86d0ca5739ce576f60a0ae0967292c

LABEL maintainer="Billy Rowell <wrowell@pacb.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG SVPACK_GIT_HASH
WORKDIR /opt
RUN git clone https://github.com/PacificBiosciences/svpack.git \
	&& cd svpack \
	&& git checkout "${SVPACK_GIT_HASH}" \
	&& chmod +x /opt/svpack/svpack

ENV PATH "${PATH}:/opt/svpack"
