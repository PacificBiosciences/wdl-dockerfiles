FROM quay.io/pacbio/pb_wdl_base@sha256:1fa8b22abd2e32afe18cb2ff41dd0426d813a0450515efc0e626cf35f5658af3

LABEL maintainer="Tom Mokveld <tmokveld@pacificbiosciences.com>"

ARG IMAGE_NAME
ENV IMAGE_NAME "${IMAGE_NAME}"
ARG IMAGE_TAG
ENV IMAGE_TAG "${IMAGE_TAG}"

ARG ISOCALL_VERSION
RUN wget https://github.com/PacificBiosciences/isocall/releases/download/${ISOCALL_VERSION}/isocall-v${ISOCALL_VERSION}-linux-x86_64.gz \
	&& gunzip -c isocall-v${ISOCALL_VERSION}-linux-x86_64.gz > /usr/local/bin/isocall \
	&& chmod +x /usr/local/bin/isocall
