FROM python:3.9-buster

LABEL maintainer="Gregory Concepcion <gconcepcion@pacificbiosciences.com>"

ARG IMAGE_NAME
ARG IMAGE_TAG
ARG TOOL_VER

ENV IMAGE_NAME="${IMAGE_NAME}" \
    IMAGE_TAG="${IMAGE_TAG}" \
    TOOL_VER="${TOOL_VER}" \
    PATH="${PATH}:/opt/scripts"


RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        tar \
        wget \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    cd /opt && \
    wget https://github.com/PacificBiosciences/pb-StarPhase/archive/refs/tags/v${TOOL_VER}.tar.gz && \
    tar zxvf v${TOOL_VER}.tar.gz && \
    ln -s pbStarPhase-${TOOL_VER} ./pbstarphase && \
    ln pb-StarPhase-${TOOL_VER}/data/v${TOOL_VER}/cpic_20240404.json.gz ./pbstarphase_db.json.gz && \
    mkdir /opt/bin && \
    cd /opt/bin && \
    wget https://github.com/PacificBiosciences/pb-StarPhase/releases/download/v${TOOL_VER}/pbstarphase-v${TOOL_VER}-x86_64-unknown-linux-gnu.tar.gz && \
    tar zxvf pbstarphase-v${TOOL_VER}-x86_64-unknown-linux-gnu.tar.gz && \
    chmod +x pbstarphase-v${TOOL_VER}-x86_64-unknown-linux-gnu/pbstarphase && \
    ln -s /opt/bin/pbstarphase-v${TOOL_VER}-x86_64-unknown-linux-gnu/pbstarphase /opt/bin

ENV PATH /opt/bin:$PATH
