SHELL := /bin/bash

BUILD_SCRIPT := ./util/build_docker_images

# Target container registry; required for push and mirror targets.
# Usage: make push-hifiasm REMOTE_REPO=quay.io/pacbio
REMOTE_REPO ?=

# ──────────────────────────────────────────────
# Docker images (subdirectories of docker/)
# ──────────────────────────────────────────────
DOCKER_IMAGES := $(notdir $(patsubst %/,%,$(wildcard docker/*/)))

# ──────────────────────────────────────────────
# Mirror images (subdirectories of mirror/)
# Each entry is: <name>:<source-registry>
# ──────────────────────────────────────────────
MIRROR_DEEPVARIANT_SRC     := google
MIRROR_GLNEXUS_SRC         := ghcr.io/dnanexus-rnd/

MIRROR_IMAGES := deepvariant glnexus

# ──────────────────────────────────────────────
# Phony declarations
# ──────────────────────────────────────────────
.PHONY: help #\
	build-all $(addprefix build-,$(DOCKER_IMAGES)) \
	push-all   $(addprefix push-,$(DOCKER_IMAGES)) \
	mirror-all $(addprefix mirror-,$(MIRROR_IMAGES))

.DEFAULT_GOAL := help

# ──────────────────────────────────────────────
# Guard macro — enforces REMOTE_REPO is set
# ──────────────────────────────────────────────
define _require_remote_repo
	@if [[ -z "$(REMOTE_REPO)" ]]; then \
		echo ""; \
		echo "  ERROR: REMOTE_REPO is not set."; \
		echo "  Usage: make $@ REMOTE_REPO=<registry>  (e.g. quay.io/pacbio)"; \
		echo ""; \
		exit 1; \
	fi
endef

# ──────────────────────────────────────────────
# Help
# ──────────────────────────────────────────────
help: ## Show this help message
	@awk 'BEGIN { \
		FS = ":.*##"; \
		printf "\nUsage:\n  make \033[36m<target>\033[0m [REMOTE_REPO=<registry>]\n" \
	} \
	/^[a-zA-Z_0-9%-]+:.*?##/ { printf "  \033[36m%-32s\033[0m %s\n", $$1, $$2 } \
	/^##@/ { printf "\n\033[1m%s\033[0m\n", substr($$0, 5) } ' \
	$(MAKEFILE_LIST)
	@echo ""

##@ Build (local only)

# build-all: $(addprefix build-,$(DOCKER_IMAGES)) ## Build all docker images locally

build-%: ## Build a single image locally
	@if [[ ! -d docker/$* ]]; then \
		echo ""; \
		echo "  ERROR: Unknown image '$*'."; \
		echo "  Valid images: $(DOCKER_IMAGES)"; \
		echo ""; \
		exit 1; \
	fi
	$(BUILD_SCRIPT) -d docker/$*

##@ Push (requires REMOTE_REPO=<registry>)

# push-all: $(addprefix push-,$(DOCKER_IMAGES)) ## Build and push all images to REMOTE_REPO

push-%: ## Build and push a single image to REMOTE_REPO
	$(call _require_remote_repo)
	@if [[ ! -d docker/$* ]]; then \
		echo ""; \
		echo "  ERROR: Unknown image '$*'."; \
		echo "  Valid images: $(DOCKER_IMAGES)"; \
		echo ""; \
		exit 1; \
	fi
	$(BUILD_SCRIPT) -d docker/$* -p -c $(REMOTE_REPO)

##@ Mirror upstream images (requires REMOTE_REPO=<registry>)

# mirror-all: $(addprefix mirror-,$(MIRROR_IMAGES)) ## Mirror all upstream images to REMOTE_REPO

mirror-deepvariant: ## Mirror google/deepvariant to REMOTE_REPO
	$(call _require_remote_repo)
	$(BUILD_SCRIPT) -s $(MIRROR_DEEPVARIANT_SRC) -d mirror/deepvariant -p -c $(REMOTE_REPO)
	$(BUILD_SCRIPT) -s $(MIRROR_DEEPVARIANT_SRC) -d mirror/deepvariant_gpu -p -c $(REMOTE_REPO)

mirror-glnexus: ## Mirror ghcr.io/dnanexus-rnd/glnexus to REMOTE_REPO
	$(call _require_remote_repo)
	$(BUILD_SCRIPT) -s $(MIRROR_GLNEXUS_SRC) -d mirror/glnexus -p -c $(REMOTE_REPO)
