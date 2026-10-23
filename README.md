# wdl-dockerfiles

This repo houses definitions for the Docker images used by PacBio workflows.

## Directory structure

Docker image definitions are found in [the docker directory](docker).

Each Docker image must minimally define a [`build.env` file](#the-buildenv-file).

Example directory structure:

```text
docker/
├── bwa_mem/
│   ├── build.env
│   └── Dockerfile
├── Conda_dockerfile # a template dockerfile, used to build multiple images
├── htslib/
│   └── build.env # Uses the Conda_dockerfile
├── pb_human_wgs_scripts/
│   ├── build.env
│   └── Dockerfile
├── samtools/
│   ├── build.env
│   ├── Dockerfile
│   └── scripts/
└── yak
    └── build.env # Uses the Conda_dockerfile
```

## The `build.env` file

Each target image is defined via the presence of a `build.env` file, which is used to specify the name and version tag for the corresponding Docker image. It must contain at minimum the following variables:

- `IMAGE_NAME`
- `IMAGE_TAG`

All variables defined in the `build.env` file will be made available as build arguments during Docker image build.

The `IMAGE_TAG` variable can be built using other variables defined in the `build.env` file, as long as those other variables are defined before `IMAGE_TAG`. For example, the following `IMAGE_TAG` would be set to `0.7.8_1.15`:

```shell
# Tool verisons
BWA_VERSION=0.7.8
SAMTOOLS_VERSION=1.15

# Image info
IMAGE_NAME=bwa_samtools
IMAGE_TAG=${BWA_VERSION}_${SAMTOOLS_VERSION}
```

The `Dockerfile` corresponding to a `build.env` must either:

1. Be named `Dockerfile` and exist in the same directory as the `build.env` file, or
2. Be specified using the `DOCKERFILE` variable in the `build.env` file (e.g. for image builds that reuse a common Dockerfile). The path to the Dockerfile should be relative to the directory containing the `build.env` file.

### Special variables

These variables have special meanings when the Docker image is being built. † specifies that the variable is required.

- `IMAGE_NAME`†: specifies the name of the built image
- `IMAGE_TAG`†: specifies the tag of the built image
- `NOBUILD`: When set to `true`, skip building this image
- `DOCKERFILE`: Specify an alternate path to a Dockerfile, relative to the `build.env` file. If left undefined, the Dockerfile must be named `Dockerfile` and be present at the same directory level as the corresponding `build.env` file.
- `CONDA_ENVIRONMENT_TEMPLATE`: For conda-based images, set this to the path (relative to the `build.env` file) to the conda environment template file. This file may have variables set in place of tool versions and will be populated with the variables from `build.env` at build time.

## Building Docker images

Docker images can be built using the [build_docker_images](util/build_docker_images) utility script.

### Container engine and platform

The script works with [podman](https://podman.io/), [docker](https://www.docker.com/), or Apple's [container](https://github.com/apple/container). By default it uses the first of `podman`, `docker`, `container` found on `PATH`. To choose one, pass `-e <engine>` or set `CONTAINER_ENGINE`.

Images are built and pulled for `linux/amd64` by default, including on Apple Silicon. To override, pass `-a <platform>` or set `CONTAINER_PLATFORM`.

Both environment variables also work through `make`, e.g. `CONTAINER_ENGINE=container make build-hifiasm`.

### Image digests

After each push, the script looks up the image's `sha256` digest from the remote registry. It prints a summary of `registry/image:tag@sha256:...` references for use in commit and PR messages. The lookup uses [crane](https://github.com/google/go-containerregistry/tree/main/cmd/crane). Without crane, podman and docker fall back to the digest recorded locally by the push; Apple `container` prints a warning and no digest.

To install crane:

- macOS: `brew install crane`
- Linux: download the release tarball for your platform from [go-containerregistry releases](https://github.com/google/go-containerregistry/releases) and put `crane` on your `PATH`. If Go is installed, `go install github.com/google/go-containerregistry/cmd/crane@latest` also works.

crane reads registry credentials from `~/.docker/config.json`. If your engine stores logins elsewhere, run `crane auth login <registry>` once.

### Apple `container` prerequisites

- The script uses `readlink -f`, which needs macOS 12.3 or later.
- Start the service with `container system start`.
- Building `linux/amd64` images on Apple Silicon needs Rosetta: `softwareupdate --install-rosetta`.
- `container registry login` stores credentials in the macOS Keychain, which crane cannot read. Run `crane auth login <registry>` once so the digest lookup works.
- Conda-based images need `envsubst`, which is not on stock macOS: `brew install gettext`.

### Build a single image

```bash
./util/build_docker_images -d docker/smrtcell_stats
```

#### Build all images in the `docker` directory

```bash
./util/build_docker_images -d docker
```

#### Build and push all images in the `docker` directory

```bash
./util/build_docker_images -d docker -p
```

#### Build and push all images in the `docker` directory, using the `pacbio` container registry

```bash
./util/build_docker_images -d docker -p -c pacbio
```
