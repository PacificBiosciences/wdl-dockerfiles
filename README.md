# wdl-dockerfiles

This repo houses definitions for the Docker images used by PacBio workflows.


## Directory structure

Docker image definitions are found in [the docker directory](docker). Docker images are split into two categories: images built from conda environments ([`docker/conda/*`](docker/conda)), which all use [the same base Dockerfile](docker/conda/Conda_dockerfile) to build images, and other images which are defined using their own Dockerfiles.

Regardless of the type of docker image being built, each docker image must minimally define a [`build.env` file](#the-build.env-file).

Example directory structure:
```
docker/
├── bwa_mem/
│   ├── build.env
│   └── Dockerfile
├── conda/
│   ├── Conda_dockerfile
│   ├── smrtcells_stats/
│   │   └── build.env
│   └── whatshap/
│       └── build.env
├── pb_human_wgs_scripts/
│   ├── build.env
│   └── Dockerfile
└── samtools/
    ├── build.env
    ├── Dockerfile
    └── scripts/
```

### The `build.env` file

The build env file is a newline-delimited file containing key=value pairs for any environment variables that should be available as `--build-arg`s during Docker image build. This file must minimally define the `IMAGE_TAG` that will be used to tag the built image.

For dockers built from conda environments, the `IMAGE_TAG` is derived as the first 7 characters of the `CONDA_ENV_REPO_HASH` and can be omitted.

Example `build.env` file:
```
IMAGE_TAG=1.14
SAMTOOLS_VERSION=0.0.19
```

The `IMAGE_NAME` is also automatically made available as a build arg, and is set to the name of the directory containing the `build.env` file.


### Images built from conda environments ([`docker/conda`](docker/conda))

These images all use [the same base Dockerfile](docker/conda/Conda_dockerfile).

Each docker image built from this base Dockerfile is defined by creating a directory named `<image_name>` that contains a `build.env` file at the path `docker/conda/<image_name>/build.env`. The `image_name` will also be used as the name of the conda environment that is installed and sourced in the image.

The `build.env` file must define the following environment variables:

- `CONDA_ENV_REPO`: The URL for the repo containing conda environment definitions
- `CONDA_ENV_REPO_HASH`: The `CONDA_ENV_REPO` hash to checkout during image build

Images will be named as `<image_name>:<image_tag>`, where `image_name` is the name of the directory where `build.env` is found and `image_tag` is the first seven characters of the `CONDA_ENV_REPO_HASH`.


### Other images

Images that are not built from conda environments are defined by creating a `Dockerfile` and a `build.env` file at the path `docker/<image_name>/{Dockerfile,build.env}`.

Images will be named `<image_name>:<image_tag>`, where `image_name` is the name of the directory where `build.env` and the `Dockerfile` are found, and `image_tag` is the value of `IMAGE_TAG` in the corresponding `build.env` file.


## Building Docker images

Docker images can be built using the [build_docker_images](util/build_docker_images) utility script. Running this script with no `-d` argument will build all conda and non-conda images in the repo.

### Build a single image

```bash
./util/build_docker_images -d docker/conda/smrtcell_stats
```

#### Build all images

```bash
./util/build_docker_images
```

#### Build and push all images

```bash
./util/build_docker_images -p
```

#### Build and push all images in the `docker/conda` directory, using the `pacbio` container registry

```bash
./util/build_docker_images -p -d docker/conda -c pacbio
```
