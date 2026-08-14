FROM quay.io/pacbio/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae

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

ARG TRUVARI_VERSION
ARG PYWFA_VERSION
ARG PYABPOA_VERSION
RUN uv pip install --no-build-isolation \
	"setuptools<82.0.0" Cython

# truvari>=5.4.0 depends directly on pywfa (not just the [bwa] extra), and
# pywfa only ships as an sdist -- it always compiles its vendored WFA2-lib
# from source. WFA2-lib's own Makefile hardcodes `-march=native`, so
# building on this host bakes in AVX-512 instructions and SIGILLs on any
# deployment node without AVX-512 (e.g. Zen 2 EPYC 7H12 nodes on compute9).
# Pre-build pywfa here with a portable x86-64-v3 (AVX2, no AVX-512)
# baseline so the truvari install below finds it already satisfied and
# doesn't rebuild it with -march=native.
RUN pywfa_url=$(curl -sL "https://pypi.org/pypi/pywfa/${PYWFA_VERSION}/json" \
		| python3 -c "import sys, json; print(json.load(sys.stdin)['urls'][0]['url'])") \
	&& mkdir -p /tmp/pywfa_src && cd /tmp/pywfa_src \
	&& curl -sL -o "pywfa-${PYWFA_VERSION}.tar.gz" "${pywfa_url}" \
	&& tar --extract --file "pywfa-${PYWFA_VERSION}.tar.gz" \
	&& sed -i 's/-march=native/-march=x86-64-v3/' \
		"pywfa-${PYWFA_VERSION}/pywfa/WFA2_lib/Makefile" \
	&& uv pip install --no-build-isolation "./pywfa-${PYWFA_VERSION}" \
	&& cd / && rm -rf /tmp/pywfa_src

# truvari also depends directly on pyabpoa (used for the poa refine align
# backend) -- same -march=native default SIMD problem, but its setup.py
# already supports picking a portable target via an env var, no Makefile
# patch needed. AVX2=1 selects -mavx2 instead of -march=native.
RUN pyabpoa_url=$(curl -sL "https://pypi.org/pypi/pyabpoa/${PYABPOA_VERSION}/json" \
		| python3 -c "import sys, json; print(json.load(sys.stdin)['urls'][0]['url'])") \
	&& mkdir -p /tmp/pyabpoa_src && cd /tmp/pyabpoa_src \
	&& curl -sL -o "pyabpoa-${PYABPOA_VERSION}.tar.gz" "${pyabpoa_url}" \
	&& tar --extract --file "pyabpoa-${PYABPOA_VERSION}.tar.gz" \
	&& AVX2=1 uv pip install --no-build-isolation "./pyabpoa-${PYABPOA_VERSION}" \
	&& cd / && rm -rf /tmp/pyabpoa_src

RUN uv pip install --no-build-isolation \
	"truvari[bwa]==${TRUVARI_VERSION}"

# truvari's phab.py (used by refine) relies on pysam's samtools.faidx()
# returning captured command stdout as a string. pysam 0.24.0 -- the only
# pysam release with Python 3.14 support, since it's the first release
# published after Python 3.14 existed -- silently stopped capturing
# faidx's output specifically (confirmed by diffing pysam's own dispatch
# code between 0.23.3 and 0.24.0: faidx isn't in the MAP_STDOUT_OPTIONS
# dict that controls which methods get their stdout redirected to a named
# file, so the real output leaks to the process's actual stdout and
# faidx() returns an empty string instead). No pysam release both
# supports Python 3.14 and predates this regression, so patch truvari
# itself: write directly via `samtools faidx`'s own -o flag instead of
# relying on captured stdout.
COPY patch_phab.py /tmp/patch_phab.py
RUN python3 /tmp/patch_phab.py && rm /tmp/patch_phab.py
