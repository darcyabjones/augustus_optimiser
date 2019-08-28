ARG IMAGE="darcyabjones/base:panann-v0.0.1"

FROM "${IMAGE}" as htslib_builder

ARG HTSLIB_TAG="1.9"
ARG SAMTOOLS_TAG="1.9"

ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"

ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"

ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       autoconf \
       build-essential \
       ca-certificates \
       git \
       libbz2-dev \
       libcurl4-gnutls-dev \
       libgsl-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       libperl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates

# CA cert stuff required for git clone https

WORKDIR /tmp/htslib
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && git clone "${HTSLIB_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${HTSLIB_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure --prefix="${HTSLIB_PREFIX}" --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install \
  && add_runtime_dep libbz2-1.0 libcurl3-gnutls libssl1.1 lzma zlib1g
# htslib: libssl (providing libcrypto) only required if need Amazon S3 support.

WORKDIR /tmp/samtools
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && git clone "${SAMTOOLS_REPO}" . \
  && git fetch --tags \
  && git checkout "tags/${SAMTOOLS_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure \
       --prefix="${SAMTOOLS_PREFIX}" \
       --with-htslib="${HTSLIB_PREFIX}" \
  && make -j $(grep -c ^processor /proc/cpuinfo) \
  && make -j $(grep -c ^processor /proc/cpuinfo) install \
  && add_runtime_dep libncurses5 zlib1g
# Samtools also depends on htslib and it's dependencies.


FROM "${IMAGE}" as genometools_builder
ARG GENOMETOOLS_VERSION="1.5.10"
ARG GENOMETOOLS_URL="http://genometools.org/pub/genometools-${GENOMETOOLS_VERSION}.tar.gz"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       libcairo2-dev \
       libpango1.0-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O genometools.tar.gz "${GENOMETOOLS_URL}" \
  && tar zxf genometools.tar.gz \
  && rm genometools.tar.gz \
  && cd genometools*/ \
  && sed -i 's/-Wall//g' Makefile \
  && make errorcheck=no \
  && make errorcheck=no prefix="${GENOMETOOLS_PREFIX}" install \
  && add_runtime_dep libcairo2 libpango-1.0-0 libpangocairo-1.0-0


FROM "${IMAGE}" as aegean_builder

ARG AEGEAN_VERSION="v0.15.0"
ARG AEGEAN_URL="https://github.com/standage/AEGeAn/archive/v0.15.0.tar.gz"
ARG AEGEAN_PREFIX_ARG="/opt/aegean/${AEGEAN_VERSION}"
ENV AEGEAN_PREFIX="${AEGEAN_PREFIX_ARG}"

ARG GENOMETOOLS_VERSION="1.5.10"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=genometools_builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=genometools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       build-essential \
       libcairo2-dev \
       libpango1.0-dev \
       python \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && wget -O aegean.tar.gz "${AEGEAN_URL}" \
  && tar zxf aegean.tar.gz \
  && rm aegean.tar.gz

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && cd AEGeAn*/ \
  && sed -i "s~/usr/local/include/genometools~${GENOMETOOLS_PREFIX}/include/genometools~" Makefile \
  && make test \
  && make prefix="${AEGEAN_PREFIX}" install install-scripts \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"


FROM "${IMAGE}" as augustus_builder

ARG AUGUSTUS_COMMIT="8b1b14a7489e4545e89c8725dc33268f6c2a9117"
ARG AUGUSTUS_REPO="https://github.com/Gaius-Augustus/Augustus.git"
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"

ARG HTSLIB_TAG="1.9"
ARG HTSLIB_REPO="https://github.com/samtools/htslib.git"

ARG SAMTOOLS_TAG="1.9"
ARG SAMTOOLS_REPO="https://github.com/samtools/samtools.git"

# Install htslib and samtools.
# Bam2wig in augustus depends on some intermediate samtools/htslib compilation
# rather than the actual headers/shared libraries, so I have to compile it
# separately.

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && apt-get update \
  && apt-get install -y \
       autoconf \
       build-essential \
       ca-certificates \
       git \
       libbz2-dev \
       libcurl4-gnutls-dev \
       liblzma-dev \
       libncurses5-dev \
       libssl-dev \
       zlib1g-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates

WORKDIR /tmp/htslib
RUN  git clone ${HTSLIB_REPO} . \
  && git fetch --tags \
  && git checkout "tags/${HTSLIB_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure --enable-libcurl \
  && make -j $(grep -c ^processor /proc/cpuinfo)

WORKDIR /tmp/samtools
RUN  git clone ${SAMTOOLS_REPO} . \
  && git fetch --tags \
  && git checkout "tags/${SAMTOOLS_TAG}" \
  && autoheader \
  && autoconf \
  && ./configure \
  && make -j $(grep -c ^processor /proc/cpuinfo)

# Install augustus

# This is for bam2wig
ENV TOOLDIR="/tmp"

WORKDIR /tmp/augustus
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       libbamtools-dev \
       libboost-all-dev \
       libboost-iostreams-dev \
       libboost-graph-dev \
       libcurl4-gnutls-dev \
       libgsl-dev \
       liblpsolve55-dev \
       libssl-dev \
       libsuitesparse-dev \
       libsqlite3-dev \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${AUGUSTUS_REPO}" . \
  && git fetch --tags \
  && git checkout "${AUGUSTUS_COMMIT}" \
  && mkdir bin \
  && sed -i "s/# SQLITE = true/SQLITE = true/g" common.mk \
  && sed -i "s/# COMPGENEPRED = true/COMPGENEPRED = true/g" common.mk \
  && sed -i 's~INSTALLDIR = .*~INSTALLDIR="${AUGUSTUS_PREFIX}"~g' Makefile \
  && cd auxprogs/bam2wig \
  && make \
  && cd /tmp/augustus \
  && make \
  && make install \
  && make test \
  && add_runtime_dep \
       libbamtools2.5.1 \
       libcurl3-gnutls \
       libgsl23 \
       libparallel-forkmanager-perl \
       libssl1.1 \
       libsqlite3-0 \
       lp-solve \
       zlib1g


FROM "${IMAGE}"

ARG GENOMETOOLS_VERSION="1.5.10"
ARG GENOMETOOLS_PREFIX_ARG="/opt/genometools/${GENOMETOOLS_VERSION}"
ENV GENOMETOOLS_PREFIX="${GENOMETOOLS_PREFIX_ARG}"
LABEL genometools.version="${GENOMETOOLS_VERSION}"

ENV PATH="${GENOMETOOLS_PREFIX}/bin:${PATH}"
ENV INCLUDE="${GENOMETOOLS_PREFIX}/include:${INCLUDE}"
ENV CPATH="${GENOMETOOLS_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${GENOMETOOLS_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=genometools_builder "${GENOMETOOLS_PREFIX}" "${GENOMETOOLS_PREFIX}"
COPY --from=genometools_builder "${APT_REQUIREMENTS_FILE}" /build/apt/genometools.txt


ARG AEGEAN_VERSION="v0.15.0"
ARG AEGEAN_PREFIX_ARG="/opt/aegean/${AEGEAN_VERSION}"
ENV AEGEAN_PREFIX="${AEGEAN_PREFIX_ARG}"
LABEL aegean.version="${AEGEAN_VERSION}"

ENV PATH="${AEGEAN_PREFIX}/bin:${PATH}"
ENV INCLUDE="${AEGEAN_PREFIX}/include:${INCLUDE}"
ENV CPATH="${AEGEAN_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${AEGEAN_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${AEGEAN_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV LD_RUN_PATH="${AEGEAN_PREFIX}/lib:${LD_RUN_PATH}"

COPY --from=aegean_builder "${AEGEAN_PREFIX}" "${AEGEAN_PREFIX}"
COPY --from=aegean_builder "${APT_REQUIREMENTS_FILE}" /build/apt/aegean.txt


ARG AUGUSTUS_COMMIT="8b1b14a7489e4545e89c8725dc33268f6c2a9117"
ARG AUGUSTUS_PREFIX_ARG="/opt/augustus/${AUGUSTUS_COMMIT}"
ENV AUGUSTUS_PREFIX="${AUGUSTUS_PREFIX_ARG}"
LABEL augustus.version="${AUGUSTUS_COMMIT}"

ENV PATH="${AUGUSTUS_PREFIX}/bin:${AUGUSTUS_PREFIX}/scripts:${PATH}"
ENV AUGUSTUS_CONFIG_PATH="${AUGUSTUS_PREFIX}/config"

COPY --from=augustus_builder "${AUGUSTUS_PREFIX}" "${AUGUSTUS_PREFIX}"
COPY --from=augustus_builder "${APT_REQUIREMENTS_FILE}" /build/apt/augustus.txt

ARG HTSLIB_TAG="1.9"
ARG SAMTOOLS_TAG="1.9"
ARG HTSLIB_PREFIX_ARG="/opt/htslib/${HTSLIB_TAG}"
ARG SAMTOOLS_PREFIX_ARG="/opt/samtools/${SAMTOOLS_TAG}"
ENV HTSLIB_PREFIX="${HTSLIB_PREFIX_ARG}"
ENV SAMTOOLS_PREFIX="${SAMTOOLS_PREFIX_ARG}"
LABEL htslib.version="${HTSLIB_TAG}"
LABEL samtools.version="${SAMTOOLS_TAG}"

ENV PATH="${SAMTOOLS_PREFIX}/bin:${BCFTOOLS_PREFIX}/bin:${HTSLIB_PREFIX}/bin:${PATH}"
ENV CPATH="${HTSLIB_PREFIX}/include:${CPATH}"
ENV LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LIBRARY_PATH}"
ENV LD_LIBRARY_PATH="${HTSLIB_PREFIX}/lib:${LD_LIBRARY_PATH}"

COPY --from=htslib_builder "${HTSLIB_PREFIX}" "${HTSLIB_PREFIX}"
COPY --from=htslib_builder "${SAMTOOLS_PREFIX}" "${SAMTOOLS_PREFIX}"
COPY --from=htslib_builder "${APT_REQUIREMENTS_FILE}" /build/apt/htslib.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && apt-get install -y --no-install-recommends libmpich-dev mpich python3 python3-dev python3-pip python3-setuptools python3-wheel \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"


COPY . /tmp
WORKDIR /tmp
RUN pip3 install .



# This is useful for testing.
# COPY --from=builder "/tmp/augustus/examples" "${AUGUSTUS_PREFIX}/examples"
# RUN augustus --species=human --UTR=on ${AUGUSTUS_PREFIX}/examples/example.fa
