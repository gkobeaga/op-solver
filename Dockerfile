FROM ubuntu:14.04

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

ARG PACKAGES="git"
ARG BUILD_DEPS="build-essential autoconf automake libtool m4 libgmp-dev"

RUN apt-get update
RUN apt-get install -y $PACKAGES
RUN apt-get install -y $BUILD_DEPS

# Install op-solver
WORKDIR /src
RUN set -x \
  && git clone --depth 1  https://github.com/gkobeaga/op-solver \
  && ( \
    cd op-solver \
    && ./autogen.sh \
    && mkdir build && cd build \
    && ../configure \
    && make \
    && make check \
    && make distcheck \
    && make install \
  )

# Download OPLib
WORKDIR /
RUN set -x \
  && git clone --depth 1  https://github.com/bcamath-ds/OPLib

WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/op-solver"]
