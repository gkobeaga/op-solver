FROM alpine:latest

ENV LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/

ARG PACKAGES="git"
ARG BUILD_DEPS="build-base autoconf automake libtool m4 gmp-dev"

RUN apk update
RUN apk --no-cache --update add $PACKAGES
RUN apk --no-cache --update add $BUILD_DEPS

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
    && make install \
  )

# Download OPLib
WORKDIR /
RUN set -x \
  && git clone --depth 1  https://github.com/bcamath-ds/OPLib

WORKDIR /tmp
ENTRYPOINT ["/usr/local/bin/op-solver"]
