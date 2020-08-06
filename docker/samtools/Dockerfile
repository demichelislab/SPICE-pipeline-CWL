FROM debian:buster-slim AS build

ARG pkg_version

WORKDIR /build

RUN apt-get update                                                                                             && \
    apt-get install -y wget                                                                                       \
                       bzip2                                                                                      \
                       build-essential                                                                            \
                       libncurses-dev                                                                             \
                       zlib1g-dev                                                                                 \
                       libbz2-dev                                                                                 \
                       liblzma-dev                                                                                \
                                                                                                               && \
    mkdir /install                                                                                             && \
    wget https://github.com/samtools/samtools/releases/download/${pkg_version}/samtools-${pkg_version}.tar.bz2 && \
    tar -xjf samtools-${pkg_version}.tar.bz2                                                                   && \
    cd samtools-${pkg_version}/                                                                                && \
    ./configure --prefix=/installed                                                                            && \
    make                                                                                                       && \
    make install

FROM debian:buster-slim

COPY --from=build /installed /usr

RUN apt-get update                    && \
    apt-get install -y libbz2-1.0        \
                       liblzma5          \
                       libncursesw6      \
                       zlib1g         && \
    rm -rf /var/lib/apt/lists/*

CMD [ "/usr/bin/samtools" ]
