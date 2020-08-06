FROM debian:buster-slim AS build

ARG pkg_version

WORKDIR /build

ENV site_perl /usr/local/lib/site_perl

RUN apt-get update                                                                    && \
    apt-get install -y                                                                   \
        wget                                                                             \
        libfindbin-libs-perl                                                             \
        libdbi-perl                                                                      \
        libhttp-tiny-perl                                                                \
        libarchive-extract-perl                                                          \
        libarchive-zip-perl                                                              \
                                                                                      && \
    wget https://github.com/Ensembl/ensembl-vep/archive/release/${pkg_version}.tar.gz && \
    tar -xzf ${pkg_version}.tar.gz                                                    && \
    cd ensembl-vep-release-${pkg_version}                                             && \
    mkdir -p ${site_perl}                                                             && \
    echo 'y' | ./INSTALL.pl -n -a a -l -d ${site_perl} --NO_TEST                      && \
    cp -a modules/Bio/EnsEMBL/VEP ${site_perl}/Bio/EnsEMBL/                           && \
    mkdir ../bin                                                                      && \
    for c in vep filter_vep variant_recoder; do cp ${c} ../bin/; done

FROM debian:buster-slim

RUN apt-get update                     && \
    apt-get install -y                    \
        libdbi-perl                       \
        libtry-tiny-perl                  \
                                       && \
    rm -rf /build /var/lib/apt/lists/*

COPY --from=build ${site_perl} ${site_perl}
COPY --from=build /build/bin/* /usr/local/bin/

CMD [ "/usr/local/bin/vep" ]
