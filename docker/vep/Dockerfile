FROM alpine:3.6 AS build

ARG pkg_version

WORKDIR /build

ENV local_perl /usr/local/share/perl5
ENV site_perl ${local_perl}/site_perl

RUN apk --no-cache add openssl                                                           \
                       curl                                                              \
                       git                                                               \
                       gcc                                                               \
                       make                                                              \
                       zlib-dev                                                          \
                       musl-dev                                                          \
                       perl                                                              \
                       perl-dev                                                          \
                       perl-dbi                                                          \
                       perl-module-build                                                 \
                                                                                      && \
    wget https://github.com/Ensembl/ensembl-vep/archive/release/${pkg_version}.tar.gz && \
    tar -xzf ${pkg_version}.tar.gz                                                    && \
    cd ensembl-vep-release-${pkg_version}                                             && \
    echo 'y' | ./INSTALL.pl -n -a a -d ${site_perl}                                   && \
    cp -a modules/Bio/EnsEMBL/VEP ${site_perl}/Bio/EnsEMBL/                           && \
    mkdir ../bin                                                                      && \
    for c in vep filter_vep variant_recoder; do cp ${c} ../bin/; done

FROM alpine:3.6

RUN apk --no-cache add perl      \
                       perl-dbi  \
                       perl-json

COPY --from=build ${site_perl} ${local_perl}
COPY --from=build /build/bin/* /usr/local/bin/

CMD [ "/usr/local/bin/vep" ]

