FROM alpine:3.5 AS build

ARG pkg_version

WORKDIR /build

ENV python_path installed/lib/python2.7/site-packages

RUN apk --no-cache add ca-certificates                                                    \
                       openssl                                                            \
                       python2                                                            \
                       python2-dev                                                        \
                       py-setuptools                                                      \
                       gcc                                                                \
                       g++                                                                \
                       make                                                               \
                       musl-dev                                                           \
                       curl-dev                                                           \
                       zlib-dev                                                           \
                                                                                       && \
    ln -s /usr/include/locale.h /usr/include/xlocale.h                                 && \
    wget https://github.com/broadinstitute/oncotator/archive/v${pkg_version}.tar.gz    && \
    tar -xzf v${pkg_version}.tar.gz                                                    && \
    cd oncotator-${pkg_version}/                                                       && \
    sed -i -re 's/pysam == 0.9.0/pysam == 0.10.0/' setup.py                            && \
    mkdir -p ${python_path}                                                            && \
    PYTHONPATH=$(pwd)/${python_path} python setup.py install --prefix $(pwd)/installed && \
    ln -s $(pwd)/installed /installed

FROM alpine:3.5

RUN apk --no-cache add python2       \
                       py-setuptools \
                       curl

COPY --from=build /installed/ /usr/

CMD [ "/bin/oncotator" ]
