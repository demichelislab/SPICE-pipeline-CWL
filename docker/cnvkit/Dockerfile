FROM alpine:3.6 AS build

ARG pkg_version

WORKDIR /build

ENV pkg_name cnvkit
ENV install_path /installed
ENV python_path ${install_path}/lib/python2.7/site-packages

RUN apk --no-cache add ca-certificates                                             \
                       openssl                                                     \
                       python2                                                     \
                       python2-dev                                                 \
                       openblas-dev                                                \
                       py-setuptools                                               \
                       py2-numpy                                                   \
                       gcc                                                         \
                       g++                                                         \
                       musl-dev                                                    \
                       make                                                        \
                       cmake                                                       \
                       py-numpy-dev                                                \
                       zlib-dev                                                    \
                       bzip2-dev                                                   \
                       xz-dev                                                      \
                       freetype-dev                                                \
                       jpeg-dev                                                    \
                       py2-cairo-dev                                               \
                       py2-numpy-f2py                                              \
                                                                                && \
    wget http://www.qhull.org/download/qhull-2015-src-7.2.0.tgz                 && \
    mkdir ${install_path}                                                       && \
    tar -xzf qhull-2015-src-7.2.0.tgz                                           && \
    cd qhull-2015.2                                                             && \
    cmake -DCMAKE_INSTALL_PREFIX=${install_path} .                              && \
    make                                                                        && \
    make install                                                                && \
    wget https://github.com/etal/${pkg_name}/archive/v${pkg_version}.tar.gz     && \
    tar -xzf v${pkg_version}.tar.gz                                             && \
    cd ${pkg_name}-${pkg_version}/                                              && \
    mkdir -p ${python_path}                                                     && \
    PYTHONPATH=${python_path} python setup.py install --prefix ${install_path}

FROM alpine:3.6

RUN apk --no-cache add python2                                                                   \
                       py-setuptools                                                             \
                       py2-numpy                                                                 \
                       curl                                                                      \
                       R                                                                         \
                                                                                              && \
    ln -s /usr/bin/cnvkit.py /usr/bin/cnvkit                                                  && \
    Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('PSCBS', 'cghFLasso')"

COPY --from=build /installed/ /usr/

CMD [ "/usr/bin/cnvkit" ]
