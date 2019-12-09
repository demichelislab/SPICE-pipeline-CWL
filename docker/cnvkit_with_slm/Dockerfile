FROM alpine:3.6 AS build

ADD rootfs/ /

WORKDIR /build

ENV slmsuite_version 1.0
ENV pkg_name cnvkit
ENV install_path /installed
ENV python_path ${install_path}/lib/python2.7/site-packages

RUN apk --no-cache add ca-certificates                                                              \
                       git                                                                          \
                       openssl                                                                      \
                       python2                                                                      \
                       python2-dev                                                                  \
                       openblas-dev                                                                 \
                       py-setuptools                                                                \
                       py2-numpy                                                                    \
                       gcc                                                                          \
                       g++                                                                          \
                       musl-dev                                                                     \
                       make                                                                         \
                       cmake                                                                        \
                       py-numpy-dev                                                                 \
                       zlib-dev                                                                     \
                       bzip2-dev                                                                    \
                       xz-dev                                                                       \
                       freetype-dev                                                                 \
                       jpeg-dev                                                                     \
                       py2-cairo-dev                                                                \
                       py2-numpy-f2py                                                               \
                       R                                                                            \
                       R-dev                                                                        \
                                                                                                 && \
    wget http://www.qhull.org/download/qhull-2015-src-7.2.0.tgz                                  && \
    mkdir ${install_path}                                                                        && \
    tar -xzf qhull-2015-src-7.2.0.tgz                                                            && \
    cd qhull-2015.2                                                                              && \
    cmake -DCMAKE_INSTALL_PREFIX=${install_path} .                                               && \
    make                                                                                         && \
    make install                                                                                 && \
    cd /build                                                                                    && \
    git clone https://github.com/hcraT/${pkg_name} -b add_slm --depth=1                          && \
    cd ${pkg_name}                                                                               && \
    mkdir -p ${python_path}                                                                      && \
    PYTHONPATH=${python_path} python setup.py install --prefix ${install_path}                   && \
    ln -s /usr/bin/cnvkit.py /installed/bin/cnvkit                                               && \
    cd /build                                                                                    && \
    mkdir ./tmp                                                                                  && \
    cd ./tmp                                                                                     && \
    wget https://downloads.sourceforge.net/project/slmsuite/slmsuite-${slmsuite_version}.tar.bz2 && \
    tar -xjf slmsuite-${slmsuite_version}.tar.bz2                                                && \
    cd ..                                                                                        && \
    ./bin/install_r_packages.R


FROM alpine:3.6

RUN apk --no-cache add python2       \
                       py-setuptools \
                       py2-numpy     \
                       R

COPY --from=build /installed/ /usr/

CMD [ "/usr/bin/cnvkit" ]
