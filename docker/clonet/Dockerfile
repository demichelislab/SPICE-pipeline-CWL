FROM alpine:3.7

ENV pkg_folder /usr/share/clonet

ADD rootfs/ /

RUN apk --no-cache add R                                 \
                       R-dev                             \
                       git                               \
                       g++                               \
                       openssl                           \
                                                      && \
    cd /build                                         && \
    ./bin/install_r_packages.R                        && \
    git clone https://bitbucket.org/deid00/clonet.git && \
    cd clonet                                         && \
    mkdir -p ${pkg_folder}/src                        && \
    mv Functions Tools CLONET.R ${pkg_folder}/src/    && \
    ln -s ${pkg_folder}/src/CLONET.R /usr/bin/        && \
    ln -s ${pkg_folder}/src/CLONET.R /usr/bin/clonet  && \
    cd /                                              && \
    rm -rf /build/

CMD [ "/usr/bin/clonet" ]
