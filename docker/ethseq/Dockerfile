FROM debian:buster-slim

ARG pkg_version

ADD rootfs/ /

RUN apt-get update                     && \
    apt-get install -y                    \
      r-base                              \
                                       && \
    /build/bin/install_r_packages.R    && \
    rm -rf /build /var/lib/apt/lists/*

CMD [ "/usr/bin/Rscript" ]
