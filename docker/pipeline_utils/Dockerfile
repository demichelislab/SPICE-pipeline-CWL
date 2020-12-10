FROM r-base

ADD rootfs/ /

RUN /build/bin/install_r_packages.R && \
    rm -rf /build

CMD [ "/usr/bin/pipeline_utils" ]
