FROM openjdk:8-jre-slim-buster

ARG pkg_version

ADD rootfs/ /

RUN apt-get update                                                                                                   && \
    apt-get install -y                                                                                                  \
      wget                                                                                                              \
      libgomp1                                                                                                          \
      zip                                                                                                               \
                                                                                                                     && \
    mkdir /build                                                                                                     && \
    cd /build                                                                                                        && \
    wget https://github.com/broadinstitute/gatk/releases/download/${pkg_version}/gatk-${pkg_version}.zip -O gatk.zip && \
    unzip gatk.zip                                                                                                   && \
    mkdir -p /usr/share/gatk/                                                                                        && \
    cp gatk-${pkg_version}/gatk-package-${pkg_version}-local.jar /usr/share/gatk/gatk.jar                            && \
    cd /                                                                                                             && \
    rm -rf /build /var/lib/apt/lists/*

CMD [ "/usr/bin/gatk" ]
