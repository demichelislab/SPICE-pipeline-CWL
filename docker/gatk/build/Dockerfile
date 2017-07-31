FROM openjdk:8-jdk-alpine AS build

ARG pkg_version

WORKDIR /build

RUN apk --no-cache add ca-certificates                                                                        \
                       openssl                                                                                \
                       grep                                                                                   \
                       maven                                                                               && \
    wget https://github.com/broadgsa/gatk-protected/archive/${pkg_version}.tar.gz                          && \
    tar -xzf ${pkg_version}.tar.gz                                                                         && \
    cd gatk-protected-${pkg_version}/                                                                      && \
    sed -i -re 's_/bin/bash_/bin/sh_' public/src/main/scripts/shell/delete_maven_links.sh                  && \
    sed -i -re '/<module>external-example<\/module>/d' public/pom.xml                                      && \
    mvn verify -P\!queue                                                                                   && \
    mv protected/gatk-package-distribution/target/gatk-package-distribution-${pkg_version}.jar ../gatk.jar

FROM openjdk:8-jre-alpine

ADD rootfs/ /
COPY --from=build /build/gatk.jar /usr/share/gatk/

CMD [ "/usr/bin/gatk" ]
