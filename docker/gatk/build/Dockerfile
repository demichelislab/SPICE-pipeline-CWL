FROM openjdk:8-jdk-alpine AS build

ARG pkg_version
ARG glibc_version="2.26-r0"

WORKDIR /build

ENV sgerrand_pub_key="https://raw.githubusercontent.com/sgerrand/alpine-pkg-glibc/master/sgerrand.rsa.pub"

RUN apk --no-cache add ca-certificates                                                                                              \
                       openssl                                                                                                      \
                       grep                                                                                                         \
                       maven                                                                                                        \
                                                                                                                                 && \
    wget -O sgerrand.rsa.pub "${sgerrand_pub_key}"                                                                               && \
    wget -O glibc.apk https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${glibc_version}/glibc-${glibc_version}.apk && \
    wget https://github.com/broadgsa/gatk-protected/archive/${pkg_version}.tar.gz                                                && \
    tar -xzf ${pkg_version}.tar.gz                                                                                               && \
    cd gatk-protected-${pkg_version}/                                                                                            && \
    sed -i -re 's_/bin/bash_/bin/sh_' public/src/main/scripts/shell/delete_maven_links.sh                                        && \
    sed -i -re '/<module>external-example<\/module>/d' public/pom.xml                                                            && \
    sed -i -re '/<version>0.5.2<\/version>/ s/2/3/' public/gatk-root/pom.xml                                                     && \
    mvn verify -P\!queue                                                                                                         && \
    mv protected/gatk-package-distribution/target/gatk-package-distribution-${pkg_version}.jar ../gatk.jar

FROM openjdk:8-jre-alpine

ADD rootfs/ /
COPY --from=build /build/gatk.jar /usr/share/gatk/
COPY --from=build /build/glibc.apk /
COPY --from=build /build/sgerrand.rsa.pub /etc/apk/keys/

RUN apk add --no-cache libgomp && apk add /glibc.apk && rm -f /glibc.apk

CMD [ "/usr/bin/gatk" ]
