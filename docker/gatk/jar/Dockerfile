FROM openjdk:8-jre-alpine

ARG gatk_jar
ARG glibc_version="2.26-r0"

ENV sgerrand_pub_key="https://raw.githubusercontent.com/sgerrand/alpine-pkg-glibc/master/sgerrand.rsa.pub"

ADD rootfs/ /
ADD ${gatk_jar} /usr/share/gatk/gatk.jar

RUN apk --no-cache add ca-certificates                                                                                 \
                       openssl                                                                                         \
                                                                                                                    && \
    wget -q -O /etc/apk/keys/sgerrand.rsa.pub "${sgerrand_pub_key}"                                                 && \
    wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/${glibc_version}/glibc-${glibc_version}.apk && \
    apk add glibc-${glibc_version}.apk


CMD [ "/usr/bin/gatk" ]
