FROM openjdk:8-jdk-alpine AS build

ARG gatk_version

WORKDIR /build

RUN apk --no-cache add ca-certificates                                                                         \
                       wget                                                                                    \
                       grep                                                                                    \
                       maven                                                                                && \
    wget https://github.com/broadgsa/gatk-protected/archive/${gatk_version}.zip                             && \
    unzip ${gatk_version}.zip                                                                               && \
    cd gatk-protected-${gatk_version}/                                                                      && \
    sed -i -re 's_/bin/bash_/bin/sh_' public/src/main/scripts/shell/delete_maven_links.sh                   && \
    mvn verify -P\!queue                                                                                    && \
    mv protected/gatk-package-distribution/target/gatk-package-distribution-${gatk_version}.jar ../gatk.jar

FROM openjdk:8-jre-alpine

ADD rootfs/ /
COPY --from=build /build/gatk.jar /usr/share/gatk/

CMD [ "/usr/bin/gatk" ]
