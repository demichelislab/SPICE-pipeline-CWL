FROM openjdk:8-jre-alpine

ARG picard_version

ADD rootfs/ /
ADD https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar /usr/share/picard/picard.jar

CMD [ "/usr/bin/gatk" ]
