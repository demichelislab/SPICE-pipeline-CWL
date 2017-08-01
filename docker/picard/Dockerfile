FROM openjdk:8-jre-alpine

ARG pkg_version

ADD rootfs/ /
ADD https://github.com/broadinstitute/picard/releases/download/${pkg_version}/picard.jar /usr/share/picard/picard.jar
RUN chmod 644 /usr/share/picard/picard.jar

CMD [ "/usr/bin/picard" ]
