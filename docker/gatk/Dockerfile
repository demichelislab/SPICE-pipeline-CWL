FROM openjdk:8-jre-alpine

ARG gatk_jar

ADD rootfs/ /
ADD ${gatk_jar} /usr/share/gatk/gatk.jar

CMD [ "/usr/bin/gatk" ]
