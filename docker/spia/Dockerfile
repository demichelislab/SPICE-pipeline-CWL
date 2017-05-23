FROM alpine:3.5

ENV pkg_folder /usr/share/spia

RUN apk --no-cache add R                                                                       \
                       openssl                                                                 \
                                                                                            && \
    mkdir /build                                                                            && \
    cd /build                                                                               && \
    wget https://demichelislab.unitn.it/lib/exe/fetch.php?media=spia.zip -O spia.zip        && \
    unzip spia.zip                                                                          && \
    cd SPIA                                                                                 && \
    sed -i -re 's|SPIAfunctions_location|"'${pkg_folder}'/src/SPIAfunctions.R"|' Bin/SPIA.R && \
    mkdir -p ${pkg_folder}/src                                                              && \
    mv Bin/SPIA.R Bin/SPIAfunctions.R ${pkg_folder}/src/                                    && \
    ln -s ${pkg_folder}/src/SPIA.R /usr/bin/                                                && \
    ln -s ${pkg_folder}/src/SPIA.R /usr/bin/spia                                            && \
    cd /                                                                                    && \
    rm -rf /build/

CMD [ "/usr/bin/SPIA.R" ]
