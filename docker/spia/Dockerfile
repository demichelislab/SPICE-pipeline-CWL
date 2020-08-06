FROM r-base

ENV pkg_folder /usr/share/spia

RUN mkdir /build                                                                            && \
    cd /build                                                                               && \
    wget http://demichelislab.unitn.it/lib/exe/fetch.php?media=spia.zip -O spia.zip         && \
    unzip spia.zip                                                                          && \
    cd SPIA                                                                                 && \
    sed -i -re 's|SPIAfunctions_location|"'${pkg_folder}'/src/SPIAfunctions.R"|' Bin/SPIA.R && \
    mkdir -p ${pkg_folder}/src                                                              && \
    mv Bin/SPIA.R Bin/SPIAfunctions.R ${pkg_folder}/src/                                    && \
    ln -s ${pkg_folder}/src/SPIA.R /usr/bin/                                                && \
    ln -s ${pkg_folder}/src/SPIA.R /usr/bin/spia                                            && \
    chmod o+x ${pkg_folder}/src/SPIA.R                                                      && \
    cd /                                                                                    && \
    rm -rf /build/

CMD [ "/usr/bin/SPIA.R" ]
