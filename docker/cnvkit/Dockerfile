FROM debian:buster-slim

ARG pkg_version

RUN apt-get update                                                                                                        && \
    apt-get install -y                                                                                                       \
        liblzma-dev                                                                                                          \
        python3-biopython                                                                                                    \
        python3-dev                                                                                                          \
        python3-matplotlib                                                                                                   \
        python3-numpy                                                                                                        \
        python3-pip                                                                                                          \
        python3-reportlab                                                                                                    \
        python3-scipy                                                                                                        \
        python3-tk                                                                                                           \
        r-base-core                                                                                                          \
        zlib1g-dev                                                                                                           \
                                                                                                                          && \
    pip3 install -U Cython                                                                                                && \
    pip3 install -U future futures pandas pomegranate pyfaidx pysam                                                       && \
    pip3 install cnvkit==${pkg_version}                                                                                   && \
    Rscript -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org/'); BiocManager::install('DNAcopy')" && \
    ln -s /usr/local/bin/cnvkit.py /usr/local/bin/cnvkit                                                                  && \
    rm -rf /var/lib/apt/lists/*

CMD [ "/usr/local/bin/cnvkit" ]
