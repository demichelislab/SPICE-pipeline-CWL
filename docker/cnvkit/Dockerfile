FROM python:slim-buster AS build

ARG pkg_version

WORKDIR /build

ENV pkg_name cnvkit
ENV install_path /installed/local

RUN apt-get update                                                          && \
    apt-get install -y                                                         \
        wget                                                                   \
                                                                            && \
    wget https://github.com/etal/${pkg_name}/archive/v${pkg_version}.tar.gz && \
    tar -xzf v${pkg_version}.tar.gz                                         && \
    cd ${pkg_name}-${pkg_version}                                           && \
    python setup.py install --prefix ${install_path}


FROM python:slim-buster


RUN apt-get update                                                                                                        && \
    apt-get install -y                                                                                                       \
        r-base                                                                                                               \
                                                                                                                          && \
    ln -s /usr/local/bin/cnvkit.py /usr/local/bin/cnvkit                                                                  && \
    Rscript -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org/'); BiocManager::install('DNAcopy')" && \
    rm -rf /build /var/lib/apt/lists/*

COPY --from=build /installed/ /usr/

CMD [ "/usr/local/bin/cnvkit" ]
