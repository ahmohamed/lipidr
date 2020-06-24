FROM bioconductor/bioconductor_docker:devel
RUN R -e "BiocManager::install('lipidr', ask=F)"

RUN mkdir /root/lipidr
COPY . /root/lipidr

RUN R -e "devtools::install('/root/lipidr', dependencies=T)"