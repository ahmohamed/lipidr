FROM bioconductor/bioconductor_docker:devel
RUN R -e "devtools::install_github('ahmohamed/lipidr', dependencies=T)"
