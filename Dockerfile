################## BASE IMAGE ######################
FROM nfcore/base

################## METADATA ######################

LABEL base_image="nfcore/base"
LABEL version="1.0"
LABEL software="rf-mut-f"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for **RF-mut-f**"
LABEL about.home="http://github.com/adigenova/RF-mut-f"
LABEL about.documentation="http://github.com/adigenova/RF-mut-f/README.md"
LABEL about.license_file="http://github.com/adigenova/RF-mut-f/LICENSE.txt"
LABEL about.license="MIT"

################## MAINTAINER ######################
MAINTAINER **digenovaa** <**digenovaac@fellows.iarc.fr**>

################## INSTALLATION ######################
COPY environment.yml /
RUN conda env update -n rf-mut-f -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rf-mut-f/bin:$PATH
RUN conda env export --name rf-mut-f > rf-mut-f-v1.0.yml
