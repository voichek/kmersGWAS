FROM ubuntu:22.10
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# install Python3 and R
RUN apt update
RUN apt install -y python3
RUN apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#RUN sudo add-apt-repository -y 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt install -y r-base
RUN apt install -y python2
RUN apt install -y kmc
RUN apt install -y gemma
RUN apt install -y wget
RUN R -e 'install.packages("MASS")'
RUN R -e 'install.packages("mvnpermute")'
RUN R -e 'install.packages("matrixcalc")'

RUN wget https://github.com/voichek/kmersGWAS/releases/download/v0.2-beta/v0_2_beta.zip
RUN unzip v0_2_beta.zip
ENV PATH="/app/bin:${PATH}"
