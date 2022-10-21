FROM ubuntu:22.04

ENV PATH="/usr/KaMRaT/apps/:$PATH"
ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive \
    && buildDeps="build-essential cmake libarmadillo-dev libensmallen-dev libmlpack-dev libboost-iostreams-dev" \
    && apt-get update \
    && apt-get -yq install $buildDeps \
    && apt -y install git
RUN cd /usr \
    && git clone --recursive https://github.com/Transipedia/KaMRaT.git \
    && cd KaMRaT/ \
    && git checkout dev \
    && cmake . \
    && make -j
RUN apt-get purge -y --auto-remove $buildDeps git \
    && echo "Creation Complete !"

ENTRYPOINT ["kamrat"]
