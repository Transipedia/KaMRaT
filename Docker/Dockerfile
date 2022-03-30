FROM ubuntu:20.04

ENV PATH="/usr/KaMRaT/bin:$PATH"
ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive \
    && buildDeps="build-essential cmake libarmadillo-dev libensmallen-dev libmlpack-dev=3.2.2-3 libboost-iostreams-dev" \
    && apt-get update \
    && apt-get -yq install $buildDeps \
    && apt -y install git
RUN cd /usr \
    # && git clone -b release/v1.0 https://github.com/Transipedia/KaMRaT.git \
    # && cd KaMRaT/related-tools/prepare_kmer_table/dekupl-joinCounts/ \
    # && git submodule update --init --recursive
    && git clone -b release/v1.0 --recursive https://github.com/Transipedia/KaMRaT.git
RUN cd /usr/KaMRaT \
    && bash compile.bash \
    && apt-get purge -y --auto-remove $buildDeps git \
    && echo "Creation Complete !"

ENTRYPOINT ["kamrat"]
