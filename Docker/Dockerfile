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
    && cmake . \
    && make -j
RUN cd /usr/KaMRaT/related-tools/prepare_kmer_table/dekupl-joinCounts \
    && make \
    && mv joinCounts /usr/KaMRaT/apps/
RUN cp /usr/KaMRaT/related-tools/prepare_kmer_table/*.bash /usr/KaMRaT/apps/ \
    && cp /usr/KaMRaT/related-tools/prepare_kmer_table/*.pl /usr/KaMRaT/apps/
RUN apt-get purge -y --auto-remove $buildDeps git \
    && echo "Creation Complete !"
