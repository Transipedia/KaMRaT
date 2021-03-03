FROM ubuntu:20.04

ENV PATH="/usr/KaMRaT/bin:$PATH"
ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get clean all && apt-get update && apt-get -y upgrade && \
    apt-get -yq install build-essential cmake libboost-math-dev libboost-program-options-dev libboost-test-dev libboost-serialization-dev libarmadillo-dev binutils-dev python3-pandas python-numpy cython python-setuptools libensmallen-dev libstb-dev libmlpack-dev libboost-iostreams-dev && \
    apt -y install git && \
    cd /usr && git clone --recursive https://github.com/Transipedia/KaMRaT.git && cd KaMRaT && bash compile.bash && \
    echo "Creation Complete !"

ENTRYPOINT ["kamrat"]

