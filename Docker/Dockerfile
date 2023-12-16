FROM ubuntu:latest

ENV PATH="/usr/KaMRaT/apps/:$PATH"
ENV LC_ALL=C

RUN export DEBIAN_FRONTEND=noninteractive \
	&& apt-get clean all && apt-get update && apt-get -y upgrade \
	&& apt-get -yq install build-essential cmake libboost-math-dev libboost-program-options-dev libboost-test-dev libboost-serialization-dev libarmadillo-dev binutils-dev python3-pandas python3-numpy python3-dev python-setuptools libensmallen-dev libstb-dev liblapack-dev libmlpack-dev libblas-dev libboost-dev libboost-iostreams-dev libarmadillo-dev pigz \
	&& apt -y install git jellyfish \
	&& pip install snakemake \
	&& cd /usr && git clone --recursive https://github.com/Transipedia/KaMRaT.git && cd KaMRaT \
	&& cmake . && make -j
RUN cd /usr/KaMRaT/related-tools/prepare_kmer_table/dekupl-joinCounts \
	&& make \
	&& mv joinCounts /usr/KaMRaT/apps/
RUN cp /usr/KaMRaT/related-tools/prepare_kmer_table/*.bash /usr/KaMRaT/apps/ \
	&& cp /usr/KaMRaT/related-tools/prepare_kmer_table/*.pl /usr/KaMRaT/apps/ \
	&& echo "Creation Complete !"
