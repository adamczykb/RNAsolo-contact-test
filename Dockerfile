FROM debian:trixie
ARG DEBIAN_FRONTEND=noninteractive
WORKDIR /opt

RUN apt-get -y  update && apt-get install -y wget python3 python3-dev python3-pip git

# FINGERRNA
RUN apt-get --no-install-recommends -y install openbabel python3-openbabel python3-pip python-is-python3 \
    python3-pandas python3-numpy python3-tqdm python3-yaml cmake g++ gcc unzip pymol vim

WORKDIR /opt

RUN git clone --depth=1 https://github.com/n-szulc/fingernat.git


RUN git clone https://github.com/mhekkel/libmcfp.git 
WORKDIR /opt/libmcfp
RUN cmake -S . -B build
RUN cmake --build build
RUN cmake --install build

WORKDIR /opt

RUN git clone https://github.com/PDB-REDO/libcifpp.git 
WORKDIR /opt/libcifpp
RUN cmake -S . -B build
RUN cmake --build build
RUN cmake --install build

WORKDIR /opt
COPY mkdssp /usr/local/bin/mkdssp
# RUN git clone https://github.com/PDB-REDO/dssp.git 
# WORKDIR /opt/dssp
# RUN cmake -S . -B build
# RUN cmake --build build
# RUN cmake --install build

WORKDIR /opt

COPY . /opt
RUN unzip hbplus.zip
WORKDIR /opt/hbplus
RUN make
WORKDIR /opt/
RUN cp /opt/hbplus/hbplus /usr/local/bin/hbplus  
RUN cp /opt/fingernat/code/ /usr/local/bin/fingeRNAt -r
RUN pip3 install -r requirements.txt --break-system-packages

ENTRYPOINT ["tail", "-f", "/dev/null"]
