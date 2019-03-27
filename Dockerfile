#Use ubuntu as base image
FROM ubuntu:xenial

# metadata
LABEL base.image="ubuntu:xenial"
LABEL version="1"
LABEL software="sniffles"
LABEL software.version="1.0"
LABEL description="influenza snp detection pipeline"
LABEL website="https://github.com/nwflorek/sniffles"

# Maintainer
MAINTAINER Nicholas Florek <nicholas.florek@slh.wisc.edu>

#install python
RUN apt-get update && apt-get install -y software-properties-common &&\
add-apt-repository -y ppa:openjdk-r/ppa && apt-get update && apt-get install -y \
 curl\
 python\
 build-essential\
 openjdk-11-jre\
 zip\
 libz-dev\
 libbz2-dev\ 
 liblzma-dev\
 libncurses5-dev

ENV PATH="/tools:${PATH}"

#install picard and varscan
RUN mkdir /tools; curl -L 'https://github.com/broadinstitute/picard/releases/download/2.18.27/picard.jar' -o /tools/picard.jar &&\
 curl -L 'https://downloads.sourceforge.net/project/varscan/VarScan.v2.3.9.jar' -o /tools/varscan.jar

#install trimmomatic
RUN curl -L 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip' -o Trimmomatic-0.38.zip &&\
 unzip Trimmomatic-0.38.zip &&\
 mv Trimmomatic-0.38/trimmomatic-0.38.jar /tools/trimmomatic.jar &&\
 mv Trimmomatic-0.38/adapters /tools/adapters &&\
 rm -r Trimmomatic-0.38 Trimmomatic-0.38.zip

#install bowtie2
RUN curl -L 'https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip' -o bowtie2-2.3.4.1-source.zip &&\
 unzip bowtie2-2.3.4.1-source.zip &&\
 cd bowtie2-2.3.4.1; make NO_TBB=1 && make install &&\
 rm -r /bowtie2-2.3.4.1 /bowtie2-2.3.4.1-source.zip

#install samtools
RUN curl -L 'https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2' -o samtools-1.9.tar.bz2 &&\
 tar jxvf samtools-1.9.tar.bz2 &&\
 cd samtools-1.9/htslib-1.9 &&\
 make &&\
 make install &&\
 cd ..; ./configure &&\
 make install &&\
 rm -r /samtools-1.9.tar.bz2 /samtools-1.9

#install bcftools
RUN curl -L 'https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2' -o bcftools-1.9.tar.bz2 &&\
 tar jxvf bcftools-1.9.tar.bz2 &&\
 cd bcftools-1.9 &&\
 ./configure &&\
 make install &&\
 rm -r /bcftools-1.9.tar.bz2 /bcftools-1.9

#install ncbi-blast
RUN curl -L 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz' -o ncbi-blast.tar.gz &&\
 tar -xzf ncbi-blast.tar.gz &&\
 mv ncbi-blast-2.7.1+/bin/* usr/bin &&\
 rm -r /ncbi-blast.tar.gz /ncbi-blast-2.7.1+/

#install bbtools
RUN curl -L 'https://downloads.sourceforge.net/project/bbmap/BBMap_38.43.tar.gz' -o BBMap_38.43.tar.gz &&\
 tar -xzf BBMap_38.43.tar.gz &&\
 mv bbmap /tools/bbmap &&\
 rm -r BBMap_38.43.tar.gz

#install lofreq
RUN curl -L 'https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.3.1_linux-x86-64.tgz' -o lofreq_star-2.1.3.1_linux-x86-64.tgz &&\
 tar -xzf lofreq_star-2.1.3.1_linux-x86-64.tgz &&\
 cd lofreq_star-2.1.3.1 && cd bin &&\
 mv * /usr/bin &&\
 rm -r /lofreq_star-2.1.3.1_linux-x86-64.tgz /lofreq_star-2.1.3.1

WORKDIR /data
