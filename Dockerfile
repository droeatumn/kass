FROM ubuntu:latest
MAINTAINER Dave Roe

# env vars
ENV NXF_OPTS="-Xms1G -Xmx200G"
ENV JAVA_OPTS="-Xms4G -Xmx200G"
ENV LD_LIBRARY_PATH=/opt/lib:$LD_LIBRARY_PATH
ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64/
ENV TMPDIR=/tmp

# apt stuff
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
    gcc g++ zip unzip maven subversion gzip openjdk-17-jdk openjdk-17-doc wget \
    zlib1g-dev gnuplot lynx libncurses5-dev libncursesw5-dev libbz2-dev \
    liblzma-dev python3-pip python-setuptools tabix bwa \
  && apt-get clean 

# install groovy with sdkman
# https://stackoverflow.com/questions/53656537/install-sdkman-in-docker-image
RUN mv /bin/sh /bin/sh.bak && ln -s /bin/bash /bin/sh
ENV SDKMAN_DIR /usr/local/sdkman
RUN set -x \
    && apt-get update \
    && apt-get install -y unzip --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*
RUN curl -s get.sdkman.io | bash \
    && chmod 700 /usr/local/sdkman/bin/sdkman-init.sh 
RUN source "/usr/local/sdkman/bin/sdkman-init.sh" \
    && sdk install groovy  \
    && rm -f /bin/sh && mv /bin/sh.bak /bin/sh
ENV GROOVY_HOME /usr/local/sdkman/candidates/groovy/current/ \
    && ENV PATH $GROOVY_HOME/bin:$PATH

# install other stuff
RUN cd /opt  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
#droe(remove?)  && cd /opt/bin && wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/linux64.table2asn.gz \
#droe(remove?)  && gunzip linux64.table2asn.gz && ln -s linux64.table2asn linux64.table2asn_GFF && chmod 755 linux64.table2asn \
  && cd /opt && wget https://sourceforge.net/projects/bbmap/files/latest/download \
  && mv download BBMap.tar && tar -xvzf BBMap.tar && rm BBMap.tar \
#  && wget https://github.com/refresh-bio/KMC/releases/download/v3.1.1/KMC3.1.1.linux.tar.gz \
#  && gunzip KMC3.1.1.linux.tar.gz && tar -xvf KMC3.1.1.linux.tar && rm -f KMC3.1.1.linux.tar \  
  && cd /opt && wget https://gite.lirmm.fr/lorma/lorma-releases/uploads/219b51b0d8d6ce378650743dc5f09024/lorma-bin_0.5_linux64.tar.gz \
  && gunzip lorma-bin_0.5_linux64.tar.gz && tar -xvf lorma-bin_0.5_linux64.tar && rm lorma-bin_0.5_linux64.tar \
  && cd /opt && wget https://gite.lirmm.fr/lordec/lordec-releases/uploads/710113d83c210b6989ccfbdbafa89234/lordec-bin_0.9_linux64.tar.bz2 \
  && bunzip2 lordec-bin_0.9_linux64.tar.bz2 && tar -xvf lordec-bin_0.9_linux64.tar && rm lordec-bin_0.9_linux64.tar \
  && cd /opt && wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
  && bunzip2 samtools-1.16.1.tar.bz2 && tar -xvf samtools-1.16.1.tar && rm samtools-1.16.1.tar \
  && cd samtools-1.16.1 && ./configure --prefix=/opt && make && make install && cd /opt && rm -rf samtools* \
  && cd /opt && wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip \
  && unzip bowtie2-2.4.5-linux-x86_64.zip && rm bowtie2-2.4.5-linux-x86_64.zip \
  && cd /opt/ && wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip \
  && unzip qualimap_v2.2.1.zip && rm qualimap_v2.2.1.zip \
  && cd /opt && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
  && unzip fastqc_v0.11.9.zip && rm fastqc_v0.11.9.zip \
  && chmod 750 /opt/FastQC/fastqc \
#  && wget https://sourceforge.net/projects/quast/files/quast-5.0.2.tar.gz \
#  && tar -zxvf quast-5.0.2.tar.gz && cd quast-5.0.2 && ./setup.py install && cd .. \
#  && rm quast-5.0.2.tar.gz \
  && wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 \
  && bunzip2 minimap2-2.24_x64-linux.tar.bz2 \
  && tar -xvf minimap2-2.24_x64-linux.tar && rm minimap2-2.24_x64-linux.tar
RUN mkdir -p /opt/jars && cd /opt/jars/ \
  && wget https://repo1.maven.org/maven2/org/slf4j/slf4j-nop/2.0.2/slf4j-nop-2.0.2.jar \
  && wget https://repo1.maven.org/maven2/org/slf4j/slf4j-api/2.0.2/slf4j-api-2.0.2.jar

# install assembler
# from Nurk et al. 2020 ("HiCanu ...")
#   "HiCanu was run using Canu branch hicanu_rc with the following commands: canu -assemble -p asm -d asm genomeSize=G -pacbio- hifi reads.fastq.gz"
# apparently incorporated into main branch before 2.2
RUN cd /opt && curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz \
  && tar -xJf canu-2.2.Linux.tar.xz && rm canu-2*.xz

RUN pip3 install numpy cython NanoPlot

# google guava
RUN mkdir -p /opt/jars && cd /opt/jars \
  && wget https://repo1.maven.org/maven2/com/google/guava/guava/31.1-jre/guava-31.1-jre.jar \
  && chmod 755 guava-31.1-jre.jar

# kass
RUN cd /opt && git clone https://github.com/droeatumn/kass.git
RUN chmod 755 /opt/kass/bin/* /opt/kass/src/*
ENV CLASSPATH /opt/kass/bin/jars/slf4j-api-2.0.2.jar:/opt/kass/bin/jars/dsh-commandline-1.1.jar:/opt/kass/bin/jars/super-csv.jar:/opt/jars/slf4j-nop-2.0.2.jar:$CLASSPATH
ENV CLASSPATH /opt/kass/bin/jars/biojava-alignment-5.4.0.jar:/opt/kass/bin/jars/biojava-core-5.4.0.jar:$CLASSPATH
ENV CLASSPATH /opt/jars/guava-31.1-jre.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH

# environment variables
ENV PATH /opt/bin:$PATH
ENV PATH /opt/lorma-bin_0.5_linux64:$PATH
ENV PATH /opt/lordec-bin_0.9_linux64:$PATH
ENV PATH /opt/bowtie2-2.4.5-linux-x86_64:$PATH
ENV PATH /opt/qualimap_v2.2.1:$PATH
ENV PATH /opt/bbmap:$PATH
ENV PATH /opt/canu-2.2/bin:$PATH
ENV PATH /opt/FastQC:$PATH
ENV PATH /root/miniconda2/bin:$PATH
ENV PATH /opt/minimap2-2.24_x64-linux:$PATH

# docker
ENV TMPDIR=/opt/kass/work
ENV TMP=/opt/kass/work

ENV PATH /opt/kass:$PATH
ENV PATH /opt/kass/bin:$PATH
ENV PATH /opt/kass/src:$PATH
CMD ["/opt/kass/main.nf"]
#CMD ["/bin/bash"]
