FROM augustus:latest

# env vars
ENV NXF_OPTS "-Xms1G -Xmx50G"
ENV JAVA_OPTS "-Xms4G -Xmx50G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
ENV TMPDIR /tmp

# apt stuff
# augustus is 3.0 in apt, 3.3.3 in github
RUN apt-get update && apt-get install -qyy curl git make vim cmake \
    gcc g++ unzip maven subversion gzip openjdk-8-jdk groovy wget \
    zlib1g-dev gnuplot lynx libncurses5-dev libncursesw5-dev libbz2-dev \
    liblzma-dev python python3-pip python-setuptools cython3 tabix bwa  \
  && apt-get clean 

# install stuff
RUN cd /opt  && mkdir -p /opt/bin \
  && cd /opt/bin && curl -fsSL get.nextflow.io | /bin/bash \
  && cd /opt/bin && wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/table2asn_GFF/linux64.table2asn_GFF.gz \
  && gunzip linux64.table2asn_GFF.gz && chmod 755 linux64.table2asn_GFF \
  && cd /opt/bin && wget https://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/scripts/gbf2tbl.pl \
  && chmod 755 gbf2tbl.pl \
  && cd /opt && wget https://sourceforge.net/projects/bbmap/files/latest/download \
  && mv download BBMap.tar && tar -xvzf BBMap.tar && rm BBMap.tar \
  && wget https://github.com/refresh-bio/KMC/releases/download/v3.1.1/KMC3.1.1.linux.tar.gz \
  && gunzip KMC3.1.1.linux.tar.gz && tar -xvf KMC3.1.1.linux.tar && rm -f KMC3.1.1.linux.tar \  
  && cd /opt && wget https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Linux-amd64.tar.xz && tar -xJf canu-2.0.*.tar.xz && rm canu-2.0*.xz \
#  && cd /opt && wget https://gite.lirmm.fr/lorma/lorma-releases/uploads/219b51b0d8d6ce378650743dc5f09024/lorma-bin_0.5_linux64.tar.gz \
#  && gunzip lorma-bin_0.5_linux64.tar.gz && tar -xvf lorma-bin_0.5_linux64.tar && rm lorma-bin_0.5_linux64.tar \
  && cd /opt && wget https://gite.lirmm.fr/lordec/lordec-releases/uploads/710113d83c210b6989ccfbdbafa89234/lordec-bin_0.9_linux64.tar.bz2 \
  && bunzip2 lordec-bin_0.9_linux64.tar.bz2 && tar -xvf lordec-bin_0.9_linux64.tar && rm lordec-bin_0.9_linux64.tar \
  && cd /opt && wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2 \
  && bunzip2 samtools-1.10.tar.bz2 && tar -xvf samtools-1.10.tar && rm samtools-1.10.tar \
  && cd samtools-1.10 && ./configure --prefix=/opt && make && make install && cd /opt && rm -rf samtools* \
  && cd /opt && wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip \
  && unzip bowtie2-2.3.5.1-linux-x86_64.zip && rm bowtie2-2.3.5.1-linux-x86_64.zip \
  && cd /opt/ && wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip \
  && unzip qualimap_v2.2.1.zip && rm qualimap_v2.2.1.zip \
  && cd /opt && wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
  && unzip fastqc_v0.11.9.zip && rm fastqc_v0.11.9.zip \
  && chmod 750 /opt/FastQC/fastqc \
  && wget https://sourceforge.net/projects/quast/files/quast-5.0.2.tar.gz \
  && tar -zxvf quast-5.0.2.tar.gz && cd quast-5.0.2 && ./setup.py install && cd .. \
  && rm quast-5.0.2.tar.gz

RUN pip3 install --upgrade NanoPlot

# google guava
RUN cd /opt \
  && git clone https://github.com/google/guava.git \
  && cd guava/guava \
  && mvn install

# kpi
#RUN cd /opt && git clone https://github.com/droeatumn/kass.git
RUN mkdir -p /opt/kass/input /opt/kass/output /opt/kass/bin /opt/kass/work /opt/kass/src
# todo: remove the 'add' when switching to git
ADD *.nf /opt/kass/
ADD input /opt/kass/input/
ADD bin /opt/kass/bin/
ADD src /opt/kass/src/
ENV CLASSPATH /opt/kass/bin/jars/slf4j-api-1.7.5.jar:/opt/kass/bin/jars/biojava4-core.jar:/opt/kass/bin/jars/dsh-commandline-1.1.jar:/opt/kass/bin/jars/super-csv.jar:$CLASSPATH
ENV CLASSPATH /opt/guava/guava/target/guava-HEAD-jre-SNAPSHOT.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH
ENV CLASSPATH /opt/jars/guava-21.0.jar:$CLASSPATH

# environment variables
ENV PATH /opt/bin:$PATH
#ENV PATH /opt/lorma-bin_0.5_linux64:$PATH
ENV PATH /opt/kass/bin/lorma-bin_0.5_linux64:$PATH
ENV PATH /opt/lordec-bin_0.9_linux64:$PATH
ENV PATH /opt/bowtie2-2.3.5.1-linux-x86_64:$PATH
ENV PATH /opt/qualimap_v2.2.1:$PATH
ENV PATH /opt/bbmap:$PATH
ENV PATH /opt/canu-2.0/Linux-amd64/bin:$PATH
ENV PATH /opt/FastQC:$PATH
ENV PATH /root/miniconda2/bin:$PATH
# docker
ENV AUGUSTUS_CONFIG_PATH /root/augustus/config
ENV PATH /root/augustus/bin:$PATH
ENV PATH /root/augustus/config:$PATH
#apt ENV PATH /usr/share/augustus/scripts:$PATH
#apt ENV AUGUSTUS_CONFIG_PATH /usr/share/augustus/config
ENV TMPDIR=/opt/kass/work
ENV TMP=/opt/kass/work

ENV PATH /opt/kass:$PATH
ENV PATH /opt/kass/bin:$PATH
ENV PATH /opt/kass/src:$PATH
CMD ["/opt/kass/main.nf"]
#CMD ["/bin/bash"]
