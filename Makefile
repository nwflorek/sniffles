# Sniffles
#Author: Nick Florek
#Performs SNP analysis of influenza genomes from raw reads.

PREFIX := $(shell pwd)
TMPDIR := build

# Determine OS
OS := $(shell uname -s)

default: install

help:
	@echo "Please see the Readme for additional help"

install: install-prerequisites
	@echo "Installation complete!"

##TODO add checking to see if these programs exist on the filesystem already
install-prerequisites: install-mkdir install-vcftools install-picard install-snpeff install-trimmomatic install-varscan \
	install-samtools install-bowtie2 install-samtools install-lofreq install-trinity install-bbmap \
	install-blast install-popoolation cleanup
	@echo "Done installing the prerequisites."

install-mkdir:
	-mkdir build lib lib/bin

install-picard:
	mkdir lib/picard
	curl -L 'https://github.com/broadinstitute/picard/releases/download/2.18.10/picard.jar' -o lib/picard/picard.jar

install-snpeff:
	mkdir lib/snpeff
	curl -L 'https://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip' -o $(TMPDIR)/snpEff_latest_core.zip
	unzip $(TMPDIR)/snpEff_latest_core.zip -d $(TMPDIR)/snpEff_latest_core
	mv $(TMPDIR)/snpEff_latest_core/snpEff/* lib/snpeff/

install-trimmomatic:
	curl -L 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip' -o $(TMPDIR)/Trimmomatic-0.36.zip
	unzip $(TMPDIR)/Trimmomatic-0.36.zip -d $(TMPDIR)
	mkdir lib/trimmomatic
	mv $(TMPDIR)/Trimmomatic-0.36/trimmomatic-0.36.jar lib/trimmomatic/trimmomatic-0.36.jar
	mv $(TMPDIR)/Trimmomatic-0.36/adapters lib/trimmomatic/

install-varscan:
	mkdir lib/varscan
	curl -L 'https://downloads.sourceforge.net/project/varscan/VarScan.v2.3.9.jar' -o lib/varscan/VarScan.v2.3.9.jar

install-bowtie2:
	curl -L 'https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip' -o $(TMPDIR)/bowtie2.3.4.1-source.zip
	unzip $(TMPDIR)/bowtie2.3.4.1-source.zip -d $(TMPDIR)
	cd $(TMPDIR)/bowtie2-2.3.4.1; $(MAKE) NO_TBB=1; $(MAKE) prefix=../../lib/ install

install-lofreq:
	@if [ "$(OS)" = "Linux" ]; then\
		curl -L 'https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.3.1_linux-x86-64.tgz' -o $(TMPDIR)/lofreq_star-2.1.3.1.tgz;\
	fi
	@if [ "$(OS)" = "Darwin" ]; then\
		curl -L 'https://github.com/CSB5/lofreq/raw/master/dist/lofreq_star-2.1.3.1_macosx.tgz' -o $(TMPDIR)/lofreq_star-2.1.3.1.tgz;\
	fi
	tar -xzf $(TMPDIR)/lofreq_star-2.1.3.1.tgz -C $(TMPDIR)
	mv $(TMPDIR)/lofreq_star-2.1.3.1/bin/* lib/bin/

install-samtools:
	curl -L 'https://downloads.sourceforge.net/project/samtools/samtools/1.9/samtools-1.9.tar.bz2' -o $(TMPDIR)/samtools-1.9.tar.bz2
	cd $(TMPDIR); tar jxvf samtools-1.9.tar.bz2
	cd $(TMPDIR)/samtools-1.9/htslib-1.9; $(MAKE); $(MAKE) prefix=../../lib/ install
	cd $(TMPDIR)/samtools-1.9; ./configure;\
	 #make DFLAGS="-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=0" LIBCURSES="";\
	 $(MAKE) prefix=../../lib/ install
	cp $(TMPDIR)/samtools-1.9/htslib-1.9/bgzip lib/bin/bgzip
	cp $(TMPDIR)/samtools-1.9/htslib-1.9/tabix lib/bin/tabix

install-trinity:
	curl -L 'https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.6.6.tar.gz' -o $(TMPDIR)/Trinity-v2.6.6.tar.gz
	tar -xzf $(TMPDIR)/Trinity-v2.6.6.tar.gz -C $(TMPDIR)
	cd $(TMPDIR)/trinityrnaseq-Trinity-v2.6.6; $(MAKE); $(MAKE) plugins
	mv $(TMPDIR)/trinityrnaseq-Trinity-v2.6.6/ lib/trinity-2.6.6/

install-bbmap:
	curl -L 'https://downloads.sourceforge.net/project/bbmap/BBMap_38.16.tar.gz' -o $(TMPDIR)/BBMap_38.16.tar.gz
	tar -xzf $(TMPDIR)/BBMap_38.16.tar.gz -C $(TMPDIR)
	mv $(TMPDIR)/bbmap lib/bbmap

install-blast:
	@if [ "$(OS)" = "Linux" ]; then\
		curl -L 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz' -o $(TMPDIR)/ncbi-blast.tar.gz;\
	fi
	@if [ "$(OS)" = "Darwin" ]; then\
		curl -L 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-macosx.tar.gz' -o $(TMPDIR)/ncbi-blast.tar.gz;\
	fi
	tar -xzf $(TMPDIR)/ncbi-blast.tar.gz -C $(TMPDIR)
	mv $(TMPDIR)/ncbi-blast-2.7.1+/bin/* lib/bin/

install-popoolation:
	curl -L 'https://downloads.sourceforge.net/project/popoolation/popoolation_1.2.2.zip' -o $(TMPDIR)/popoolation_1.2.2.zip
	unzip $(TMPDIR)/popoolation_1.2.2.zip -d $(TMPDIR)
	mv $(TMPDIR)/popoolation_1.2.2 lib/poppoolation

install-vcftools:
	curl -L 'https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz' -o $(TMPDIR)/vcftools.tar.gz
	tar -xzf $(TMPDIR)/vcftools.tar.gz -C $(TMPDIR)
	cd $(TMPDIR)/vcftools-0.1.16; ./configure --prefix=$(PWD)/lib
	$(MAKE) -C $(TMPDIR)/vcftools-0.1.16
	$(MAKE) -C $(TMPDIR)/vcftools-0.1.16 install

cleanup:
	rm -r build
