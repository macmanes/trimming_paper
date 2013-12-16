#!/usr/bin/make -rRsf

all:eXpress trimmomatic Pfam GRCm38
clean_all:clean clean_check
	
eXpress:	
	wget http://bio.math.berkeley.edu/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz
	tar -zxf express-1.5.1-linux_x86_64.tgz
	cp express-1.5.1-linux_x86_64/express .
	rm express-1.5.1-linux_x86_64.tgz
	mkdir plugins/
	mv express-1.5.1-linux_x86_64/ plugins/
trimmomatic:
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
	unzip Trimmomatic-0.32.zip
	cp Trimmomatic-0.32/trimmomatic-0.32.jar .
	rm Trimmomatic-0.32.zip
	mv Trimmomatic-0.32/ plugins/
Pfam:
	wget http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-x86_64.tar.gz
	tar -zxf hmmer-3.1b1-linux-intel-x86_64.tar.gz
	cd hmmer-3.1b1-linux-intel-x86_64
	./configure && $(MAKE -j4)
	cp src/hmmpress ../.
	cp src/hmmscan ../.
	cp src/hmmsearch ../.
	cp src/hmmstat ../.
	wget http://sourceforge.net/projects/trinotate/files/TRINOTATE_RESOURCES/Pfam-A.hmm.gz
	gzip -d Pfam-A.hmm.gz
	hmmpress Pfam-A.hmm
GRCm38:
	wget ftp://ftp.ensembl.org/pub/release-71/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.71.cdna.all.fa.gz
	gzip -d Mus_musculus.GRCm38.71.cdna.all.fa.gz
clean:
	rm trimmomatic-0.32.jar express
	rm -fr plugins/*
clean_check:
	rm -fr unit_test/*