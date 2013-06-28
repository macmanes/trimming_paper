#!/usr/bin/make

all: out_1.fastq.quality out_1.fastq raw.Trinity.fasta real.0.Trinity.fasta real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta

out_1.fastq.quality:out_1.fastq
	perl $solexaqa/SolexaQA.pl -p 0.01 out_1.fastq
	cp out_1.fastq.quality ~/Dropbox/
	cp out_1.fastq.quality.pdf ~/Dropbox/
	
raw.Trinity.fasta:out_1.fastq out_2.fastq	
	$trinity/Trinity.pl --full_cleanup --seqType fq --JM 30G --left out_1.fastq  --right out_2.fastq  --CPU 8 --output raw

real.0.Trinity.fasta real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta: out_1.fastq out_2.fastq
	for TRIM in 0 1 2 5 10 15 20; do \
		java -Xmx30g -jar ~/software/trimmomatic-0.30.jar PE \
		-phred33 -threads 12 \
		out_1.fastq \
		out_2.fastq \
		T.$$TRIM.pp.1.fq \
		T.$$TRIM.up.1.fq \
		T.$$TRIM.pp.2.fq \
		T.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		perl /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1/SolexaQA.pl -p 0.01 T.$$TRIM.pp.1.fq ; \
		cp T.$$TRIM.pp.1.fq.quality ~/Dropbox/ ; \
		cp T.$$TRIM.pp.1.fq.quality.pdf ~/Dropbox/ ;\
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > right.$$TRIM.fq ; \
		/home/macmanes/trinityrnaseq_r2013-02-25/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.$$TRIM.fq --right right.$$TRIM.fq --CPU 8 --output real.$$TRIM ; \
	done