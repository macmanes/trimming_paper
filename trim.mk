#!/usr/bin/make

SOLEXA=/home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1
TRINITY=/home/macmanes/trinityrnaseq_r2013-02-25
TRIMMOMATIC=/home/macmanes/software

real: out_1.fastq.quality out_1.fastq raw.Trinity.fasta right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta
sim: right.fastq.quality sim.Trinity.fasta sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq sim.1.Trinity.fasta sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta sim.20.Trinity.fasta

out_1.fastq.quality:out_1.fastq
	perl $SOLEXA/SolexaQA.pl -p 0.01 out_1.fastq
	cp out_1.fastq.quality ~/Dropbox/
	cp out_1.fastq.quality.pdf ~/Dropbox/
	
raw.Trinity.fasta:out_1.fastq
	$trinity/Trinity.pl --full_cleanup --seqType fq --JM 30G --left out_1.fastq  --right out_2.fastq  --CPU 8 --output raw

right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq: out_1.fastq
	for TRIM in 1 2 5 10 15 20; do \
		java -Xmx30g -jar $$(TRIMMMOMATIC)/trimmomatic-0.30.jar PE \
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
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta: out_1.fastq
	for TRIM in 1 2 5 10 15 20; do \
		/home/macmanes/trinityrnaseq_r2013-02-25/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.$$TRIM.fq --right right.$$TRIM.fq --CPU 8 --output real.$$TRIM ; \
	done
	
##Make sim	
	
right.fastq.quality:right.fastq
	perl $solexaqa/SolexaQA.pl -p 0.01 right.fastq
	cp right.fastq.quality ~/Dropbox/
	cp right.fastq.quality.pdf ~/Dropbox/
	
sim.Trinity.fasta:out_1.fastq
	$trinity/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.fastq  --right right.fastq  --CPU 8 --output sim

sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq: right.fastq
	for TRIM in 1 2 5 10 15 20; do \
		java -Xmx30g -jar ~/software/trimmomatic-0.30.jar PE \
		-phred33 -threads 12 \
		left.fastq \
		right.fastq \
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
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > sim.left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > sim.right.$$TRIM.fq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

sim.1.Trinity.fasta sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta sim.20.Trinity.fasta: out_1.fastq
	for TRIM in 1 2 5 10 15 20; do \
		/home/macmanes/trinityrnaseq_r2013-02-25/Trinity.pl --full_cleanup --seqType fq --JM 30G --left sim.left.$$TRIM.fq --right sim.right.$$TRIM.fq --CPU 8 --output sim.$$TRIM ; \
	done
