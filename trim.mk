#!/usr/bin/make -rRsf

MEM=5
TRIM=5
CPU=5
RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')
BCODES := /home/macmanes/Dropbox/barcodes.fa
SOLEXA := /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1
TRINITY := /home/macmanes/trinityrnaseq-code/trunk
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
PFAM := /media/macmanes/raid/blastdb/Pfam-A.hmm


all: subsamp trim trim1 trim2 trim3 trim4 trin trin1 trin2 trin3 trin4 pslx orf


$(READ1).quality:
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 $(READ1)
	cp $(READ1).quality ~/Dropbox/
	cp $(READ1).quality.pdf ~/Dropbox/
#run.Trinity.fasta:SRR449363_1.fastq SRR449363_2.fastq
#	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
#	--left $(RUN) --right $(RUN) --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)


subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 1000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)		
subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 10000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)	
subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 20000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.20M.$(READ1)
	mv subsamp_2.fastq raw.20M.$(READ2)	
subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 50000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.50M.$(READ1)
	mv subsamp_2.fastq raw.50M.$(READ2)	
subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 75000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.75M.$(READ1)
	mv subsamp_2.fastq raw.75M.$(READ2)	
subsamp : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 100000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.100M.$(READ1)
	mv subsamp_2.fastq raw.100M.$(READ2)	

trim: 
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
		raw.10M.$(READ1) \
		raw.10M.$(READ2) \
		10M.$$TRIM.pp.1.fq \
		10M.$$TRIM.up.1.fq \
		10M.$$TRIM.pp.2.fq \
		10M.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		cat 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq > 10M.left.$$TRIM.fq ; \
		cat 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq > 10M.right.$$TRIM.fq ; \
		rm 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq ; done
trin: 
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 10M.left.$$TRIM.fq --right 10M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 10M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done

trim1: 
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
		raw.20M.$(READ1) \
		raw.20M.$(READ2) \
		20M.$$TRIM.pp.1.fq \
		20M.$$TRIM.up.1.fq \
		20M.$$TRIM.pp.2.fq \
		20M.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		cat 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq > 20M.left.$$TRIM.fq ; \
		cat 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq > 20M.right.$$TRIM.fq ; \
		rm 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq ; done
trin1: 
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 20M.left.$$TRIM.fq --right 20M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 20M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done

trim2: 
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
		raw.50M.$(READ1) \
		raw.50M.$(READ2) \
		50M.$$TRIM.pp.1.fq \
		50M.$$TRIM.up.1.fq \
		50M.$$TRIM.pp.2.fq \
		50M.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		cat 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq > 50M.left.$$TRIM.fq ; \
		cat 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq > 50M.right.$$TRIM.fq ; \
		rm 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq ; done
trin2: 
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 50M.left.$$TRIM.fq --right 50M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 50M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done

trim3: 
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
		raw.75M.$(READ1) \
		raw.75M.$(READ2) \
		75M.$$TRIM.pp.1.fq \
		75M.$$TRIM.up.1.fq \
		75M.$$TRIM.pp.2.fq \
		75M.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		cat 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq > 75M.left.$$TRIM.fq ; \
		cat 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq > 75M.right.$$TRIM.fq ; \
		rm 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq ; done
trin3: 
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 75M.left.$$TRIM.fq --right 75M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 75M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done
trim4: 
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
		raw.100M.$(READ1) \
		raw.100M.$(READ2) \
		100M.$$TRIM.pp.1.fq \
		100M.$$TRIM.up.1.fq \
		100M.$$TRIM.pp.2.fq \
		100M.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		cat 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq > 100M.left.$$TRIM.fq ; \
		cat 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq > 100M.right.$$TRIM.fq ; \
		rm 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq ; done
trin4: 
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 100M.left.$$TRIM.fq --right 100M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 100M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done

pslx: 
	for i in `*Trinity.fasta`; do \
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done

orf: 
	for i in `*Trinity.fasta`; do \
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done



















