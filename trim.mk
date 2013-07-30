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

.PHONY: check

check:
	@echo "\n\n\n"###I am checking to see if you have all the dependancies installed.### "\n"
	command -v trimmomatic-0.30.jar >/dev/null 2>&1 || { echo >&2 "I require Trimmomatic but it's not installed.  Aborting."; exit 1; }
	@echo Trimmomatic is Installed
	command -v $(TRINITY/Trinity.pl) >/dev/null 2>&1 || { echo >&2 "I require Trinity but it's not installed.  Aborting."; exit 1; }
	@echo Trinity is Installed
	command -v fastq-converter-v2.0.pl >/dev/null 2>&1 || { echo >&2 "I require fastq-converter-v2.0.pl (Reptile package) but it's not installed.  Aborting."; exit 1; }
	command -v reptile-omp >/dev/null 2>&1 || { echo >&2 "I require reptile-omp but it's not installed.  Aborting."; exit 1; }
	command -v reptile_merger >/dev/null 2>&1 || { echo >&2 "I require reptile_merger but it's not installed.  Aborting."; exit 1; }
	@echo Reptile is installed"\n"




all: 10M.SRR449363_1.fastq 10M.SRR449363_2.fastq 10M.Trinity.fasta 10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq \
	10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta \
	right.fq.quality 20M.Trinity.fasta 20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq \
	20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta \
	

#sub: SRR449363_1.fastq SRR449363_1.fastq
#10M: right.fq.quality 10M.Trinity.fasta 10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq \
#	10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta 
#20M: right.fq.quality 20M.Trinity.fasta 20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq \
#	20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta 
#pslx: 
#orf : 

10M.SRR449363_1.fastq : SRR449363_1.fastq
	python ~/error_correction/scripts/subsampler.py 10000000 SRR449363_1.fastq SRR449363_2.fastq
20M.SRR449363_1.fastq : SRR449363_1.fastq
	python ~/error_correction/scripts/subsampler.py 20000000 SRR449363_1.fastq SRR449363_2.fastq

SRR449363_1.fastq.quality:SRR449363_1.fastq
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 SRR449363_1.fastq
	cp SRR449363_1.fastq.quality ~/Dropbox/
	cp SRR449363_1.fastq.quality ~/Dropbox/
	
SRR449363_1.Trinity.fasta:SRR449363_1.fastq SRR449363_2.fastq
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left $(RUN) --right $(RUN) --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)

right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq: out_1.fastq
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
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
		perl $(SOLEXA)/SolexaQA.pl -p 0.01 T.$$TRIM.pp.1.fq ; \
		cp T.$$TRIM.pp.1.fq.quality ~/Dropbox/ ; \
		cp T.$$TRIM.pp.1.fq.quality.pdf ~/Dropbox/ ;\
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > right.$$TRIM.fq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta: out_1.fastq
	for TRIM in 1 2 5 10 15 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.$$TRIM.fq --right right.$$TRIM.fq --CPU $(CPU) --output real.$$TRIM ; \
	done
	
##Make sim	
	
right.fq.quality:right.fq
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 right.fq
	cp right.fq.quality ~/Dropbox/
	cp right.fq.quality.pdf ~/Dropbox/
	
sim.Trinity.fasta:right.fq.fastq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.fq  --right right.fq  --CPU $(CPU) --output sim

pslx: 
	for i in `*Trinity.fasta`; do \
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done

orf: 
	for i in `*Trinity.fasta`; do \
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done



















