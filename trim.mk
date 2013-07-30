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

## Need to add 10M.Trinity.fasta, run.Trinity.fasta and SRR449363_1.fastq.quality back into all when its working fine.
all: 10M.$(READ1) 10M.$(READ2)  \
	10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq \
	10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq \
	10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta


#SRR449363_1.fastq.quality:SRR449363_1.fastq
#	perl $(SOLEXA)/SolexaQA.pl -p 0.01 SRR449363_1.fastq
#	cp SRR449363_1.fastq.quality ~/Dropbox/
#	cp SRR449363_1.fastq.quality ~/Dropbox/
#run.Trinity.fasta:SRR449363_1.fastq SRR449363_2.fastq
#	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
#	--left $(RUN) --right $(RUN) --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)


10M.$(READ1) 10M.$(READ2) : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 10000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)		

10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq: 10M.$(READ1) 10M.$(READ2)
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
		rm 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq ; \
	done

10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta: 10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 10M.left.$$TRIM.fq --right 10M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 10M.$$TRIM; \
	done
	
#pslx: 
#	for i in `*Trinity.fasta`; do \
#		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done

#orf: 
#	for i in `*Trinity.fasta`; do \
#		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done



















