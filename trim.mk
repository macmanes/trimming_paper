#!/usr/bin/make -rRsf

MEM=5
TRIM=5
CPU=5
RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')
BCODES := /home/macmanes/Dropbox/barcodes.fa
SOLEXA := /home/macmanes/software/
TRINITY := /home/macmanes/trinityrnaseq-code/trunk
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
PFAM := /media/macmanes/raid/blastdb/Pfam-A.hmm


all: subsamp1 trim trin subsamp2 trim1 trin1 trin1a subsamp3 trim2 trin2 subsamp4 trim3 trin3 subsamp5 trim4 trin4 pslx orf
quality: $(READ1).quality

subsamp10:raw.10M.$(READ1) raw.10M.$(READ2)
trim10:10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq
trin10:10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta raw.10M.Trinity.fasta
subsamp20:raw.20M.$(READ1) raw.20M.$(READ2)
trim20:20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq
trin20:20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta raw.20M.Trinity.fasta
subsamp50:raw.50M.$(READ1) raw.50M.$(READ2)
trim50:50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq
trin50:50M.2.Trinity.fasta 50M.5.Trinity.fasta 50M.10.Trinity.fasta 50M.20.Trinity.fasta raw.50M.Trinity.fasta
subsamp75:raw.75M.$(READ1) raw.75M.$(READ2)
trim75:75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq
trin75:75M.2.Trinity.fasta 75M.5.Trinity.fasta 75M.10.Trinity.fasta 75M.20.Trinity.fasta raw.75M.Trinity.fasta

subsamp100:raw.100M.$(READ1) raw.100M.$(READ2)
trim10:100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq
trin10:100M.2.Trinity.fasta 100M.5.Trinity.fasta 100M.10.Trinity.fasta 100M.20.Trinity.fasta raw.100M.Trinity.fasta


orf: 10M.2.Trinity.fasta.pep 10M.5.Trinity.fasta.pep 10M.10.Trinity.fasta.pep 10M.20.Trinity.fasta.pep 20M.2.Trinity.fasta.pep 20M.5.Trinity.fasta.pep 20M.10.Trinity.fasta.pep 20M.20.Trinity.fasta.pep 50M.2.Trinity.fasta.pep 50M.5.Trinity.fasta.pep 10M.50.Trinity.fasta.pep 50M.20.Trinity.fasta.pep 75M.2.Trinity.fasta.pep 75M.5.Trinity.fasta.pep 75M.10.Trinity.fasta.pep 75M.20.Trinity.fasta.pep 100M.2.Trinity.fasta.pep 100M.5.Trinity.fasta.pep 100M.10.Trinity.fasta.pep 100M.20.Trinity.fasta.pep

$(READ1).quality: $(READ1) $(READ2)
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 $(READ1)
	cp $(READ1).quality ~/Dropbox/
	cp $(READ1).quality.pdf ~/Dropbox/
#run.Trinity.fasta:SRR449363_1.fastq SRR449363_2.fastq
#	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
#	--left $(RUN) --right $(RUN) --group_pairs_distance 999 --CPU $(CPU) --output $(RUN)


raw.10M.$(READ1) raw.10M.$(READ2): $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 10000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)	

10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq: raw.10M.$(READ1) raw.10M.$(READ2)
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
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
		MINLEN:25 2>> trim10; \
		cat 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq > 10M.left.$$TRIM.fq ; \
		cat 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq > 10M.right.$$TRIM.fq ; \
		rm 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq ; done
10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta raw.10M.Trinity.fasta: 10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 10M.left.$$TRIM.fq --right 10M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 10M.$$TRIM; rm 10M.left.$$TRIM.fq 10M.right.$$TRIM.fq; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left raw.10M.$(READ1) --right raw.10M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.10M; rm raw.10M.$(READ1) raw.10M.$(READ2)

raw.20M.$(READ1) raw.20M.$(READ2): $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 20000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.20M.$(READ1)
	mv subsamp_2.fastq raw.20M.$(READ2)	
20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq: raw.20M.$(READ1) raw.20M.$(READ2)
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
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
		MINLEN:25 2>> trim20 ; \
		cat 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq > 20M.left.$$TRIM.fq ; \
		cat 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq > 20M.right.$$TRIM.fq ; \
		rm 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq ; done
20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta raw.20M.Trinity.fasta: 20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 20M.left.$$TRIM.fq --right 20M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 20M.$$TRIM; rm 20M.left.$$TRIM.fq 20M.right.$$TRIM.fq; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left raw.20M.$(READ1) --right raw.20M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.20M; rm raw.20M.$(READ1) raw.20M.$(READ2)

raw.50M.$(READ1) raw.50M.$(READ2): $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 50000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.50M.$(READ1)
	mv subsamp_2.fastq raw.50M.$(READ2)	
50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq: raw.50M.$(READ1) raw.50M.$(READ2)
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
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
		MINLEN:25 2>> trim50 ; \
		cat 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq > 50M.left.$$TRIM.fq ; \
		cat 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq > 50M.right.$$TRIM.fq ; \
		rm 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq ; done

50M.2.Trinity.fasta 50M.5.Trinity.fasta 50M.10.Trinity.fasta 50M.20.Trinity.fasta raw.50M.Trinity.fasta: 50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 50M.left.$$TRIM.fq --right 50M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 50M.$$TRIM; rm 50M.left.$$TRIM.fq 50M.right.$$TRIM.fq; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left raw.50M.$(READ1) --right raw.50M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.50M; rm raw.50M.$(READ1) raw.50M.$(READ2)

raw.75M.$(READ1) raw.75M.$(READ2) : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 75000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.75M.$(READ1)
	mv subsamp_2.fastq raw.75M.$(READ2)	
75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq: $raw.75M.$(READ1) raw.75M.$(READ2)
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
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
		MINLEN:25 2>> trim75 ; \
		cat 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq > 75M.left.$$TRIM.fq ; \
		cat 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq > 75M.right.$$TRIM.fq ; \
		rm 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq ; done
75M.2.Trinity.fasta 75M.5.Trinity.fasta 75M.10.Trinity.fasta 75M.20.Trinity.fasta raw.75M.Trinity.fasta: 75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 75M.left.$$TRIM.fq --right 75M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 75M.$$TRIM; rm 75M.left.$$TRIM.fq 75M.right.$$TRIM.fq; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left raw.75M.$(READ1) --right raw.75M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.75M; rm raw.75M.$(READ1) raw.75M.$(READ2)

raw.100M.$(READ1) raw.100M.$(READ2) : $(READ1) $(READ2)
	python ~/error_correction/scripts/subsampler.py 100000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.100M.$(READ1)
	mv subsamp_2.fastq raw.100M.$(READ2)	
100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq: raw.100M.$(READ1) raw.100M.$(READ2)
	@echo About to start trimming
	for TRIM in 2 5 10 20; do \
		java -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
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
		MINLEN:25 2>> trim100 ; \
		cat 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq > 100M.left.$$TRIM.fq ; \
		cat 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq > 100M.right.$$TRIM.fq ; \
		rm 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq ; done
100M.2.Trinity.fasta 100M.5.Trinity.fasta 100M.10.Trinity.fasta 100M.20.Trinity.fasta raw.100M.Trinity.fasta: 100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq
	for TRIM in 2 5 10 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
		--left 100M.left.$$TRIM.fq --right 100M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 100M.$$TRIM; rm 100M.left.$$TRIM.fq 100M.right.$$TRIM.fq; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 2 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G  \
	--left raw.100M.$(READ1) --right raw.100M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.100M; rm raw.100M.$(READ1) raw.100M.$(READ2)

pslx: 
#	for i in `ls *Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done
	for i in `ls raw.50M.Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done
	for i in `ls 75M.20.Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done
	for i in `ls raw.75M.Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done
	for i in `ls 100M.20.Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done
	for i in `ls raw.100M.Trinity.fasta`; do $(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query $$i; rm *maps *selected *summary; done

10M.2.Trinity.fasta.pep 10M.5.Trinity.fasta.pep 10M.10.Trinity.fasta.pep 10M.20.Trinity.fasta.pep 20M.2.Trinity.fasta.pep 20M.5.Trinity.fasta.pep 20M.10.Trinity.fasta.pep 20M.20.Trinity.fasta.pep 50M.2.Trinity.fasta.pep 50M.5.Trinity.fasta.pep 10M.50.Trinity.fasta.pep 50M.20.Trinity.fasta.pep 75M.2.Trinity.fasta.pep 75M.5.Trinity.fasta.pep 75M.10.Trinity.fasta.pep 75M.20.Trinity.fasta.pep 100M.2.Trinity.fasta.pep 100M.5.Trinity.fasta.pep 100M.10.Trinity.fasta.pep 100M.20.Trinity.fasta.pep: 
#	for i in `ls *Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done
	for i in `ls raw.50M.Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done
	for i in `ls 75M.20.Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done
	for i in `ls raw.75M.Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done
	for i in `ls 100M.20.Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done
	for i in `ls raw.100M.Trinity.fasta`; do $(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t $$i --search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $$i.pep; done














