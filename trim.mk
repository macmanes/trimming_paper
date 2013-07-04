#!/usr/bin/make -rRsf

SOLEXA ?= $(shell which 'SolexaQA.pl')
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
TRIMMOMATIC := /home/macmanes/software
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa

real: out_1.fastq.quality out_1.fastq raw.Trinity.fasta right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq real.1.Trinity.fasta \
	real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta

sim: right.fq.quality sim.Trinity.fasta sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq sim.1.Trinity.fasta \
	sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta sim.20.Trinity.fasta 

pslx: 

out_1.fastq.quality:out_1.fastq
	perl $(SOLEXA) -p 0.01 out_1.fastq
	cp out_1.fastq.quality ~/Dropbox/
	cp out_1.fastq.quality.pdf ~/Dropbox/
	
raw.Trinity.fasta:out_1.fastq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left out_1.fastq  --right out_2.fastq  --CPU 8 --output raw

right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq: out_1.fastq
	@echo About to start trimming
	for TRIM in 1 2 5 10 15 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC)/trimmomatic-0.30.jar PE \
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
		perl $(SOLEXA)/SolexaQA.pl -p 0.01 T.$$TRIM.pp.1.fq ; \
		cp T.$$TRIM.pp.1.fq.quality ~/Dropbox/ ; \
		cp T.$$TRIM.pp.1.fq.quality.pdf ~/Dropbox/ ;\
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > right.$$TRIM.fq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta: out_1.fastq
	for TRIM in 1 2 5 10 15 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.$$TRIM.fq --right right.$$TRIM.fq --CPU 8 --output real.$$TRIM ; \
	done
	
##Make sim	
	
right.fq.quality:right.fq
	perl $(SOLEXA) -p 0.01 right.fq
	cp right.fq.quality ~/Dropbox/
	cp right.fq.quality.pdf ~/Dropbox/
	
sim.Trinity.fasta:right.fq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.fq  --right right.fq  --CPU 8 --output sim

sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq: right.fq
	for TRIM in 1 2 5 10 15 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC)/trimmomatic-0.30.jar PE \
		-phred33 -threads 12 \
		left.fq \
		right.fq \
		T.$$TRIM.pp.1.fq \
		T.$$TRIM.up.1.fq \
		T.$$TRIM.pp.2.fq \
		T.$$TRIM.up.2.fq \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 ; \
		perl $(SOLEXA) -p 0.01 T.$$TRIM.pp.1.fq ; \
		cp T.$$TRIM.pp.1.fq.quality ~/Dropbox/ ; \
		cp T.$$TRIM.pp.1.fq.quality.pdf ~/Dropbox/ ;\
		rm rm *matrix *segments ;\
		cat T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq > sim.left.$$TRIM.fq ; \
		cat T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq > sim.right.$$TRIM.fq ; \
		rm T.$$TRIM.pp.2.fq T.$$TRIM.up.2.fq T.$$TRIM.pp.1.fq T.$$TRIM.up.1.fq ; \
	done

sim.1.Trinity.fasta sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta sim.20.Trinity.fasta: right.fq
	for TRIM in 1 2 5 10 15 20; do \
		$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left sim.left.$$TRIM.fq --right sim.right.$$TRIM.fq --CPU 8 --output sim.$$TRIM ; \
	done



pslx: 
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.2.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.1.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.5.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.10.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.15.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.20.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.1.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.2.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.5.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.10.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.15.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.20.Trinity.fasta
	rm *maps *selected *summary

orf: 
	for i in `ls *Trinity.fasta`; do \
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU 4 -t $$i --search_pfam /media/macmanes/raid/blastdb/Pfam-A.hmm; rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep $i.pep; done



















