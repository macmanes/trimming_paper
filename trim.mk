#!/usr/bin/make -rRsf


SOLEXA := /home/macmanes/SolexaQA_v.2.1/SolexaQA_v.2.1
TRINITY := /home/macmanes/trinityrnaseq_r2013-02-25
MUS := /media/macmanes/hd/flux/genomes/mus/Mus_musculus.GRCm38.71.cdna.all.fa
PFAM := /media/macmanes/raid/blastdb/Pfam-A.hmm
CPU=2
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')

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
real: out_1.fastq.quality out_1.fastq raw.Trinity.fasta right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq real.1.Trinity.fasta \
	real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta
sim: right.fq.quality sim.Trinity.fasta sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq sim.1.Trinity.fasta \
	sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta sim.20.Trinity.fasta 
<<<<<<< HEAD
pslx: sim.1.Trinity.fasta.pslx sim.2.Trinity.fasta.pslx sim.5.Trinity.fasta.pslx sim.10.Trinity.fasta.pslx sim.15.Trinity.fasta.pslx sim.20.Trinity.fasta.pslx \
	sim.Trinity.fasta.pslx real.1.Trinity.fasta.pslx real.2.Trinity.fasta.pslx real.5.Trinity.fasta.pslx real.10.Trinity.fasta.pslx real.15.Trinity.fasta.pslx \
	real.20.Trinity.fasta.pslx raw.Trinity.fasta.pslx
full:
=======

pslx: 
>>>>>>> 799b5d39f365af57cc69a6c5acc6ba1f78a79f48

out_1.fastq.quality:out_1.fastq
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 out_1.fastq
	cp out_1.fastq.quality ~/Dropbox/
	cp out_1.fastq.quality.pdf ~/Dropbox/
	
raw.Trinity.fasta:out_1.fastq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left out_1.fastq  --right out_2.fastq  --CPU $(CPU) --output raw

right.1.fq right.2.fq right.5.fq right.10.fq right.15.fq right.20.fq: out_1.fastq
	@echo About to start trimming
	for TRIM in 1 2 5 10 15 20; do \
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
	
<<<<<<< HEAD
sim.Trinity.fasta:right.fq.fastq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.fq  --right right.fq  --CPU $(CPU) --output sim
=======
sim.Trinity.fasta:right.fq
	$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left left.fq  --right right.fq  --CPU 8 --output sim
>>>>>>> 799b5d39f365af57cc69a6c5acc6ba1f78a79f48

sim.left.1.fq sim.left.2.fq sim.left.5.fq sim.left.10.fq sim.left.15.fq sim.left.20.fq: right.fq
	for TRIM in 1 2 5 10 15 20; do \
		java -Xmx30g -jar $(TRIMMOMATIC) PE \
		-phred33 -threads $(CPU) \
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
		$(TRINITY)/Trinity.pl --full_cleanup --seqType fq --JM 30G --left sim.left.$$TRIM.fq --right sim.right.$$TRIM.fq --CPU $(CPU) --output sim.$$TRIM ; \
	done



<<<<<<< HEAD
sim.1.Trinity.fasta.pslx sim.2.Trinity.fasta.pslx sim.5.Trinity.fasta.pslx sim.10.Trinity.fasta.pslx sim.15.Trinity.fasta.pslx sim.20.Trinity.fasta.pslx \
sim.Trinity.fasta.pslx real.1.Trinity.fasta.pslx real.2.Trinity.fasta.pslx real.5.Trinity.fasta.pslx real.10.Trinity.fasta.pslx real.15.Trinity.fasta.pslx \
real.20.Trinity.fasta.pslx raw.Trinity.fasta.pslx: sim.1.Trinity.fasta sim.2.Trinity.fasta sim.5.Trinity.fasta sim.10.Trinity.fasta sim.15.Trinity.fasta \
sim.20.Trinity.fasta sim.Trinity.fasta real.1.Trinity.fasta real.2.Trinity.fasta real.5.Trinity.fasta real.10.Trinity.fasta real.15.Trinity.fasta real.20.Trinity.fasta raw.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.2.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.1.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.5.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.10.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.15.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.20.Trinity.fasta
=======
pslx: 
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.2.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.1.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.5.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.10.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.15.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query sim.20.Trinity.fasta
>>>>>>> 799b5d39f365af57cc69a6c5acc6ba1f78a79f48
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.1.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.2.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.5.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.10.Trinity.fasta
#	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.15.Trinity.fasta
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query real.20.Trinity.fasta
	rm *maps *selected *summary


$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.Trinity.fasta --search_pfam $(PFAM)





































