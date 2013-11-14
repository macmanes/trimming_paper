#!/usr/bin/make -rRsf

#USAGE:

MEM=5
TRIM=5
CPU=5
RUN=run
READ1=left.fastq
READ2=right.fastq
TRIMMOMATIC ?= $(shell which 'trimmomatic-0.30.jar')
BCODES := /usr/users/1/macmanes/software/barcodes.fa
TRINITY := /usr/local/packages/trinity/r2013-08-14/
MUS := /brashear/macmanes/runs/trim/Mus_musculus.GRCm38.71.cdna.all.fa
PFAM := /brashear/macmanes/runs/trim/pfam/Pfam-A.hmm

<<<<<<< HEAD
raw.50M.$(READ1) raw.50M.$(READ2):
	python ~/error_correction/scripts/subsampler.py 50000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.50M.$(READ1)
	mv subsamp_2.fastq raw.50M.$(READ2)
=======
raw.10M.$(READ1) raw.10M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 10000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)	
trim10:
	@echo About to start trimming
	java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
	-phred33 -threads $(CPU) \
	../raw.10M.$(READ1) \
	../raw.10M.$(READ2) \
	10M.$$TRIM.pp.1.fq \
	10M.$$TRIM.up.1.fq \
	10M.$$TRIM.pp.2.fq \
	10M.$$TRIM.up.2.fq \
	ILLUMINACLIP:$(BCODES):2:40:15 \
	LEADING:$$TRIM \
	TRAILING:$$TRIM \
	SLIDINGWINDOW:4:$$TRIM \
	MINLEN:25 2>> trim10.log; \
	cat 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq > 10M.left.$$TRIM.fq ; \
	cat 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq > 10M.right.$$TRIM.fq ; \
	rm 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq
trin10:
	$(TRINITY)/Trinity.pl --bflyGCThread 25 --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left 10M.left.$$TRIM.fq --right 10M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 10M.$$TRIM >> 10trin$$TRIM.out
pslx10:
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 10M.$$TRIM.Trinity.fasta; rm *maps *selected *summary
pep10:
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 10M.$$TRIM.Trinity.fasta \
	--search_pfam $(PFAM) >>pfam10.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 10M.$$TRIM.Trinity.fasta.pep
map10:
	bowtie2-build -q 10M.$$TRIM.Trinity.fasta index; echo -e '\n' Mapping at PHRED=$$TRIM '\n' >> 10M.$$TRIM.mapping.log; \
	bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) 2>>10M.$$TRIM.mapping.log | express -o 10M.$$TRIM.xprs -p8 10M.$$TRIM.Trinity.fasta 2>>10M.$$TRIM.mapping.log; rm index*

raw.20M.$(READ1) raw.20M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 20000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.20M.$(READ1)
	mv subsamp_2.fastq raw.20M.$(READ2)	
trim20:
	@echo About to start trimming
	java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
	-phred33 -threads $(CPU) \
	../raw.20M.$(READ1) \
	../raw.20M.$(READ2) \
	20M.$$TRIM.pp.1.fq \
	20M.$$TRIM.up.1.fq \
	20M.$$TRIM.pp.2.fq \
	20M.$$TRIM.up.2.fq \
	ILLUMINACLIP:$(BCODES):2:40:15 \
	LEADING:$$TRIM \
	TRAILING:$$TRIM \
	SLIDINGWINDOW:4:$$TRIM \
	MINLEN:25 2>> trim20.log; \
	cat 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq > 20M.left.$$TRIM.fq ; \
	cat 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq > 20M.right.$$TRIM.fq ; \
	rm 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq
trin20:
	$(TRINITY)/Trinity.pl --bflyGCThread 25 --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left 20M.left.$$TRIM.fq --right 20M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 20M.$$TRIM >> 20trin$$TRIM.out
pslx20:
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 20M.$$TRIM.Trinity.fasta; rm *maps *selected *summary
pep20:
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 20M.$$TRIM.Trinity.fasta \
	--search_pfam $(PFAM) >>pfam20.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 20M.$$TRIM.Trinity.fasta.pep
map20:
	bowtie2-build -q 20M.$$TRIM.Trinity.fasta index; echo -e '\n' Mapping at PHRED=$$TRIM '\n' >> 20M.$$TRIM.mapping.log; \
	bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) 2>>20M.$$TRIM.mapping.log | express -o 20M.$$TRIM.xprs -p8 20M.$$TRIM.Trinity.fasta 2>>20M.$$TRIM.mapping.log; rm index*



raw.50M.$(READ1) raw.50M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 50000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.50M.$(READ1)
	mv subsamp_2.fastq raw.50M.$(READ2)	
>>>>>>> a147beb60c68e837a1c79fb6213480f746e753fa
trim50:
	@echo About to start trimming
	java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
	-phred33 -threads $(CPU) \
	../raw.50M.$(READ1) \
	../raw.50M.$(READ2) \
	50M.$$TRIM.pp.1.fq \
	50M.$$TRIM.up.1.fq \
	50M.$$TRIM.pp.2.fq \
	50M.$$TRIM.up.2.fq \
	ILLUMINACLIP:$(BCODES):2:40:15 \
	LEADING:$$TRIM \
	TRAILING:$$TRIM \
	SLIDINGWINDOW:4:$$TRIM \
	MINLEN:25 2>> trim50.log; \
	cat 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq > 50M.left.$$TRIM.fq ; \
	cat 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq > 50M.right.$$TRIM.fq ; \
	rm 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq
trin50:
	$(TRINITY)/Trinity.pl --bflyGCThread 25 --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left 50M.left.$$TRIM.fq --right 50M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 50M.$$TRIM >> 50trin$$TRIM.out
pslx50:
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 50M.$$TRIM.Trinity.fasta; rm *maps *selected *summary
pep50:
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 50M.$$TRIM.Trinity.fasta \
	--search_pfam $(PFAM) >>pfam50.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 50M.$$TRIM.Trinity.fasta.pep
map50:
	bowtie2-build -q 50M.$$TRIM.Trinity.fasta index; echo -e '\n' Mapping at PHRED=$$TRIM '\n' >> 50M.$$TRIM.mapping.log; \
	bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) 2>>50M.$$TRIM.mapping.log | express -o 50.$$TRIM.xprs -p8 50M.$$TRIM.Trinity.fasta 2>>50M.$$TRIM.mapping.log; rm index*

<<<<<<< HEAD
raw.75M.$(READ1) raw.75M.$(READ2):
	python ~/error_correction/scripts/subsampler.py 75000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.75M.$(READ1)
	mv subsamp_2.fastq raw.75M.$(READ2)
=======
raw.75M.$(READ1) raw.75M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 75000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.75M.$(READ1)
	mv subsamp_2.fastq raw.75M.$(READ2)	
>>>>>>> a147beb60c68e837a1c79fb6213480f746e753fa
trim75:
	@echo About to start trimming
	mkdir $$TRIM.trim75
	cd $$TRIM.trim75
	java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
	-phred33 -threads $(CPU) \
	../raw.75M.$(READ1) \
	../raw.75M.$(READ2) \
	75M.$$TRIM.pp.1.fq \
	75M.$$TRIM.up.1.fq \
	75M.$$TRIM.pp.2.fq \
	75M.$$TRIM.up.2.fq \
	ILLUMINACLIP:$(BCODES):2:40:15 \
	LEADING:$$TRIM \
	TRAILING:$$TRIM \
	SLIDINGWINDOW:4:$$TRIM \
	MINLEN:25 2>> trim75.log; \
	cat 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq > 75M.left.$$TRIM.fq ; \
	cat 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq > 75M.right.$$TRIM.fq ; \
	rm 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq
trin75:
	$(TRINITY)/Trinity.pl --bflyGCThread 25 --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left 75M.left.$$TRIM.fq --right 75M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 75M.$$TRIM >> 75trin$$TRIM.out
pslx75:
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 75M.$$TRIM.Trinity.fasta; rm *maps *selected *summary
pep75:
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 75M.$$TRIM.Trinity.fasta \
	--search_pfam $(PFAM) >> pfam75.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 75M.$$TRIM.Trinity.fasta.pep
map75:
	bowtie2-build -q 75M.$$TRIM.Trinity.fasta index; echo -e '\n' Mapping at PHRED=$$TRIM '\n' >> 75M.$$TRIM.mapping.log; \
<<<<<<< HEAD
	bowtie2 -p 12 -X 999 -k 30 -x index -1 ../$(READ1) -2 ../$(READ2) 2>>75M.$$TRIM.mapping.log | express -o 75.$$TRIM.xprs -p8 75M.$$TRIM.Trinity.fasta 2>>75M.$$TRIM.mapping.log ; rm index*

raw.100M.$(READ1) raw.100M.$(READ2):
	python ~/error_correction/scripts/subsampler.py 100000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.100M.$(READ1)
	mv subsamp_2.fastq raw.100M.$(READ2)
=======
	bowtie2 -p 12 -X 999 -k 30 -x index -1 ../$(READ1) -2 ../$(READ2) 2>>75M.$$TRIM.mapping.log | express -o 75.$$TRIM.xprs -p8 75M.$$TRIM.Trinity.fasta 2>>75M.$$TRIM.mapping.log ; rm index* 


raw.100M.$(READ1) raw.100M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 100000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.100M.$(READ1)
	mv subsamp_2.fastq raw.100M.$(READ2)	
>>>>>>> a147beb60c68e837a1c79fb6213480f746e753fa
trim100:
	@echo About to start trimming
	java -XX:ParallelGCThreads=32 -Xmx$(MEM)g -jar $(TRIMMOMATIC) PE \
	-phred33 -threads $(CPU) \
	../raw.100M.$(READ1) \
	../raw.100M.$(READ2) \
	100M.$$TRIM.pp.1.fq \
	100M.$$TRIM.up.1.fq \
	100M.$$TRIM.pp.2.fq \
	100M.$$TRIM.up.2.fq \
	ILLUMINACLIP:$(BCODES):2:40:15 \
	LEADING:$$TRIM \
	TRAILING:$$TRIM \
	SLIDINGWINDOW:4:$$TRIM \
	MINLEN:25 2>> trim100.log; \
	cat 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq > 100M.left.$$TRIM.fq ; \
	cat 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq > 100M.right.$$TRIM.fq ; \
<<<<<<< HEAD
	rm 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq
=======
	rm 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq 
>>>>>>> a147beb60c68e837a1c79fb6213480f746e753fa
trin100:
	$(TRINITY)/Trinity.pl --bflyGCThread 25 --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left 100M.left.$$TRIM.fq --right 100M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 100M.$$TRIM >> 100trin$$TRIM.out
pslx100:
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 100M.$$TRIM.Trinity.fasta; rm *maps *selected *summary
pep100:
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 100M.$$TRIM.Trinity.fasta \
	--search_pfam $(PFAM) >> pfam100.log; \
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 100M.$$TRIM.Trinity.fasta.pep
map100:
	bowtie2-build -q 100M.$$TRIM.Trinity.fasta index; echo -e '\n' Mapping at PHRED=$$TRIM '\n' >> 100M.$$TRIM.mapping.log
<<<<<<< HEAD
	bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) 2>>100M.$$TRIM.mapping.log | express -o 100.$$TRIM.xprs -p8 100M.$$TRIM.Trinity.fasta 2>>100M.$$TRIM.mapping.log ; rm index*
=======
	bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) 2>>100M.$$TRIM.mapping.log | express -o 100.$$TRIM.xprs -p8 100M.$$TRIM.Trinity.fasta 2>>100M.$$TRIM.mapping.log ; rm index*
>>>>>>> a147beb60c68e837a1c79fb6213480f746e753fa
