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



subsamp10:raw.10M.$(READ1) raw.10M.$(READ2)
trim10:10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq
trin10:10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta raw.10M.Trinity.fasta \
10M.2.Trinity.fasta.pslx 10M.5.Trinity.fasta.pslx 10M.10.Trinity.fasta.pslx 10M.20.Trinity.fasta.pslx raw.10M.Trinity.fasta.pslx \
10M.2.Trinity.fasta.pep 10M.5.Trinity.fasta.pep 10M.10.Trinity.fasta.pep 10M.20.Trinity.fasta.pep raw.10M.Trinity.fasta.pep \
10M.2.xprs 10M.5.xprs 10M.10.xprs 10M.20.xprs raw.10M.xprs

subsamp20:raw.20M.$(READ1) raw.20M.$(READ2)
trim20:20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq
trin20:20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta raw.20M.Trinity.fasta \
20M.2.Trinity.fasta.pslx 20M.5.Trinity.fasta.pslx 20M.10.Trinity.fasta.pslx 20M.20.Trinity.fasta.pslx raw.20M.Trinity.fasta.pslx \
20M.2.Trinity.fasta.pep 20M.5.Trinity.fasta.pep 20M.10.Trinity.fasta.pep 20M.20.Trinity.fasta.pep raw.20M.Trinity.fasta.pep \
20M.2.xprs 20M.5.xprs 20M.10.xprs 20M.20.xprs raw.20M.xprs

subsamp50:raw.50M.$(READ1) raw.50M.$(READ2)
trim50:50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq
trin50:50M.2.Trinity.fasta 50M.5.Trinity.fasta 50M.10.Trinity.fasta 50M.20.Trinity.fasta raw.50M.Trinity.fasta \
50M.2.Trinity.fasta.pslx 50M.5.Trinity.fasta.pslx 50M.10.Trinity.fasta.pslx 50M.20.Trinity.fasta.pslx raw.50M.Trinity.fasta.pslx \
50M.2.Trinity.fasta.pep 50M.5.Trinity.fasta.pep 50M.10.Trinity.fasta.pep 50M.20.Trinity.fasta.pep raw.50M.Trinity.fasta.pep \
50M.2.xprs 50M.5.xprs 50M.10.xprs 50M.20.xprs raw.50M.xprs

subsamp75:raw.75M.$(READ1) raw.75M.$(READ2)
trim75:75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq
trin75:75M.2.Trinity.fasta 75M.5.Trinity.fasta 75M.10.Trinity.fasta 75M.20.Trinity.fasta raw.75M.Trinity.fasta \
75M.2.Trinity.fasta.pslx 75M.5.Trinity.fasta.pslx 75M.10.Trinity.fasta.pslx 75M.20.Trinity.fasta.pslx raw.75M.Trinity.fasta.pslx \
75M.2.Trinity.fasta.pep 75M.5.Trinity.fasta.pep 75M.10.Trinity.fasta.pep 75M.20.Trinity.fasta.pep raw.75M.Trinity.fasta.pep \
75M.2.xprs 75M.5.xprs 75M.10.xprs 75M.20.xprs raw.75M.xprs

subsamp100:raw.100M.$(READ1) raw.100M.$(READ2)
trim100:100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq
trin100:100M.2.Trinity.fasta 100M.5.Trinity.fasta 100M.10.Trinity.fasta 100M.20.Trinity.fasta raw.100M.Trinity.fasta \
100M.2.Trinity.fasta.pslx 100M.5.Trinity.fasta.pslx 100M.10.Trinity.fasta.pslx 100M.20.Trinity.fasta.pslx raw.100M.Trinity.fasta.pslx \
100M.2.Trinity.fasta.pep 100M.5.Trinity.fasta.pep 100M.10.Trinity.fasta.pep 100M.20.Trinity.fasta.pep raw.100M.Trinity.fasta.pep \
100M.2.xprs 100M.5.xprs 100M.10.xprs 100M.20.xprs raw.100M.xprs


$(READ1).quality: $(READ1) $(READ2)
	perl $(SOLEXA)/SolexaQA.pl -p 0.01 $(READ1)
	cp $(READ1).quality ~/Dropbox/
	cp $(READ1).quality.pdf ~/Dropbox/


raw.10M.$(READ1) raw.10M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 10000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.10M.$(READ1)
	mv subsamp_2.fastq raw.10M.$(READ2)	
10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq: 
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
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 1>> trim10.log; \
		cat 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq > 10M.left.$$TRIM.fq ; \
		cat 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq > 10M.right.$$TRIM.fq ; \
		rm 10M.$$TRIM.pp.2.fq 10M.$$TRIM.up.2.fq 10M.$$TRIM.pp.1.fq 10M.$$TRIM.up.1.fq ; done
10M.2.Trinity.fasta 10M.5.Trinity.fasta 10M.10.Trinity.fasta 10M.20.Trinity.fasta raw.10M.Trinity.fasta \
10M.2.Trinity.fasta.pslx 10M.5.Trinity.fasta.pslx 10M.10.Trinity.fasta.pslx 10M.20.Trinity.fasta.pslx raw.10M.Trinity.fasta.pslx \
10M.2.Trinity.fasta.pep 10M.5.Trinity.fasta.pep 10M.10.Trinity.fasta.pep 10M.20.Trinity.fasta.pep raw.10M.Trinity.fasta.pep \
10M.2.xprs 10M.5.xprs 10M.10.xprs 10M.20.xprs raw.10M.xprs: 10M.left.2.fq 10M.left.5.fq 10M.left.10.fq 10M.left.20.fq 10M.right.2.fq 10M.right.5.fq 10M.right.10.fq 10M.right.20.fq
	for TRIM in 20 2 5 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
		--left 10M.left.$$TRIM.fq --right 10M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 10M.$$TRIM; \
##FL Reconstruction
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 10M.$$TRIM.Trinity.fasta; rm *maps *selected *summary; \
##ORF ID
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 10M.$$TRIM.Trinity.fasta \
		--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; \
		rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 10M.$$TRIM.Trinity.fasta.pep; \
##Mapping and eXpress
		bowtie2-build -q 10M.$$TRIM.Trinity.fasta index; \
		bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o 10.$$TRIM.xprs -p8 10M.$$TRIM.Trinity.fasta >>10M.$$TRIM.mapping.log ; rm index* ; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left raw.10M.$(READ1) --right raw.10M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.10M
##FL Reconstruction
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.10M.Trinity.fasta; rm *maps *selected *summary
##ORF ID
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.10M.Trinity.fasta \
	--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep raw.10M.Trinity.fasta.pep
##mapping and eXpress
	bowtie2-build -q raw.10M.Trinity.fasta index;bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o raw.10.TRIM.xprs -p8 raw.10M.Trinity.fasta >>raw.10M.mapping.log; rm index*



raw.20M.$(READ1) raw.20M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 20000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.20M.$(READ1)
	mv subsamp_2.fastq raw.20M.$(READ2)	
20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq: 
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
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 1>> trim10.log; \
		cat 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq > 20M.left.$$TRIM.fq ; \
		cat 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq > 20M.right.$$TRIM.fq ; \
		rm 20M.$$TRIM.pp.2.fq 20M.$$TRIM.up.2.fq 20M.$$TRIM.pp.1.fq 20M.$$TRIM.up.1.fq ; done
20M.2.Trinity.fasta 20M.5.Trinity.fasta 20M.10.Trinity.fasta 20M.20.Trinity.fasta raw.20M.Trinity.fasta \
20M.2.Trinity.fasta.pslx 20M.5.Trinity.fasta.pslx 20M.10.Trinity.fasta.pslx 20M.20.Trinity.fasta.pslx raw.20M.Trinity.fasta.pslx \
20M.2.Trinity.fasta.pep 20M.5.Trinity.fasta.pep 20M.10.Trinity.fasta.pep 20M.20.Trinity.fasta.pep raw.20M.Trinity.fasta.pep \
20M.2.xprs 20M.5.xprs 20M.10.xprs 20M.20.xprs raw.20M.xprs: 20M.left.2.fq 20M.left.5.fq 20M.left.10.fq 20M.left.20.fq 20M.right.2.fq 20M.right.5.fq 20M.right.10.fq 20M.right.20.fq
	for TRIM in 20 2 5 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
		--left 20M.left.$$TRIM.fq --right 20M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 20M.$$TRIM; \
##FL Reconstruction
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 20M.$$TRIM.Trinity.fasta; rm *maps *selected *summary; \
##ORF ID
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 20M.$$TRIM.Trinity.fasta \
		--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; \
		rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 20M.$$TRIM.Trinity.fasta.pep; \
##Mapping and eXpress
		bowtie2-build -q 20M.$$TRIM.Trinity.fasta index; \
		bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o 10.$$TRIM.xprs -p8 20M.$$TRIM.Trinity.fasta >>20M.$$TRIM.mapping.log ; rm index* ; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left raw.20M.$(READ1) --right raw.20M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.20M
##FL Reconstruction
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.20M.Trinity.fasta; rm *maps *selected *summary
##ORF ID
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.20M.Trinity.fasta \
	--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep raw.20M.Trinity.fasta.pep
##mapping and eXpress
	bowtie2-build -q raw.20M.Trinity.fasta index;bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o raw.10.TRIM.xprs -p8 raw.20M.Trinity.fasta >>raw.20M.mapping.log; rm index*





raw.50M.$(READ1) raw.50M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 50000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.50M.$(READ1)
	mv subsamp_2.fastq raw.50M.$(READ2)	
50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq: 
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
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 1>> trim10.log; \
		cat 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq > 50M.left.$$TRIM.fq ; \
		cat 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq > 50M.right.$$TRIM.fq ; \
		rm 50M.$$TRIM.pp.2.fq 50M.$$TRIM.up.2.fq 50M.$$TRIM.pp.1.fq 50M.$$TRIM.up.1.fq ; done
50M.2.Trinity.fasta 50M.5.Trinity.fasta 50M.10.Trinity.fasta 50M.20.Trinity.fasta raw.50M.Trinity.fasta \
50M.2.Trinity.fasta.pslx 50M.5.Trinity.fasta.pslx 50M.10.Trinity.fasta.pslx 50M.20.Trinity.fasta.pslx raw.50M.Trinity.fasta.pslx \
50M.2.Trinity.fasta.pep 50M.5.Trinity.fasta.pep 50M.10.Trinity.fasta.pep 50M.20.Trinity.fasta.pep raw.50M.Trinity.fasta.pep \
50M.2.xprs 50M.5.xprs 50M.10.xprs 50M.20.xprs raw.50M.xprs: 50M.left.2.fq 50M.left.5.fq 50M.left.10.fq 50M.left.20.fq 50M.right.2.fq 50M.right.5.fq 50M.right.10.fq 50M.right.20.fq
	for TRIM in 20 2 5 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
		--left 50M.left.$$TRIM.fq --right 50M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 50M.$$TRIM; \
##FL Reconstruction
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 50M.$$TRIM.Trinity.fasta; rm *maps *selected *summary; \
##ORF ID
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 50M.$$TRIM.Trinity.fasta \
		--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; \
		rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 50M.$$TRIM.Trinity.fasta.pep; \
##Mapping and eXpress
		bowtie2-build -q 50M.$$TRIM.Trinity.fasta index; \
		bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o 10.$$TRIM.xprs -p8 50M.$$TRIM.Trinity.fasta >>50M.$$TRIM.mapping.log ; rm index* ; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left raw.50M.$(READ1) --right raw.50M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.50M
##FL Reconstruction
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.50M.Trinity.fasta; rm *maps *selected *summary
##ORF ID
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.50M.Trinity.fasta \
	--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep raw.50M.Trinity.fasta.pep
##mapping and eXpress
	bowtie2-build -q raw.50M.Trinity.fasta index;bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o raw.10.TRIM.xprs -p8 raw.50M.Trinity.fasta >>raw.50M.mapping.log; rm index*


raw.75M.$(READ1) raw.75M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 75000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.75M.$(READ1)
	mv subsamp_2.fastq raw.75M.$(READ2)	
75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq: 
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
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 1>> trim10.log; \
		cat 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq > 75M.left.$$TRIM.fq ; \
		cat 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq > 75M.right.$$TRIM.fq ; \
		rm 75M.$$TRIM.pp.2.fq 75M.$$TRIM.up.2.fq 75M.$$TRIM.pp.1.fq 75M.$$TRIM.up.1.fq ; done
75M.2.Trinity.fasta 75M.5.Trinity.fasta 75M.10.Trinity.fasta 75M.20.Trinity.fasta raw.75M.Trinity.fasta \
75M.2.Trinity.fasta.pslx 75M.5.Trinity.fasta.pslx 75M.10.Trinity.fasta.pslx 75M.20.Trinity.fasta.pslx raw.75M.Trinity.fasta.pslx \
75M.2.Trinity.fasta.pep 75M.5.Trinity.fasta.pep 75M.10.Trinity.fasta.pep 75M.20.Trinity.fasta.pep raw.75M.Trinity.fasta.pep \
75M.2.xprs 75M.5.xprs 75M.10.xprs 75M.20.xprs raw.75M.xprs: 75M.left.2.fq 75M.left.5.fq 75M.left.10.fq 75M.left.20.fq 75M.right.2.fq 75M.right.5.fq 75M.right.10.fq 75M.right.20.fq
	for TRIM in 20 2 5 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
		--left 75M.left.$$TRIM.fq --right 75M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 75M.$$TRIM; \
##FL Reconstruction
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 75M.$$TRIM.Trinity.fasta; rm *maps *selected *summary; \
##ORF ID
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 75M.$$TRIM.Trinity.fasta \
		--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; \
		rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 75M.$$TRIM.Trinity.fasta.pep; \
##Mapping and eXpress
		bowtie2-build -q 75M.$$TRIM.Trinity.fasta index; \
		bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o 10.$$TRIM.xprs -p8 75M.$$TRIM.Trinity.fasta >>75M.$$TRIM.mapping.log ; rm index* ; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left raw.75M.$(READ1) --right raw.75M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.75M
##FL Reconstruction
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.75M.Trinity.fasta; rm *maps *selected *summary
##ORF ID
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.75M.Trinity.fasta \
	--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep raw.75M.Trinity.fasta.pep
##mapping and eXpress
	bowtie2-build -q raw.75M.Trinity.fasta index;bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o raw.10.TRIM.xprs -p8 raw.75M.Trinity.fasta >>raw.75M.mapping.log; rm index*


raw.100M.$(READ1) raw.100M.$(READ2): 
	python ~/error_correction/scripts/subsampler.py 100000000 $(READ1) $(READ2)
	mv subsamp_1.fastq raw.100M.$(READ1)
	mv subsamp_2.fastq raw.100M.$(READ2)	
100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq: 
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
		ILLUMINACLIP:$(BCODES):2:40:15 \
		LEADING:$$TRIM \
		TRAILING:$$TRIM \
		SLIDINGWINDOW:4:$$TRIM \
		MINLEN:25 1>> trim10.log; \
		cat 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq > 100M.left.$$TRIM.fq ; \
		cat 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq > 100M.right.$$TRIM.fq ; \
		rm 100M.$$TRIM.pp.2.fq 100M.$$TRIM.up.2.fq 100M.$$TRIM.pp.1.fq 100M.$$TRIM.up.1.fq ; done
100M.2.Trinity.fasta 100M.5.Trinity.fasta 100M.10.Trinity.fasta 100M.20.Trinity.fasta raw.100M.Trinity.fasta \
100M.2.Trinity.fasta.pslx 100M.5.Trinity.fasta.pslx 100M.10.Trinity.fasta.pslx 100M.20.Trinity.fasta.pslx raw.100M.Trinity.fasta.pslx \
100M.2.Trinity.fasta.pep 100M.5.Trinity.fasta.pep 100M.10.Trinity.fasta.pep 100M.20.Trinity.fasta.pep raw.100M.Trinity.fasta.pep \
100M.2.xprs 100M.5.xprs 100M.10.xprs 100M.20.xprs raw.100M.xprs: 100M.left.2.fq 100M.left.5.fq 100M.left.10.fq 100M.left.20.fq 100M.right.2.fq 100M.right.5.fq 100M.right.10.fq 100M.right.20.fq
	for TRIM in 20 2 5 10; do \
		$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
		--left 100M.left.$$TRIM.fq --right 100M.right.$$TRIM.fq --group_pairs_distance 999 --CPU $(CPU) --output 100M.$$TRIM; \
##FL Reconstruction
		$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query 100M.$$TRIM.Trinity.fasta; rm *maps *selected *summary; \
##ORF ID
		$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t 100M.$$TRIM.Trinity.fasta \
		--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm; \
		rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep 100M.$$TRIM.Trinity.fasta.pep; \
##Mapping and eXpress
		bowtie2-build -q 100M.$$TRIM.Trinity.fasta index; \
		bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o 10.$$TRIM.xprs -p8 100M.$$TRIM.Trinity.fasta >>100M.$$TRIM.mapping.log ; rm index* ; done
	$(TRINITY)/Trinity.pl --full_cleanup --min_kmer_cov 1 --seqType fq --JM $(MEM)G --bflyHeapSpaceMax $(MEM)G \
	--left raw.100M.$(READ1) --right raw.100M.$(READ2) --group_pairs_distance 999 --CPU $(CPU) --output raw.100M
##FL Reconstruction
	$(TRINITY)/Analysis/FL_reconstruction_analysis/FL_trans_analysis_pipeline.pl --target $(MUS) --query raw.100M.Trinity.fasta; rm *maps *selected *summary
##ORF ID
	$(TRINITY)/trinity-plugins/transdecoder/transcripts_to_best_scoring_ORFs.pl --CPU $(CPU) -t raw.100M.Trinity.fasta \
	--search_pfam /media/macmanes/hd2/pfam/Pfam-A.hmm
	rm longest_orfs* *gff3 *dat *scores *cds *bed *inx; mv best_candidates.eclipsed_orfs_removed.pep raw.100M.Trinity.fasta.pep
##mapping and eXpress
	bowtie2-build -q raw.100M.Trinity.fasta index;bowtie2 -p 12 -X 999 -k 30 -x index -1 $(READ1) -2 $(READ2) | express -o raw.10.TRIM.xprs -p8 raw.100M.Trinity.fasta >>raw.100M.mapping.log; rm index*
