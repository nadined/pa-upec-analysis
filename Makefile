
RAW_DATA_DIR = /home/nadine/PA_UPEC/reads
TRIMMED_DATA_DIR = $(RAW_DATA_DIR)/trimmed.thu-sep-27-12-35-14-2012
OUTPUT_DIR  = /home/nadine/BINF7000

VELVET_OPTIMISER = /home/nadine/VelvetOptimiser-2.2.5/VelvetOptimiser.pl

# older samtools version for use with srst2
#
# when calling this version use $(SAMTOOLS118)
#
# example:
# $(OUTPUT_DIR)/%_Assembly/align.bam: \
# 		$(OUTPUT_DIR)/%_Assembly/align.sam
# 	$(SAMTOOLS118) view -bS $< > $@
SAMTOOLS118 = /home/nadine/local/samtools-0.1.18/samtools

SRST2 = /home/nadine/srst2/scripts/srst2.py 

GETMLST = /home/nadine/srst2/scripts/getmlst.py

HASH_LENGTH = 49

START_HASH_LENGTH = 35
END_HASH_LENGTH = 71

CASE_IDS_FILE = cases.dat
PREFIXES      = $(shell cat $(CASE_IDS_FILE))

SPECIFIED_ASSEMBLIES = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly_$(HASH_LENGTH))

OPTIMISED_ASSEMBLIES = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly)

SORTED_BAMS = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/align_sorted.bam)

BAM_QC = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/align_sorted_stats)

KRAKEN_CLASS = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/contigs.kraken)

KRAKEN_LABEL =  $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/contigs.labels)

ANNOTATE = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/Annotation)

MLST = $(foreach case, $(PREFIXES), \
	$(OUTPUT_DIR)/$(case)_Assembly/MLST)

# preprocess assemblies (separate reads)
$(OUTPUT_DIR)/%_Assembly_Raw_$(HASH_LENGTH): \
		$(RAW_DATA_DIR)/%_R1.fastq.gz \
		$(RAW_DATA_DIR)/%_R2.fastq.gz
	velveth $@ $(HASH_LENGTH) -shortPaired \
		-fastq.gz -separate $< $(word 2, $^)

# preprocess assemblies (trimmed interleaved reads)
$(OUTPUT_DIR)/%_Assembly_$(HASH_LENGTH): \
		$(TRIMMED_DATA_DIR)/%_trimmed_ill.fastq.gz 
	velveth $@ $(HASH_LENGTH) -shortPaired \
		-fastq.gz -interleaved $< $(word 2, $^)

# create assemblies
$(OUTPUT_DIR)/%_Assembly_$(HASH_LENGTH)/contigs.fa: \
		$(OUTPUT_DIR)/%_Assembly_$(HASH_LENGTH)
	velvetg $<

# preprocess assemblies (trimmed reads)
$(OUTPUT_DIR)/%_Assembly/contigs.fa: \
		$(TRIMMED_DATA_DIR)/%_trimmed_ill.fastq.gz
	perl $(VELVET_OPTIMISER) \
		-d $(dir $@) \
		-s $(START_HASH_LENGTH) -e $(END_HASH_LENGTH) \
		-f '-fastq.gz -shortPaired $<' 

# create indexes
$(OUTPUT_DIR)/%_Assembly/contigs.fa.amb: \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa
	bwa index $<

# create sai
$(OUTPUT_DIR)/%_Assembly/align.sai: \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa \
		$(TRIMMED_DATA_DIR)/%_trimmed_ill.fastq.gz \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa.amb
	bwa aln -t 8 $< $(word 2, $^) > $@

# create sam
$(OUTPUT_DIR)/%_Assembly/align.sam: \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa \
		$(OUTPUT_DIR)/%_Assembly/align.sai \
		$(TRIMMED_DATA_DIR)/%_trimmed_ill.fastq.gz
	bwa samse $< $(word 2, $^) $(word 3, $^) > $@

# create bam
$(OUTPUT_DIR)/%_Assembly/align.bam: \
		$(OUTPUT_DIR)/%_Assembly/align.sam
	samtools view -bS $< > $@

# create sorted bam
$(OUTPUT_DIR)/%_Assembly/align_sorted.bam: \
		$(OUTPUT_DIR)/%_Assembly/align.bam
	samtools sort -m 6G $< $(basename $@)

#run qualimap bamqc 
$(OUTPUT_DIR)/%_Assembly/align_sorted_stats: \
		$(OUTPUT_DIR)/%_Assembly/align_sorted.bam
	qualimap bamqc -bam $<


#classify and label sequences with kraken
$(OUTPUT_DIR)/%_Assembly/contigs.kraken: \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa
	kraken --preload --db minikraken_20141208 $< \
		&& kraken --db minikraken_20141208 $< > $@

$(OUTPUT_DIR)/%_Assembly/contigs.labels: \
		$(OUTPUT_DIR)/%_Assembly/contigs.kraken
	kraken-translate --db minikraken_20141208 $< > $@
#annotate with Prokka
$(OUTPUT_DIR)/%_Assembly/Annotation: \
		$(OUTPUT_DIR)/%_Assembly/contigs.fa
	prokka --outdir $@ $<
#mlst with srst2
$(OUTPUT_DIR)/%_Assembly/MLST:\
		$(RAW_DATA_DIR)/%_R1.fastq.gz \
                $(RAW_DATA_DIR)/%_R2.fastq.gz
	python $(GETMLST) --species "Escherichia coli#1" && python $(SRST2) --output $@ --input_pe $< $(word 2, $^) --forward _R1 --reverse _R2 --mlst_db Escherichia_coli#1.fasta --mlst_definitions ecoli.txt --mlst_delimiter '-'


#find SNPs with kSNP3.0

#first use MakeFasta to combine genomes into single fasta
$(OUTPUT_DIR)/SNPs/combined_assemblies.fa:\
		$(OUTPUT_DIR)/SNPs/genomeslist
	MakeFasta $< $@ 

#use Kchooser to determine optimum k
$(OUTPUT_DIR)/SNPs/Kchooser.report:\
		$(OUTPUT_DIR)/SNPs/combined_assemblies.fa
	Kchooser $< $@

#run kSNP3
$(OUTPUT_DIR)/SNPs/testrun:\
		$(OUTPUT_DIR)/SNPs/genomeslist \
		$(OUTPUT_DIR)/SNPs/genomenames
	kSNP3 -in $< -k 19 -outdir $@ -annotate $(word 2, $^)

#find virulence genes with SeqFindr
$(OUTPUT_DIR)/virulence/VFDB_results:\
		$(OUTPUT_DIR)/virulence/Schembri_VFDB.fa \
		$(OUTPUT_DIR)/assemblies
	SeqFindr $<  $(word 2, $^)  -o $@ -l

$(OUTPUT_DIR)/virulence/plasmids_results:\
		$(OUTPUT_DIR)/virulence/plasmid_replicons.fa \
                $(OUTPUT_DIR)/assemblies	
	SeqFindr $<  $(word 2, $^) -o $@ -l

$(OUTPUT_DIR)/virulence/Islands_200_results: \
		$(OUTPUT_DIR)/virulence/Islands_200bp_chunks.fa \
		$(OUTPUT_DIR)/assemblies
	SeqFindr $< $(word 2, $^) -o $@ -l

$(OUTPUT_DIR)/virulence/Islands_500_results: \
                $(OUTPUT_DIR)/virulence/Islands_500bp_chunks.fa \
                $(OUTPUT_DIR)/assemblies
	SeqFindr $< $(word 2, $^) -o $@ -l

all: $(MLST)

test: $(OUTPUT_DIR)/SNPs/Kchooser.report

.INTERMEDIATE: \
	$(OUTPUT_DIR)/*_Assembly/contigs.fa.amb \
	$(OUTPUT_DIR)/*_Assembly/align.sam \
	$(OUTPUT_DIR)/*_Assembly/align.sai \
	$(OUTPUT_DIR)/*_Assembly/align.bam

clean:
	rm -rf $(SPECIFIED_ASSEMBLIES)
	rm -rf $(OPTIMISED_ASSEMBLIES)
	rm -f *Logfile.txti
