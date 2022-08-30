"""
Author: J. Hughes and Salvo Camiolo
Affiliation: CVR Bioinformatics
Aim: A Snakemake workflow to detect low frequency mutations in HCMV genome sequenced using Illumina.
Date: 26 Aug 2022
Run: snakemake -s HCMVmut.smk
Latest modification: 
  - todo
  
  View the DAG: snakemake --forceall --rulegraph -s HCMVmut.smk | dot -Tpng > dag.png
"""

configfile: "config.yml"

# Pipeline for quality filtering and low frequency mutation detection from Illumina reads
# 300 base paired-end

# you can override PREFIX and OUTPUT by running: snakemake -C prefix=PREFIX output=OUTPUT.


# pull in all files with .fastq on the end in the 'data' directory.
PREFIX = config["fastq_directory"]
print(f"Looking for FASTQ files under '{PREFIX}'/")

OUTPUT = config["output"]
print(f"Results will go under '{OUTPUT}'/")

FILES = glob_wildcards(f'{PREFIX}/{{name}}_R1_001.fastq')

# extract the {name} values into a list
NAMES = FILES.name

# for creating the different folders with the different trimming settings
FILTER = ["tr","tr_dd","tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN", "markedDuplicates"]

# request the output files
rule all:
    input:
        #expand("{output}/RawReads/{name}_R1_001.fastq", output=OUTPUT, name=NAMES),
        #expand("{output}/tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN/{name}_1.fastq", output=OUTPUT, name=NAMES),
        #expand("{output}/tr_dd/{name}_1.fastq", output=OUTPUT, name=NAMES),
        #expand("{output}/{filter}/{name}.bam", output=OUTPUT, filter=FILTER, name=NAMES),
        expand("{output}/{filter}/{name}_lofreq.vcf", output=OUTPUT, filter=FILTER, name=NAMES),
        expand("{output}/{filter}/{name}_vardict_4_2.txt", output=OUTPUT, filter=FILTER, name=NAMES),
        expand("{output}/markedDuplicates/{name}.bam", output=OUTPUT, name=NAMES),
         
rule organise:
    input:
        fq1=f"{PREFIX}/{{n}}_R1_001.fastq",
        fq2=f"{PREFIX}/{{n}}_R2_001.fastq"
    output:
        out1="{output}/RawReads/{n}_R1_001.fastq",
        out2="{output}/RawReads/{n}_R2_001.fastq"
    shell: """ 
        cp {input.fq1} {output.out1}
        cp {input.fq2} {output.out2}
    """

rule Trim:
    input:
        fq1=rules.organise.output.out1,
        fq2=rules.organise.output.out2
    output:
        trim1="{output}/tr/{n}_1.fastq",
        trim2="{output}/tr/{n}_2.fastq"
    params:
        outdir="{output}/tr/",
        val1="{output}/tr/{n}_R1_001_val_1.fq",
        val2="{output}/tr/{n}_R2_001_val_2.fq"
    shell:"""
      trim_galore --paired {input.fq1} {input.fq2} -o {params.outdir} 
      mv {params.val1} {output.trim1}
      mv {params.val2} {output.trim2}
    """

# prefix_nodup_1.fastq and prefix_nodup_2.fastq
rule Dedup1:
    input:
        fq1=rules.Trim.output.trim1,
        fq2=rules.Trim.output.trim2
    output:
        ddfq1="{output}/tr_dd/{n}_1.fastq",
        ddfq2="{output}/tr_dd/{n}_2.fastq"
    params:
        outdir="{output}/tr_dd/",
        list=temp("{output}/tr_dd/fastuniqInput")
    shell:"""
      echo {input.fq1} > {params.list}
      echo {input.fq2} >> {params.list}
      fastuniq -i {params.list} -o {output.ddfq1} -p {output.ddfq2}
      """

# third step (Qual) - stringent quality filter 
# PRINSEQ with: 
# (i) reads retained if their average quality score was >25, 
# (ii) the 3′ end of each read trimmed if the mean quality score was <30, 
# using a sliding window (--trim_qual_window) of 5 nt and a step size (--trim_qual_step) of 1 nt, 
# (iii) homopolymeric sequences of >20 nt trimmed from the 3’ end of reads, and 
# (iv) only reads with a residual length of ≥80 nt were retained. 
# The filtered reads should be called prefix_1.fastq and 
# prefix_2.fastq but this time should be placed in a folder called 
# tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN
rule Qual:
    input:
        fq1=rules.Dedup1.output.ddfq1,
        fq2=rules.Dedup1.output.ddfq2
    output:
        prinseq1="{output}/tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN/{n}_1.fastq",
        prinseq2="{output}/tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN/{n}_2.fastq"
    params:
        outdir="{output}/tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN/{n}"
    shell:"""
      prinseq-lite.pl -fastq {input.fq1} -fastq2 {input.fq2} -out_format 3 \
        -min_qual_mean 25 \
        -trim_qual_right 30 --trim_qual_window 5 --trim_qual_step 1 \
        -trim_ns_right 20 \
        -min_len 80 \
        -out_good {params.outdir}
      head -n1 {output.prinseq1}
      head -n1 {output.prinseq2} 
    """

#Calculate SNPs for the different trimmed/filtered reads: tr tr_dd and tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN

rule bowtie2:
  input:
      fq1="{output}/{filter}/{smp}_1.fastq",
      fq2="{output}/{filter}/{smp}_2.fastq"      
  output:
      bam="{output}/{filter}/{smp}.bam"
  params:
      bowtie_index=config["index"]
  shell:"""
     echo "*************"
     echo {input.fq1} 
     echo "*************"
     echo {input.fq2} 
     echo "*************"
     echo {output.bam}
     bowtie2 --end-to-end -1 {input.fq1} -2 {input.fq2} -x {params.bowtie_index} | samtools view -h -F4 -bS | samtools sort -o {output.bam}
     samtools index {output.bam}
  """
      
rule lofreq:
  input:
      bam="{output}/{filter}/{smp}.bam"
  output:
      vcf="{output}/{filter}/{smp}_lofreq.vcf"
  params:
      ref=config["reference"]
  threads: 2
  shell:"""
     lofreq call-parallel --pp-threads {threads}  -f {params.ref} -Q 30 -q 30 -o {output.vcf} {input.bam}
  """

rule vardict:
  input:
      bam="{output}/{filter}/{smp}.bam"
  output:
      txt="{output}/{filter}/{smp}_vardict.txt",
      txt2="{output}/{filter}/{smp}_vardict_4_2.txt"
  params:
      ref=config["reference"],
      ref_id=config["ref_id"]
  threads: 2
  shell:"""
     vardict -th {threads} -G {params.ref} -f 0.001 -q 30 -N identifier -b {input.bam} -R {params.ref_id}:0-236000 > {output.txt}
     awk '$12>=2 && $13>=2 && $34==\"SNV\"' {output.txt} > {output.txt2}
  """

rule picard:
  input:
      prinseqbam="{output}/tr_dd_pr_MQ_25_TQ_30_TQW_5_TQS_1_polyN/{smp}.bam"
  output:
      rgbam=temp("{output}/markedDuplicates/{smp}_rg_added.bam"),
      markeddupbam=temp("{output}/markedDuplicates/{smp}_markedDup.bam"),
      nodupbam="{output}/markedDuplicates/{smp}.bam"
  shell:"""
      picard AddOrReplaceReadGroups I={input.prinseqbam}  O={output.rgbam} SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus
      picard  MarkDuplicates I={output.rgbam} O={output.markeddupbam} CREATE_INDEX=tr_dd_pr_polyNue VALIDATION_STRINGENCY=SILENT M=output.metr_dd_pr_polyNics
      samtools view -b -h -F 1024 {output.markeddupbam} > {output.nodupbam}
      samtools index {output.nodupbam}
  """

#Calculate SNPs for tr_dd_pr_polyN reads
# os.system("bowtie2 --end-to-end -1 "+read1_tr_dd_pr_polyN+" -2 "+read2_tr_dd_pr_polyN+" -x "+referenceIndex+" -S "+prefix+"_tr_dd_pr_polyN_alignment.sam")
# os.system("samtools view -bS -h -F 4 "+prefix+"_tr_dd_pr_polyN_alignment.sam >"+prefix +"_tr_dd_pr_polyN_alignment.bam")
# os.system("samtools sort -o "+prefix + "_tr_dd_pr_polyN_alignment_sorted.bam "+prefix +"_tr_dd_pr_polyN_alignment.bam")
# os.system("~/Software/jre1.8.0_191/bin/java -jar -XX:ParallelGCThreads=2 ~/Software/mySoftware/old/GRACy/resources/picard.jar  AddOrReplaceReadGroups I="+prefix + "_tr_dd_pr_polyN_alignment_sorted.bam  O="+prefix+"_rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=Ilumina RGPU=machine RGSM=Consensus")
# os.system("~/Software/jre1.8.0_191/bin/java -jar -XX:ParallelGCThreads=2  ~/Software/mySoftware/old/GRACy/resources/picard.jar MarkDuplicates I="+prefix+"_rg_added_sorted.bam O="+prefix+"_tr_dd_pr_polyN_markedDuplicates.bam CREATE_INDEX=tr_dd_pr_polyNue VALIDATION_STRINGENCY=SILENT M=output.metr_dd_pr_polyNics")
# os.system("samtools index "+prefix + "_tr_dd_pr_polyN_alignment_sorted.bam")
# os.system("samtools index "+prefix +"_tr_dd_pr_polyN_markedDuplicates.bam")
# os.system("samtools view -b -h -F 1024 "+prefix +"_tr_dd_pr_polyN_markedDuplicates.bam > "+prefix+"_tr_dd_pr_polyN_markedDuplicates_noDup.bam")
# os.system("samtools index "+prefix+"_tr_dd_pr_polyN_markedDuplicates_noDup.bam")
# os.system("rm -f "+prefix+"_tr_dd_pr_polyN_alignment.sam")
# os.system("rm -f "+prefix+"_tr_dd_pr_polyN_alignment.bam")
# os.system("lofreq call-parallel --pp-threads "+threads+" -f "+reference+" -Q 30 -q 30 -o "+prefix+"_tr_dd_pr_polyN_SNPs.vcf "+prefix + "_tr_dd_pr_polyN_alignment_sorted.bam")
# os.system("lofreq call-parallel --pp-threads "+threads+" -f "+reference+" -Q 30 -q 30 -o "+prefix+"_tr_dd_pr_polyN_markedDup_SNPs.vcf "+prefix +"_tr_dd_pr_polyN_markedDuplicates.bam")
# os.system("~/Software/VarDictJava/build/install/VarDict/bin/VarDict -th "+threads+" -G "+reference+" -f 0.001 -q 30 -N identifier -b "+prefix + "_tr_dd_pr_polyN_alignment_sorted.bam -R "+referenceID+":0-236000 > "+prefix+"_tr_dd_pr_polyN_noFilter.txt")
# os.system("awk '$12>=2 && $13>=2 && $34==\"SNV\"' "+prefix+"_tr_dd_pr_polyN_noFilter.txt > "+prefix+"_tr_dd_pr_polyN_noFilter_4_2.txt")
# os.system("~/Software/VarDictJava/build/install/VarDict/bin/VarDict -th "+threads+" -G "+reference+" -f 0.001 -q 30 -N identifier -b "+prefix + "_tr_dd_pr_polyN_markedDuplicates_noDup.bam -R "+referenceID+":0-236000 > "+prefix+"_tr_dd_pr_polyN_noDup_noFilter.txt")
# os.system("awk '$12>=2 && $13>=2 && $34==\"SNV\"' "+prefix+"_tr_dd_pr_polyN_noDup_noFilter.txt > "+prefix+"_tr_dd_pr_polyN_noDup_noFilter_4_2.txt")
# os.system("rm -f "+prefix+"_tr_dd_pr_polyN_alignment_sorted.bam")
# os.system("rm -f "+prefix+"_tr_dd_pr_polyN_markedDuplicates.bam")
# os.system("rm -f "+prefix+"_tr_dd_pr_polyN_markedDuplicates_noDup.bam")
