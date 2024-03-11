process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"

	input:
	val(sample)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	"""
	echo COLLECT_BASECALLED $sample.name
	cp  /mnt/shared/MedGen/sequencing_results/primary_data/*${sample.run}/raw_fastq/${sample.name}*R{1,2}* ./
	"""
} 


process TRIMMING {
	tag "trimming on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"
	
	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("*.fastq.gz")

	script:
	""" 
	cutadapt -m 50 -o ${sample.name}.trimmed.R1.fastq.gz -p ${sample.name}.trimmed.R2.fastq.gz $reads
	"""
}

process FIRST_ALIGN_BAM {
	tag "first align on $sample.name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/mapped/", mode:'copy'
	label "medium_cpu"

	input:
	tuple val(sample), path(reads)

	output:
	tuple val(sample), path("${sample.name}.sorted.bam")
	tuple val(sample), path("${sample.name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${sample.name}\\tSM:${sample.name}\\tLB:${sample.name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4 ${params.refindex} $reads \
	| samtools view -Sb -o - - | samtools sort -o ${sample.name}.sorted.bam
	samtools index ${sample.name}.sorted.bam ${sample.name}.sorted.bai	
	"""
}

process FIRST_QC {
	tag "first QC on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "smallest_process"
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	
	
	input:
	tuple val(sample), path(bam)

	output:
	path "*"

	script:
	"""
	samtools flagstat $bam > ${sample.name}.flagstat
	samtools stats $bam > ${sample.name}.samstats
	picard BedToIntervalList I=${params.covbed} O=${sample.name}.interval_list SD=${params.ref}.dict
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=${sample.name}.interval_list TARGET_INTERVALS=${sample.name}.interval_list R=${params.ref}.fa O=${sample.name}.aln_metrics
	"""
}

process MARK_DUPLICATES {
	tag "Mark duplicates on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "small_process"

	input:
	tuple val(sample), path(bam)

	output:
	path "*.txt"
	tuple val(sample), path("*first.md.bam")
	path "*.bai"

	script:
	"""
	picard MarkDuplicates I=$bam M=${sample.name}.MD.metrics.txt O=${sample.name}.first.md.bam
	samtools index ${sample.name}.first.md.bam
	"""
}

process MULTIQC {
	tag "MultiQC using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/multiqc_reports/", mode:'copy'
	label "smallest_process"

	input:
	path '*'

	output:
	path '*.html'

	script:
	"""
	export LC_ALL=C.UTF-8
	export LANG=C.UTF-8
	multiqc . -n MultiQC-"`date +"%d-%m-%Y"`".html
	"""
}

process MUTECT2 {
	tag "MUTECT2 on $sample.name using $task.cpus CPUs and $task.memory memory"
	
	input:
	tuple val(sample), path(bam)
	
	output:
	tuple val(sample), path ('*.vcf')
	path '*'

	script:
	"""
	gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --annotation StrandArtifact --min-base-quality-score 10 --intervals $params.covbed -bamout ${sample.name}.bamout.bam --output ${sample.name}.mutect.vcf
	"""
}

process FILTER_MUTECT {
	tag "filter mutect on $sample.name using $task.cpus CPUs and $task.memory memory"
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path ('*.vcf')

	script:
	"""
	gatk FilterMutectCalls -V $vcf_input -O ${sample.name}.mutect.filt.vcf
	"""

}

process NORMALIZE_MUTECT {
	tag "normalize filtered mutect on $sample.name using $task.cpus CPUs $task.memory"
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path ('*.vcf')

	script:
	"""
	bcftools norm -m-both $vcf_input > ${sample.name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "annotate mutect on $sample.name using $task.cpus CPUs $task.memory"
	container "ensemblorg/ensembl-vep:release_108.0"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path('*.vcf')

	script:
	"""
	vep -i $vcf_input --cache --cache_version 108 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --mane_select --offline --vcf --everything -o ${sample.name}.mutect2.filt.norm.vep.vcf
	"""	
}

process JOIN_VARS {
	tag "JOIN_VARS on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/joined/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(vcf), path(toMerge)
	
	output:
	tuple val(sample), path("${sample.name}.joinedvariants.tsv")

	script:
	"""
	for vcf_file in $toMerge; do bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE]\\n' "\$vcf_file" >> ${sample.name}.joinedvariants.tsv; done
	"""	
}

process JOIN_VARS_ALL {
	tag "JOIN_VARS_ALL  using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/joined/", mode:'copy'
	label "smallest_process"

	input:
	path "VcfsToMerge"
	
	output:
	path "joinedAllVariants.tsv"

	script:
	"""
	for vcf_file in $VcfsToMerge; do bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE]\\n' "\$vcf_file" >> joinedAllVariants.tsv; done
	"""	
}

process CREATE_FULL_TABLE {
	tag "creating full table on $sample.name using $task.cpus CPUs $task.memory"
		publishDir "${params.outDirectory}/${sample.run}/create_full_table/", mode:'copy'

	label "smallest_process"

	input:
	tuple val(sample), path(vcf_input)
	
	output:
	tuple val(sample), path("${sample.name}.mutect2.filt.norm.vep.full.csv")

	script:
	"""
	python $params.vcftbl simple --build GRCh38 -i $vcf_input -t ${sample.name} -o ${sample.name}.mutect2.filt.norm.vep.full.csv
	"""	
}


process ALTER_FULL_TABLE {
	tag "ALTER_FULL_TABLE on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(variantsTableCsv), path(joinedTsv), val(ntotal)
	
	output:
	tuple val(sample), path("${sample.name}.final.csv")

	script:
	"""
	echo alter full table on $sample.name
	python ${params.mergetables} --table $variantsTableCsv --varlist $joinedTsv --n $ntotal --outname ${sample.name}.final.csv
	"""	
}

process COVERAGE {
	tag "calculating coverage on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'

	input:
	tuple val(sample), path(bam)
	
	output:
	tuple val(sample), path('*.PBcov.txt')

	script:
	"""
	bedtools coverage -abam $params.covbed -b $bam -d > ${sample.name}.PBcov.txt
	"""	
}



process COVERAGE_R {
	tag "R coverage on $sample.name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/coverage/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(sample), path(pbcov)

	script:
	"""
	Rscript --vanilla $params.coverstat $pbcov ${params.outDirectory}/${sample.run}/coverage/${sample.name}.perexon_stat.txt
	"""
}

workflow {

runlist = channel.fromList(params.samples)
rawfastq = COLLECT_BASECALLED(runlist)

trimmed	= TRIMMING(rawfastq)
sortedbam	= FIRST_ALIGN_BAM(trimmed)
qc_files	= FIRST_QC(sortedbam[0])
qcdup_file	= MARK_DUPLICATES(sortedbam[0])
MULTIQC(qc_files.collect())
raw_vcf	= MUTECT2(qcdup_file[1]) //markdup.bam 
filtered	= FILTER_MUTECT(raw_vcf[0])
normalized	= NORMALIZE_MUTECT(filtered)
annotated	= ANNOTATE_MUTECT(normalized)
full_table	= CREATE_FULL_TABLE(annotated)
Vcf_paths = normalized.map({it -> [it[1]]})
Vcf_paths_collected = Vcf_paths.collect()
Combined_collected_vcf = normalized.combine(Vcf_paths.collect().map({it -> [it]}))
Combined_filtered = Combined_collected_vcf.map({
	row ->
	def name = row[0]
	def vcf = row[1]
	def filtered = removeSame(vcf, row[2])
	[name,vcf, filtered]	
	})
Joined_all_vars = JOIN_VARS_ALL(Vcf_paths_collected)
Joined__all_total = full_table.combine(Joined_all_vars.combine(full_table.count()))
ALTER_FULL_TABLE(Joined__all_total)

pbcov = COVERAGE(sortedbam[0])
COVERAGE_R(pbcov)
}


def removeSame(nm,lst) {
	def list=[]
	for (int i = 0; i < lst.size(); i++) {
		if (lst[i] != nm) {
			list.add(lst[i])
			}
		}
	return(list)
 }
