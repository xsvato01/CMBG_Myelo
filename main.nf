nextflow.enable.dsl = 2
launchDir = "${launchDir}/${params.datain.replaceAll(".*/","")}"

process TRIMMING {
	tag "trimming on $name using $task.cpus CPUs and $task.memory memory"
	publishDir  "${launchDir}/trimmed/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
	tuple val(name), path("*.fastq.gz")

	script:
	""" 
	cutadapt -m 50 -o ${name}.trimmed.R1.fastq.gz -p ${name}.trimmed.R2.fastq.gz $reads
	"""
}

process FIRST_ALIGN_BAM {
	tag "first align on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/mapped/", mode:'copy'
	
	input:
	tuple val(name), path(reads)

	output:
        tuple val(name), path("${name}.sorted.bam")
	tuple val(name), path("${name}.sorted.bai")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	bwa mem -R ${rg} -t 4 ${params.refindex} $reads \
	| samtools view -Sb -o - -| samtools sort -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai	
	"""
}

process FIRST_QC {
	tag "first QC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/first_bam_qc/", mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	path "*"

	script:
	"""
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
        picard BedToIntervalList I=${params.covbed} O=${name}.interval_list SD=${params.ref}.dict
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=${name}.interval_list TARGET_INTERVALS=${name}.interval_list R=${params.ref}.fa O=${name}.aln_metrics
	"""
}

process MARK_DUPLICATES {
	tag "Mark duplicates on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/first_bam_qc/", pattern: '*.txt', mode:'copy'
	publishDir "${launchDir}/mapped/", pattern: '*.md.ba*', mode:'copy'
	
	input:
	tuple val(name), path(bam)

	output:
	path "*.txt"
	tuple val(name), path("*first.md.bam")
	path "*.bai"

	script:
	"""
	picard MarkDuplicates I=$bam M=${name}.MD.metrics.txt O=${name}.first.md.bam
	samtools index ${name}.first.md.bam
	"""
}

process MULTIQC {
	tag "MultiQC on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/multiqc_reports/", mode:'copy'
	
	input:
	path '*'

	output:
	path 'first_report.html'

	script:
	"""
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
	multiqc . -n first_report.html
	"""
}

process MUTECT2 {
	tag "MUTECT2 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	
	input:
	tuple val(name), path(bam)
	
	output:
	tuple val(name), path ('*.vcf')

	script:
	"""
	gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --annotation StrandArtifact --min-base-quality-score 30 --intervals $params.varbed  --output ${name}.mutect.vcf
	"""

}

process FILTER_MUTECT {
	tag "filter mutect on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	
	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path ('*.vcf')

	script:
	"""
	gatk FilterMutectCalls -V $vcf_input -O ${name}.mutect.filt.vcf
	"""

}

process NORMALIZE_MUTECT {
	tag "normalize filtered mutect on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	
	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path ('*.vcf')

	script:
	"""
	bcftools norm -f ${params.ref}.fa $vcf_input -o ${name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "annotate mutect on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	
	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path('*.vcf')
	path '*.csv'
	script:
	"""
	vep -i $vcf_input --cache --cache_version 90 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --everything -o ${name}.mutect2.filt.norm.vep.vcf

	bcftools view -f 'PASS,clustered_events' ${name}.mutect2.filt.norm.vep.vcf \
	| python $params.vcftbl simple --build GRCh38 -i /dev/stdin -o ${name}.mutect2.filt.norm.vep.csv

	"""	
}

process CREATE_FULL_TABLE {
	tag "creating full table on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	
	input:
	tuple val(name), path(vcf_input)
	
	output:
	path '*'

	script:
	"""
	python $params.vcftbl simple --build GRCh38 -i $vcf_input -o ${name}.mutect2.filt.norm.vep.full.csv
	"""	
}

process COVERAGE {
	tag "calculating coverage on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(bam)
	
	output:
	tuple val(name), path('*.PBcov.txt')

	script:
	"""
	bedtools coverage -abam $params.covbed -b $bam -d > ${name}.PBcov.txt
	"""	
}


process COVERAGE_R {
	tag "R coverage on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/coverage/", mode:'copy'
	
	input:
	tuple val(name), path(pbcov)
	
	output:
	path '*'

	script:
	"""
	Rscript --vanilla $params.coverstat $pbcov >> ${name}_CXCR4.txt
	"""
}


 
workflow {
 rawfastq = channel.fromFilePairs("${params.datain}/*R{1,2}*", checkIfExists: true)
	
	trimmed		= TRIMMING(rawfastq)
	sortedbam	= FIRST_ALIGN_BAM(trimmed)
	qc_files	= FIRST_QC(sortedbam[0])
	qcdup_file	= MARK_DUPLICATES(sortedbam[0])
	MULTIQC(qc_files.mix(qcdup_file[0]).collect())
//qc_files.mix(qcdup_file[0]).collect().view()
//qcdup_file[1].view()
	raw_vcf         = MUTECT2(qcdup_file[1]) //markdup.bam 
    filtered        = FILTER_MUTECT(raw_vcf)
    normalized      = NORMALIZE_MUTECT(filtered)
    annotated       = ANNOTATE_MUTECT(normalized)
    CREATE_FULL_TABLE(normalized[0])
    pbcov           = COVERAGE(sortedbam[0])
pbcov.view()   
 COVERAGE_R(pbcov)
}
