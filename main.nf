run = "${params.datain}".split("/")
run = run[run.size()-1]
//launchDir = "${params.outDirectory}/${run}"
launchDir = "${params.outDirectory}/${params.runname}"

process COLLECT_BASECALLED {
	tag "COLLECT_BASECALLED on $name using $task.cpus CPUs and $task.memory memory"
	label "small_process"
	debug true

	input:
	val "name"

	output:
	tuple val(name), path("*.fastq.gz")

	script:
	"""
 echo $name
	echo ${params.runname}
	## this might cause an error if there are several basecalled folders:
	cp  /mnt/shared/MedGen/sequencing_results/primary_data/*${params.runname}/raw_fastq/${name}*R{1,2}* ./
	ls -al
	"""
} 


process TRIMMING {
	tag "trimming on $name using $task.cpus CPUs and $task.memory memory"
	label "small_process"
	
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
	label "medium_cpu"

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
	label "smallest_process"
	container 'registry.gitlab.ics.muni.cz:443/450402/btk_k8s:16'	
	
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
	label "small_process"

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
	label "smallest_process"

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
	
	input:
	tuple val(name), path(bam)
	
	output:
	tuple val(name), path ('*.vcf')
	path '*'

	script:
	"""
	gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --annotation StrandArtifact --min-base-quality-score 10 --intervals $params.covbed -bamout ${name}.bamout.bam --output ${name}.mutect.vcf
	"""

}

process FILTER_MUTECT {
	tag "filter mutect on $name using $task.cpus CPUs and $task.memory memory"
	label "smallest_process"

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
	label "smallest_process"

	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path ('*.vcf')

	script:
	"""
	bcftools norm -m-both $vcf_input > ${name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "annotate mutect on $name using $task.cpus CPUs $task.memory"
 container "ensemblorg/ensembl-vep:release_108.0"
	label "smallest_process"

	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path('*.vcf')

	script:
	"""
	vep -i $vcf_input --cache --cache_version 108 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --mane_select --offline --vcf --everything -o ${name}.mutect2.filt.norm.vep.vcf
	"""	
}

process JOIN_VARS {
	tag "JOIN_VARS on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/joined/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(name), path(vcf), path(toMerge)
	
	output:
	tuple val(name), path("${name}.joinedvariants.tsv")

	script:
	"""
 for vcf_file in $toMerge; do bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%SAMPLE]\\n' "\$vcf_file" >> ${name}.joinedvariants.tsv; done
	"""	
}

process CREATE_FULL_TABLE {
	tag "creating full table on $name using $task.cpus CPUs $task.memory"
	label "smallest_process"

	input:
	tuple val(name), path(vcf_input)
	
	output:
	tuple val(name), path("${name}.mutect2.filt.norm.vep.full.csv")

	script:
	"""
	python $params.vcftbl simple --build GRCh38 -i $vcf_input -t ${name} -o ${name}.mutect2.filt.norm.vep.full.csv
	"""	
}


process ALTER_FULL_TABLE {
	tag "ALTER_FULL_TABLE on $name using $task.cpus CPUs $task.memory"
	publishDir "${launchDir}/variants/", mode:'copy'
	label "smallest_process"

	input:
	tuple val(name), path(vcf), path(joinedTsv), val(ntotal)
	
	output:
	tuple val(name), path("${name}.final.csv")

	script:
	"""
	echo alter full table on $name
	python ${params.mergetables} --table $vcf --varlist $joinedTsv --n $ntotal --outname ${name}.final.csv
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
	label "smallest_process"

	input:
	tuple val(name), path(pbcov)

	script:
	"""
	Rscript --vanilla $params.coverstat $pbcov ${launchDir}/coverage/${name}.perexon_stat.txt
	"""
}
 
workflow {
	runlist = channel.fromPath(params.samplesheet).splitCsv().flatten().view{"runlist: $it"}
 rawfastq = COLLECT_BASECALLED(runlist).view{"rawfastq: $it"}

 trimmed		= TRIMMING(rawfastq)
 sortedbam	= FIRST_ALIGN_BAM(trimmed)
 qc_files	= FIRST_QC(sortedbam[0])
 qcdup_file	= MARK_DUPLICATES(sortedbam[0])
 MULTIQC(qc_files.collect())
 raw_vcf         = MUTECT2(qcdup_file[1]) //markdup.bam 
 filtered        = FILTER_MUTECT(raw_vcf[0])
 normalized      = NORMALIZE_MUTECT(filtered)
 annotated       = ANNOTATE_MUTECT(normalized)
 full_table = CREATE_FULL_TABLE(annotated)
	Vcf_paths = normalized.map({it -> [it[1]]})
	Combined_collected_vcf = normalized.combine(Vcf_paths.collect().map({it -> [it]}))
 Combined_filtered = Combined_collected_vcf.map({
		row ->
		def name = row[0]
		def vcf = row[1]
		def filtered = removeSame(vcf, row[2])
		[name,vcf, filtered]	
	})//.view{"____________Combined_filtered____________: $it"}

 joined_vars = JOIN_VARS(Combined_filtered)
	joined_total = joined_vars.combine(joined_vars.count()-1)
 ALTER_FULL_TABLE(full_table.join(joined_total))

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