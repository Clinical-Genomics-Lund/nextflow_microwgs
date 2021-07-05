#!/usr/bin/env nextflow

OUTDIR = params.outdir+'/'+params.subdir

// Get CSV input
//csv = file(params.csv)
Channel
    .fromPath(params.csv).splitCsv(header:true)
    .map{ row-> tuple(row.id, row.species, row.platform, file(row.read1), file(row.read2)) }
    .into { fastq_bwa; fastq_spades; fastq_kraken; fastq_ariba; fastq_maskpolymorph; fastq_register; fastq_cgviz }


process bwa_align {
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_bwa

	output:
		set id, species, platform, file("${id}.sam") into bwa_sam

	script:
		fasta_ref = params.refpath+'/species/'+species+'/ref.fasta'
		read2 = fastq_r2.name == 'SINGLE_END' ? '' : "$fastq_r2"

	"""
	bwa mem  \\
	-R '@RG\\tID:${id}\\tSM:${id}\\tPL:illumina' \\
	-M \\
	-t ${task.cpus} \\
	$fasta_ref $fastq_r1 $read2 > ${id}.sam
	"""
}


process convert_to_bam {
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file("${id}.sam") from bwa_sam

	output:
		set id, species, platform, file("${id}.bwa.sort.bam"), file("${id}.bwa.sort.bam.bai") into bam_markdup, bam_postalignqc
	"""
	samtools view -Sb ${id}.sam \\
		| samtools sort -o ${id}.bwa.sort.bam -
	samtools index ${id}.bwa.sort.bam
	"""
}


process bam_markdup {
	publishDir "${OUTDIR}/bam", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '32 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(bam), file(bai) from bam_markdup

	output:
		set id, species, platform, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") into bam_qc

	"""
	sambamba markdup --tmpdir tmp -t ${task.cpus} $bam ${id}.dedup.bam
	"""
}


process kraken {
	publishDir "${OUTDIR}/kraken", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '48 GB'
	time '2h'
	tag "$id"

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_kraken

	output:
		set id, species, platform, file("${id}.kraken.report") into kraken_output

	script:
		read_params = fastq_r2.name == 'SINGLE_END' ?
			"$fastq_r1" :
			"--paired $fastq_r1 $fastq_r2"

	"""
	kraken2 \\
		--gzip-compressed \\
		--db ${params.krakendb} \\
		--threads ${task.cpus} \\
		--output ${id}.kraken.out \\
		--report ${id}.kraken.report \\
		$read_params
	"""
}


process bracken {
	publishDir "${OUTDIR}/kraken", mode: 'copy', overwrite: true
	cpus params.cpu_many
	memory '48 GB'
	time '2h'
	tag "$id"

	input:
		set id, species, platform, file("${id}.kraken.report") from kraken_output

	output:
		set id, species, platform, file("${id}.bracken") into kraken_export

	"""
	bracken -d ${params.krakendb} -r 150 -i ${id}.kraken.report -o ${id}.bracken
	"""
}

process spades_assembly {
	publishDir "${OUTDIR}/assembly", mode: 'copy', overwrite: true
	cpus params.cpu_spades
	memory '64 GB'
	time '2h'
	tag "$id"

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_spades

	output:
		set id, species, platform, file("${id}.fasta") into asm_quast, asm_mlst, asm_chewbbaca, asm_maskpolymorph

	script:
		opt_platform = platform == 'iontorrent' ? '--iontorrent --careful --sc' : '--only-assembler'
		opt_reads = fastq_r2.name != 'SINGLE_END' ? "-1 $fastq_r1 -2 $fastq_r2" : "-s $fastq_r1"

	"""
	spades.py -k 21,33,55,77 -t ${task.cpus} \\
		$opt_platform \\
		$opt_reads \\
		-o spades
	mv spades/contigs.fasta ${id}.fasta
	"""
}


process quast {
	publishDir "${OUTDIR}/qc", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta) from asm_quast

	output:
		set id, species, platform, file("${id}.quast.tsv") into quast_register, quast_export

	script:
		fasta_ref = params.refpath+'/species/'+species+'/ref.fasta'

	"""
	quast.py $asm_fasta -R $fasta_ref -o quast_outdir
	cp quast_outdir/transposed_report.tsv ${id}.quast.tsv
	"""
}


process mlst {
	publishDir "${OUTDIR}/mlst", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta) from asm_mlst

	output:
		set id, species, platform, file("${id}.mlst.json") into mlst_export
		file("${id}.mlst.novel") optional true


	"""
	mlst --scheme ${species} \\
		--json ${id}.mlst.json --novel ${id}.mlst.novel \\
		${asm_fasta}
	"""
}


process ariba {
	publishDir "${OUTDIR}/ariba", mode: 'copy', overwrite: true
	// cache 'deep'
	cpus params.cpu_many
	memory '16 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2) from fastq_ariba

	output:
		set id, species, platform, file("${id}.ariba.summary.csv"), file("${id}_report.tsv") into ariba_report

	when:
		fastq_r2.toString() != "SINGLE_END"

	"""
	ariba run --force --threads ${task.cpus} ${params.refpath}/species/${species}/ariba \\
		$fastq_r1 $fastq_r2 ariba.outdir
	# rename report
	cp ariba.outdir/report.tsv ${id}_report.tsv
	# make summary
	ariba summary --col_filter n --row_filter n "${id}.ariba.summary" ${id}_report.tsv
	"""
}


process ariba_export {
	publishDir "${OUTDIR}/ariba", mode: 'copy', overwrite: true
	// cache 'deep'
	cpus params.cpu_many
	memory '1 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file("${id}.ariba.summary.csv"), file("${id}_report.tsv") from ariba_report

	output:
		set id, species, platform, file("${id}.ariba.json") into ariba_export

	"""
	ariba2json.pl ${params.refpath}/species/${species}/ariba/02.cdhit.all.fa \\
		 "${id}.ariba.summary.csv" ${id}_report.tsv > ${id}.ariba.json
	"""
}


process bwa_maskpolymorph {
	publishDir "${OUTDIR}/maskepolymorph", mode: 'copy', overwrite: true
	cpus params.cpu_bwa
	memory '32 GB'
	time '1h'
	// cache 'deep'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta), file(fastq_r1), file(fastq_r2) from asm_maskpolymorph.join(fastq_maskpolymorph, by:[0,1,2])

	output:
		set id, species, platform, file(asm_fasta), file("${id}_contigs.fasta.sort.sam") into maskpoly_sam

	script:
		read2 = fastq_r2.name == 'SINGLE_END' ? '' : "$fastq_r2"
		faidx = asm_fasta+'.fai'

	"""
	bwa index $asm_fasta
	bwa mem -t ${task.cpus} -f ${id}_contigs.fasta.sort.sam -R '@RG\\tID:${id}\\tSM:${id}\\tPL:${platform}' $asm_fasta $fastq_r1 $read2
	"""
}


process samtools_maskpolymorph {
	publishDir "${OUTDIR}/maskepolymorph", mode: 'copy', overwrite: true
	memory '2 GB'
	time '1h'
	// cache 'deep'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta), file("${id}_contigs.fasta.sort.sam") from maskpoly_sam

	output:
		set id, species, platform, file(asm_fasta), file(asm_fasta_index), file("${id}_contigs.fasta.sort.bam") into maskpoly_bam
	script:
	    asm_fasta_index = asm_fasta+'.fai'

	"""
	samtools view -b -@ ${task.cpus} ${id}_contigs.fasta.sort.sam \\
	    | samtools sort -@ ${task.cpus} - -o "${id}_contigs.fasta.sort.bam"
	samtools index "${id}_contigs.fasta.sort.bam"
	samtools faidx "$asm_fasta"
	"""
}


process freebayes_maskpolymorph {
	publishDir "${OUTDIR}/maskepolymorph", mode: 'copy', overwrite: true
	memory '8 GB'
	time '1h'
	// cache 'deep'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta), file(asm_fasta_index) ,file("${id}_contigs.fasta.sort.bam") from maskpoly_bam

	output:
		set id, species, platform, file(asm_fasta), file("${id}_contigs.fasta.vcf") into maskpoly_vcf
	"""
	freebayes -f $asm_fasta "${id}_contigs.fasta.sort.bam" -C 2 -F 0.2 --pooled-continuous > "${id}_contigs.fasta.vcf"
	"""
}


process maskpolymorph_error_corr {
	publishDir "${OUTDIR}/maskepolymorph", mode: 'copy', overwrite: true
	memory '1 GB'
	time '1h'
	// cache 'deep'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta), file("${id}_contigs.fasta.vcf") from maskpoly_vcf

	output:
		set id, species, platform, file("${id}.spades.masked.fasta") into maskpoly_chewbbaca
	"""
	error_corr_assembly.pl $asm_fasta ${id}_contigs.fasta.vcf > ${id}.spades.masked.fasta
	"""
}


process chewbbaca {
	publishDir "${OUTDIR}/chewbbaca", mode: 'copy', overwrite: true, pattern: '*.chewbbaca'
	cpus params.cpu_chewbbaca
	memory '8 GB'
	time '1h'
	// cache 'deep'
	tag "$id"

	input:
		set id, species, platform, file(asm_fasta) from maskpoly_chewbbaca.groupTuple( by: 1)
	output:
		file("results_alleles.tsv") into chewbbaca_export
		file(missingloci) into chewbbaca_register
	script:
		missingloci = "chewbbaca.missingloci"
	"""
	ls $asm_fasta > batch_input.list
	flock -e ${params.local_tmp}/chewbbaca.lock \\
		  chewBBACA.py AlleleCall --fr -i ./ -g $params.cgmlst_db --ptf Staphylococcus_aureus.trn --cpu ${task.cpus} \\
		  -o chewbbaca.folder

	cp chewbbaca.folder/results_*/results_alleles.tsv results_alleles.tsv > chewbbaca.missingloci
	bash parse_missing_loci.sh batch_input.list results_alleles.tsv $missingloci
	"""
}


process split_chewbbaca_results {
	publishDir "${OUTDIR}/chewbbaca", mode: 'copy', overwrite: true, pattern: '*.chewbbaca'
	cpus 1
	memory '1 GB'
	time '1h'
	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2), file(chewbbaca) from fastq_cgviz.combine(chewbbaca_export)
	output:
		set id, species, platform, file(fastq_r1), file(fastq_r2), file(chewbacca_result) into chewbacca_split
	script:
		chewbacca_result = id+'.chewbbaca'
	"""
	head -1 ${chewbbaca} > ${chewbacca_result}
	grep ${id} ${chewbbaca} >> ${chewbacca_result}
	"""
}


process split_chewbbaca_missing_loci {
	publishDir "${OUTDIR}/chewbbaca", mode: 'copy', overwrite: true, pattern: '*.chewbbaca'
	cpus 1
	memory '1 GB'
	time '1h'
	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2), file(chewbbaca) from fastq_register.combine(chewbbaca_register)
	output:
		set id, species, platform, file(fastq_r1), file(fastq_r2), file(chewbacca_loci) into chewbacca_missing_loci
	script:
		chewbacca_loci = id+'.chewbbaca'
	"""
	grep ${id} ${chewbbaca} > ${chewbacca_loci}
	"""
}


process postalignqc {
	publishDir "${OUTDIR}/postalignqc", mode: 'copy', overwrite: true
	cpus 4
	memory '8 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, file(bam), file(bai) from bam_postalignqc

	output:
		set id, species, platform, file("${id}.bwa.QC") into postqc_register

	script:
		cgmlst_bed = params.refpath+'/species/'+species+'/cgmlst.bed'

	"""
	postaln_qc.pl $bam $cgmlst_bed ${id} ${task.cpus} > ${id}.bwa.QC

	"""
}


process to_cdm {
	publishDir "${params.crondir}/qc", mode: 'copy', overwrite: true
	cpus 1
	memory '8 GB'
	time '1h'
	tag "$id"

	input:
		set id, species, platform, fastq_r1, fastq_r2, \
			file(missingloci), \
			file(quast), \
			file(postalignqc) \
		from chewbacca_missing_loci\
			 .join(quast_register,by:[0,1,2])\
			 .join(postqc_register,by:[0,1,2])

	output:
		set id, species, platform, file("${id}.cdm")
		set id, species, platform, rundir into register_export

	when:
		!params.test

	script:
		parts = fastq_r1.toString().split('/')
		parts.println()
		idx =  parts.findIndexOf {it ==~ /\d{6}_.{6,8}_.{4}_.{10}/}
		rundir = parts[0..idx].join("/")
		postalignqc = params.outdir+'/'+params.subdir+'/postalignqc/'+postalignqc
		quast  = params.outdir+'/'+params.subdir+'/qc/'+quast

	"""
	read missingloci < $missingloci
	echo --run-folder ${rundir} \\
         --sample-id ${id} \\
         --assay microbiology \\
         --qc ${postalignqc} \\
         --asmqc $quast \\
         --micmisloc \$missingloci > ${id}.cdm
	"""
}


process to_cgviz {
	publishDir "${params.crondir}/cgviz", mode: 'copy', overwrite: true
	cpus 1
	memory '1 GB'
	time '1h'
	input:
		set id, species, platform, file(fastq_r1), file(fastq_r2), \
		    file(chewbbaca),  \
			file(quast),  \
			file(mlst),   \
			file(kraken), \
			file(ariba)   \
		from chewbacca_split\
			.join(quast_export,by:[0,1,2])\
			.join(mlst_export,by:[0,1,2])\
			.join(kraken_export,by:[0,1,2])\
			.join(ariba_export,by:[0,1,2])
	output:
		set species, platform, file("${id}.cgviz")
	when:
		params.cgviz

	script:
		parts = fastq_r1.toString().split('/')
		idx =  parts.findIndexOf {it ==~ /\d{6}_.{6,8}_.{4}_.{10}/}
		rundir = parts[0..idx].join("/")

		chewbbaca = params.outdir+'/'+params.subdir+'/chewbbaca/'+chewbbaca
		quast = params.outdir+'/'+params.subdir+'/qc/'+quast
		mlst = params.outdir+'/'+params.subdir+'/mlst/'+mlst
		kraken = params.outdir+'/'+params.subdir+'/kraken/'+kraken
		ariba = params.outdir+'/'+params.subdir+'/ariba/'+ariba
		if( params.test ) {
			species = "test"
		}

	"""
	echo "--overwrite \\
		 --in $chewbbaca  \\
		 --id $id \\
		 --species $species \\
		 --run $rundir \\
		 --quast $quast \\
		 --mlst $mlst \\
		 --kraken $kraken \\
		 --aribavir $ariba" > ${id}.cgviz
	"""
}
