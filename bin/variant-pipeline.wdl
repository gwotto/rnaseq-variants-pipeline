workflow RNAseq {

  File inputBam
  File inputBamIndex
  String sampleName = basename(inputBam,".bam")
  
  File refFasta
  File refFastaIndex
  File refDict

  String gatk_path
  ## String? gatk4_docker_override
  ## String gatk4_docker = select_first([gatk4_docker_override, "broadinstitute/gatk:latest"])
  ## String? gatk_path_override
  ## String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])

  String java_opt
  
  Array[File] knownVcfs
  Array[File] knownVcfsIndices

  File dbSnpVcf
  File dbSnpVcfIndex

  Int? minConfidenceForVariantCalling

  File annotationsGTF
  ## File interval_list
  
  ## Optional user optimizations
  Int? haplotypeScatterCount
  Int scatterCount = select_first([haplotypeScatterCount, 6])

  Int? preemptible_tries
  Int preemptible_count = select_first([preemptible_tries, 3])


  call gtfToCallingIntervals {
    input:
    gtf = annotationsGTF,
    ref_dict = refDict,
    preemptible_count = preemptible_count,
    gatk_path = gatk_path
    ## docker = gatk4_docker
  }

  
  call SplitNCigarReads {
    input:

    input_bam = inputBam,
    input_bam_index = inputBamIndex,
    base_name = sampleName + ".split",
    ref_fasta = refFasta,
    ref_fasta_index = refFastaIndex,
    ref_dict = refDict,
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }


  call BaseRecalibrator {
    input:
    input_bam = SplitNCigarReads.output_bam,
    input_bam_index = SplitNCigarReads.output_bam_index,
    recal_output_file = sampleName + ".recal_data.csv",
    dbSNP_vcf = dbSnpVcf,
    dbSNP_vcf_index = dbSnpVcfIndex,
    known_indels_sites_VCFs = knownVcfs,
    known_indels_sites_indices = knownVcfsIndices,
    ref_dict = refDict,
    ref_fasta = refFasta,
    ref_fasta_index = refFastaIndex,
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }

  call ApplyBQSR {
    input:
    input_bam =  SplitNCigarReads.output_bam,
    input_bam_index = SplitNCigarReads.output_bam_index,
    base_name = sampleName + ".aligned.duplicates_marked.recalibrated",
    ref_fasta = refFasta,
    ref_fasta_index = refFastaIndex,
    ref_dict = refDict,
    recalibration_report = BaseRecalibrator.recalibration_report,
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }


  call ScatterIntervalList {
    input:
    interval_list = gtfToCallingIntervals.interval_list,
    scatter_count = scatterCount,
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }


  scatter (interval in ScatterIntervalList.out) {
    call HaplotypeCaller {
      input:
      input_bam = ApplyBQSR.output_bam,
      input_bam_index = ApplyBQSR.output_bam_index,
      base_name = sampleName + ".hc",
      interval_list = interval,
      ref_fasta = refFasta,
      ref_fasta_index = refFastaIndex,
      ref_dict = refDict,
      dbSNP_vcf = dbSnpVcf,
      dbSNP_vcf_index = dbSnpVcfIndex,
      stand_call_conf = minConfidenceForVariantCalling,
      preemptible_count = preemptible_count,
      ## docker = gatk4_docker,
      gatk_path = gatk_path
    }

    File HaplotypeCallerOutputVcf = HaplotypeCaller.output_vcf
    File HaplotypeCallerOutputVcfIndex = HaplotypeCaller.output_vcf_index
  }


  call MergeVCFs {
    input:
    input_vcfs = HaplotypeCallerOutputVcf,
    input_vcfs_indexes =  HaplotypeCallerOutputVcfIndex,
    output_vcf_name = sampleName + ".g.vcf.gz",
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }
  

  call VariantFiltration {
    input:
    input_vcf = MergeVCFs.output_vcf,
    input_vcf_index = MergeVCFs.output_vcf_index,
    base_name = sampleName + ".variant_filtered.vcf.gz",
    ref_fasta = refFasta,
    ref_fasta_index = refFastaIndex,
    ref_dict = refDict,
    preemptible_count = preemptible_count,
    ## docker = gatk4_docker,
    gatk_path = gatk_path
  }

  output {
  File recalibrated_bam = ApplyBQSR.output_bam
  File recalibrated_bam_index = ApplyBQSR.output_bam_index
  File merged_vcf = MergeVCFs.output_vcf
  File merged_vcf_index = MergeVCFs.output_vcf_index
  File variant_filtered_vcf = VariantFiltration.output_vcf
  File variant_filtered_vcf_index = VariantFiltration.output_vcf_index
  }
	
}


task gtfToCallingIntervals {
  File gtf
  File ref_dict
  
  String output_name = basename(gtf, ".gtf") + ".exons.interval_list"

  ## String docker
  String gatk_path
  Int preemptible_count

  command <<<

  set -e

  Rscript --no-save -<<'RCODE'
            gtf = read.table("${gtf}", sep="\t")
            gtf = subset(gtf, V3 == "exon")
            write.table(data.frame(chrom=gtf[,'V1'], start=gtf[,'V4'], end=gtf[,'V5']), "exome.bed", quote = F, sep="\t", col.names = F, row.names = F)
RCODE

  awk '{print $1 "\t" ($2 - 1) "\t" $3}' exome.bed > exome.fixed.bed

  ${gatk_path} \
      BedToIntervalList \
      -I exome.fixed.bed \
      -O ${output_name} \
      -SD ${ref_dict}
  >>>

  output {
    File interval_list = "${output_name}"
  }

  ## runtime {
  ##   docker: docker
  ##   preemptible: preemptible_count
  ## }
}



task SplitNCigarReads {

  File input_bam
  File input_bam_index
  String base_name

  File ref_fasta
  File ref_fasta_index
  File ref_dict

  String gatk_path
  ## String docker
  Int preemptible_count

  command <<<
  ${gatk_path} \
      SplitNCigarReads \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${base_name}.bam 
  >>>

  output {
    File output_bam = "${base_name}.bam"
    File output_bam_index = "${base_name}.bai"
  }

  ## runtime {
  ##   disks: "local-disk " + sub(((size(input_bam,"GB")+1)*5 + size(ref_fasta,"GB")),"\\..*","") + " HDD"
  ##   docker: docker
  ##   memory: "4 GB"
  ##   preemptible: preemptible_count
  ## }
}

task BaseRecalibrator {

  File input_bam
  File input_bam_index
  String recal_output_file
  
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices

  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path

  ## String docker
  Int preemptible_count

  ## VM option 'PrintGCTimeStamps' does not work with java 11
  command <<<
  ${gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
   -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
               BaseRecalibrator \
               -R ${ref_fasta} \
               -I ${input_bam} \
               --use-original-qualities \
               -O ${recal_output_file} \
               -known-sites ${dbSNP_vcf} \
               -known-sites ${sep=" --known-sites " known_indels_sites_VCFs}
  >>>

  output {
    File recalibration_report = recal_output_file
  }

  ## runtime {
  ##   memory: "6 GB"
  ##   disks: "local-disk " + sub((size(input_bam,"GB")*3)+30, "\\..*", "") + " HDD"
  ##   docker: docker
  ##   preemptible: preemptible_count
  ## }
}


task ApplyBQSR {
  
  File input_bam
  File input_bam_index
  String base_name
  File recalibration_report

  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path

  ## String docker
  Int preemptible_count

  command <<<
  ${gatk_path} \
      --java-options "-XX:+PrintFlagsFinal \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
      ApplyBQSR \
      --add-output-sam-program-record \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --use-original-qualities \
      -O ${base_name}.bam \
      --bqsr-recal-file ${recalibration_report}
  >>>

  output {
    File output_bam = "${base_name}.bam"
    File output_bam_index = "${base_name}.bai"
  }

  ## runtime {
  ##   memory: "3500 MB"
  ##   disks: "local-disk " + sub((size(input_bam,"GB")*4)+30, "\\..*", "") + " HDD"
  ##   preemptible: preemptible_count
  ##   docker: docker
  ## }
}

task HaplotypeCaller {
  
  File input_bam
  File input_bam_index
  String base_name

  File interval_list

  File ref_dict
  File ref_fasta
  File ref_fasta_index
  
  File dbSNP_vcf
  File dbSNP_vcf_index

  String gatk_path
  ## String docker
  Int preemptible_count

  Int? stand_call_conf

  command <<<
  ${gatk_path} --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
	       HaplotypeCaller \
	       -R ${ref_fasta} \
	       -I ${input_bam} \
	       -L ${interval_list} \
	       -O ${base_name}.vcf.gz \
	       -dont-use-soft-clipped-bases \
	       --standard-min-confidence-threshold-for-calling ${default=20 stand_call_conf} \
	       --dbsnp ${dbSNP_vcf}
  >>>

  output {
    File output_vcf = "${base_name}.vcf.gz"
    File output_vcf_index = "${base_name}.vcf.gz.tbi"
  }

  ## runtime {
  ##   docker: docker
  ##   memory: "6.5 GB"
  ##   disks: "local-disk " + sub((size(input_bam,"GB")*2)+30, "\\..*", "") + " HDD"
  ##   preemptible: preemptible_count
  ## }
}

task VariantFiltration {

  File input_vcf
  File input_vcf_index
  String base_name

  File ref_dict
  File ref_fasta
  File ref_fasta_index

  String gatk_path
  ## String docker
  Int preemptible_count

  command <<<
  ${gatk_path} \
      VariantFiltration \
      --R ${ref_fasta} \
      --V ${input_vcf} \
      --window 35 \
      --cluster 3 \
      --filter-name "FS" \
      --filter "FS > 30.0" \
      --filter-name "QD" \
      --filter "QD < 2.0" \
      -O ${base_name}
  >>>

  output {
    File output_vcf = "${base_name}"
    File output_vcf_index = "${base_name}.tbi"
  }

  ## runtime {
  ##   docker: docker
  ##   memory: "3 GB"
  ##   disks: "local-disk " + sub((size(input_vcf,"GB")*2)+30, "\\..*", "") + " HDD"
  ##   preemptible: preemptible_count
  ## }

}


task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String output_vcf_name

    Int? disk_size = 5

    String gatk_path

    ## String docker
    Int preemptible_count

    # Using MergeVcfs instead of GatherVcfs so we can create indices
    # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
    command <<<
        ${gatk_path} --java-options "-Xms2000m"  \
            MergeVcfs \
            --INPUT ${sep=' --INPUT ' input_vcfs} \
            --OUTPUT ${output_vcf_name}
    >>>

    output {
        File output_vcf = output_vcf_name
        File output_vcf_index = "${output_vcf_name}.tbi"
    }

    ## runtime {
    ##     memory: "3 GB"
    ##     disks: "local-disk " + disk_size + " HDD"
    ##     docker: docker
    ##     preemptible: preemptible_count
    ## }
}


task ScatterIntervalList {

	File interval_list
	Int scatter_count
	String gatk_path
	## String docker
	Int preemptible_count

    command <<<
        set -e
        mkdir out
        ${gatk_path} --java-options "-Xms1g" \
            IntervalListTools \
            --SCATTER_COUNT ${scatter_count} \
            --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
            --UNIQUE true \
            --SORT true \
            --INPUT ${interval_list} \
            --OUTPUT out
	
        python3 <<CODE
        import glob, os
        # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
        intervals = sorted(glob.glob("out/*/*.interval_list"))
        for i, interval in enumerate(intervals):
          (directory, filename) = os.path.split(interval)
          newName = os.path.join(directory, str(i + 1) + filename)
          os.rename(interval, newName)
        print(len(intervals))
        if len(intervals) == 0:
          raise ValueError("Interval list produced 0 scattered interval lists. Is the gtf or input interval list empty?")
        f = open("interval_count.txt", "w+")
        f.write(str(len(intervals)))
        f.close()
        CODE
    >>>

    output {
        Array[File] out = glob("out/*/*.interval_list")
        Int interval_count = read_int("interval_count.txt")
    }

    ## runtime {
    ##     disks: "local-disk 1 HDD"
    ##     memory: "2 GB"
    ##     docker: docker
    ##     preemptible: preemptible_count
    ## }
}
