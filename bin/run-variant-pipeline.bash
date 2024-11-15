## needs gatk and Rscript

module load java/11.0.2
module load python/3.6.1
module load gatk/4.1.4.1
module load R/4.0.2


## firts argument is options file
## second argument is input directory, where bam files are located
## second argument is output directory

opt_file=$1

in_dir=$2

out_dir=$3

log_dir=${in_dir}/logs

for file_path in $(ls ${in_dir}/*.bam)
      do
      in_file=$(basename $file_path)
      file_base=${in_file%%.*}
      
      cat > ${file_base}_temp.json <<EOF      
{
      "RNAseq.inputBam": "${file_path}",
      "RNAseq.inputBamIndex": "${file_path}.bai",
      "RNAseq.annotationsGTF": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/gencode.v38.annotation.gtf",
      "RNAseq.refFasta": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
      "RNAseq.refFastaIndex": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai",
      "RNAseq.refDict": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict",
      "RNAseq.dbSnpVcf": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/dbsnp_138.hg38.vcf.gz",
      "RNAseq.dbSnpVcfIndex": "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/dbsnp_138.hg38.vcf.gz.tbi",
      "RNAseq.knownVcfs": [
			  "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
			  "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/Homo_sapiens_assembly38.known_indels.vcf.gz"
    			  ],
      "RNAseq.knownVcfsIndices": [
      			   "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
			   "/SAN/CTO/reference/GRCh38/GRCh38_no_alt_analysis_set/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
    			   ],
      "RNAseq.java_opt": "String",
      "RNAseq.gatk_path": "/share/apps/cto/packages/GATK/4.1.4.1/gatk",
      "##RNAseq.interval_list": "/home/gwo/project/rnaseq-variants-pipeline/data/reference/gencode.v38.exon.ctnnb1.interval_list",
      "##RNAseq.minConfidenceForVariantCalling": "(optional) Int?",
      "##RNAseq.haplotypeScatterCount": "(optional) Int?",
      "##RNAseq.preemptible_tries": "(optional) Int?"
}
EOF


      qsub -V -cwd -l tmem=32G,h_rt=240:00:00 -o $log_dir -e $log_dir -N ${file_base}_variants <<EOF
      module list
      java -jar /share/apps/cto/packages/wdl/66/cromwell-66.jar run /share/apps/cto/rnaseq-variants-pipeline/bin/variant-pipeline.wdl -i ${file_base}_temp.json -o $opt_file
EOF

done
