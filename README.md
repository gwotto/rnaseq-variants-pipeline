# RNA-Seq variants pipeline

WDL pipeline to call variants on RNA-seq data. Based on
https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels



## Run the pipeline


Validate pipeline, write input json file, run pipeline


```r

  java -jar ~/apps/womtool-66.jar validate rna-seq.wdl 

  java -jar ~/apps/womtool-66.jar inputs rna-seq.wdl > rna-seq.inputs.json

  java -jar ~/apps/cromwell-66.jar run rna-seq.wdl  -i rna-seq.inputs.json

```

## Release notes 

* 0.0.2 (2021-08-05)
  * options file
	included a workflow options file
  * running with cromwell 66 and java 11
    cromwell >= 66 requires java 11, which does not support some
    garbage collector logging arguments anymore, e.g.
	-XX:+PrintGCTimeStamps  
	-XX:+PrintGCDateStamps
	These were removed from the wdl script
	
