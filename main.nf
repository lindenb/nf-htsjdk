include { dictionary; samples } from 'plugin/nf-htsjdk'

htsfiles_ch = Channel.fromPath(["data/*.bam","data/*.sam","data/*.bcf","data/*.vcf","data/*.vcf.gz","data/*.fai"])


htsfiles_ch. dictionary(withLength:true).
	view()

htsfiles_ch.filter{!it.name.endsWith(".fai")}.
	samples(defaultName:"NO_SAMPLE").
	view()

