include { dictionary; samples ; build } from 'plugin/nf-htsjdk'

htsfiles_ch = Channel.fromPath(["data/*.bam","data/*.sam","data/*.bcf","data/*.vcf","data/*.vcf.gz","data/*.fai"]).
	mix(Channel.of(
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam",
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz",
	"https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam",
	"https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram"
	))



htsfiles_ch.dictionary(withLength:true).
	view()

htsfiles_ch.build().
        view()


htsfiles_ch.
	filter{!(it.toString().endsWith(".fai") || it.toString().endsWith(".fa"))}.
	samples(defaultName:"NO_SAMPLE").
	view()

