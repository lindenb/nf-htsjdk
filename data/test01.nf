include { dictionary ; build } from 'plugin/nf-htsjdk'

Channel.fromPath(
		"data/rotavirus_rf.fa",
		).
	dictionary().
	view()

ch1 = Channel.of(
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam",
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz",
	"https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam",
	"https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram"
	)


ch1.dictionary().view()
ch1.build().view()
