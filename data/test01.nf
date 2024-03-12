include { dictionary ; build ; samples ; mappedContigs } from 'plugin/nf-htsjdk'

ch0 = Channel.fromPath([
		"data/rotavirus_rf.fa",
		"data/rotavirus_rf.bed.gz",
		"data/rotavirus_rf.vcf.gz",
		//"data/rotavirus_rf.bcf",
		"data/S1.rota.bam",
		"data/S1.rota.cram"
		])


ch1 = Channel.of(
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam",
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz",
	"https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam",
	"https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram"
	)


ch1.dictionary().
	filter{it[2].matches("((chr)?1|ref|RF01)")}.
	view{"DICT:$it"}

ch1.build().view{"BUILD:$it"}


// filter{!it.name.endsWith(".bed.gz")}.
// 	dictionary().

ch0.filter{!it.name.endsWith(".bed.gz")}.
	dictionary().
	filter{it[2].matches("((chr)?1|ref|RF01)")}.
	view{"LOCAL-DICT: ${it}"}


ch0.filter{!it.name.endsWith(".bed.gz")}.
	filter{!it.name.endsWith(".fa")}.
	samples().
	view{"SAMPLES: ${it}"}

ch0.	filter{!it.name.endsWith(".fa")}.
	mappedContigs().
	view{"MAPPED: ${it}"}
