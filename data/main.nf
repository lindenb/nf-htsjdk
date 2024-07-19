include { dictionary; samples ; build } from 'plugin/nf-htsjdk'

/** local files */
local_files_ch = Channel.fromPath([
	"data/*.bam",
	"data/*.sam",
	"data/*.bcf",
	"data/*.vcf",
	"data/*.vcf.gz",
	"data/*.fai",
	"data/*.interval_list",
	"data/*.interval_list.gz"
	])

/** remote files */
remote_files_ch = Channel.of(
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam",
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz",
	"https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam",
	"https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram"
	)

htsfiles_ch= local_files_ch.mix(remote_files_ch)



htsfiles_ch.
	map{[build(it),it]}.
	map{[it[0]==null?"NO_BUILD":it[0].getId(),it[1]]}.
	view{"BUILD: ${it}"}


htsfiles_ch.                
	map{[dictionary(it),it]}.
        flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}.
	filter{it[0].matches("((chr)?[12]|RF.*)")}.
        view{"DICTIONARY: ${it}"}

htsfiles_ch.
	filter{!(it.toString().endsWith(".fai") || it.toString().endsWith(".fa"))}.
        map{[samples(it,["defaultName":"NO_SAMPLE"]),it]}.
        flatMap{row->row[0].collect{sn->[sn,row[1]]}}.
	view{"SAMPLE: ${it}"}

