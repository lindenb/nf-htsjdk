# nf-htsjdk plugin 
 
This project contains a simple Nextflow plugin called `nf-htsjdk` which provides operators to handle HTS files (BAM, VCF, FASTA) in nextflow using HTSJDK ,
the Java API for high-throughput sequencing data (HTS) formats ( https://github.com/samtools/htsjdk/ ).

## Unit testing 

To run your unit tests, run the following command in the project root directory (ie. where the file `settings.gradle` is located):

```bash
./gradlew check
```

## Usage

first include the new operators in the header of your nextflow script:

```nextflow
include { dictionary; samples ; build } from 'plugin/nf-htsjdk'
```

Inputs are common local or remote HTS files:

```nextflow

local_files_ch = Channel.fromPath(["data/*.bam","data/*.sam","data/*.bcf","data/*.vcf","data/*.vcf.gz","data/*.fai"])

remote_files_ch = Channel.of(
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa",
	"http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam",
	"https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz",
	"https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam",
	"https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram"
	)

htsfiles_ch= local_files_ch.mix(remote_files_ch)
```


### dictionary

The `dictionary` operator extracts the SAM Sequence dictionary from BAM/CRAM/SAM/VCF/BCF/FAI/DICT file and generate a new item for each chromosome in the dictionary.

```nextflow
htsfiles_ch.dictionary(withLength:true).
	view{"DICTIONARY: ${it}"}

```

output

```
DICTIONARY: [10, 135534747, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
DICTIONARY: [9, 141213431, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
DICTIONARY: [chr10, 133797422, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr10, 133797422, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
DICTIONARY: [chr10, 133797422, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
DICTIONARY: [chr10_GL383545v1_alt, 179254, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr11, 135086622, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr11, 135086622, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
DICTIONARY: [chr11, 135086622, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
DICTIONARY: [chr11_KI270721v1_random, 100316, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr11_KI270721v1_random, 100316, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
DICTIONARY: [chr11_KI270902v1_alt, 106711, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr11_KI270903v1_alt, 214625, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
(...)
```
## samples

The `samples` operator  extract the samples' name from BAM/CRAM/SAM/VCF/BCF files and generate a new item for each sample.

```nextflow
htsfiles_ch.
	filter{!(it.toString().endsWith(".fai") || it.toString().endsWith(".fa"))}.
	samples(defaultName:"NO_SAMPLE").
	view{"SAMPLE: ${it}"}
```


output:

```
SAMPLE: [HG00101, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
SAMPLE: [NO_SAMPLE, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
SAMPLE: [NO_SAMPLE, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
SAMPLE: [S1, https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram]
SAMPLE: [S1, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
SAMPLE: [NO_SAMPLE, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.noRG.sam]
SAMPLE: [S1, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
SAMPLE: [S2, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
SAMPLE: [S3, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
SAMPLE: [S4, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
SAMPLE: [S5, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
SAMPLE: [NO_SAMPLE, /home/lindenb/src/nf-htsjdk/data/vcf_no_genotype.vcf]
SAMPLE: [S1, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
SAMPLE: [S2, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
SAMPLE: [S3, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
SAMPLE: [S4, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
SAMPLE: [S5, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
```

# build

The `build` operator adds a new column with the name of the build found in the dictionary (buils are currently defined in `./plugins/nf-htsjdk/build/resources/main/META-INF/builds.xml`

```nextflow
htsfiles_ch.build().
        view{"BUILD: ${it}"}
```

output

```
BUILD: [GRCh38, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
BUILD: [GRCh37, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
BUILD: [GRCh38, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
BUILD: [GRCh38, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
BUILD: [samtools, https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.noRG.sam]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/vcf_no_genotype.vcf]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.fa.fai]
```



