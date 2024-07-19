# nf-htsjdk plugin 
 
This project contains a simple Nextflow plugin called `nf-htsjdk` which provides **functions** to handle HTS files (BAM, VCF, FASTA, GATK interval_list) in nextflow using HTSJDK , the Java API for high-throughput sequencing data (HTS) formats ( https://github.com/samtools/htsjdk/ ).

## Unit testing 

To run your unit tests, run the following command in the project root directory (ie. where the file `settings.gradle` is located):

```bash
./gradlew check
```

## Usage

first include the new **functions** in the header of your nextflow script:

```nextflow
include { dictionary; samples ; build } from 'plugin/nf-htsjdk'
```

Inputs are common local or remote HTS files:

```nextflow

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

The `dictionary` function extracts the SAM Sequence dictionary from BAM/CRAM/SAM/VCF/BCF/FAI/DICT file returns a `htsjdk.samtools.SAMSequenceDictionary` [https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/SAMSequenceDictionary.java](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/SAMSequenceDictionary.java)

```nextflow
htsfiles_ch.                
    map{[dictionary(it),it]}.
    flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}.
    filter{it[0].matches("((chr)?[12]|RF.*)")}.
    view{"DICTIONARY: ${it}"}
```

output

```
DICTIONARY: [chr1, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [chr2, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
DICTIONARY: [1, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
DICTIONARY: [2, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
DICTIONARY: [chr1, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
DICTIONARY: [chr2, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
DICTIONARY: [chr1, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
DICTIONARY: [chr2, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
DICTIONARY: [RF01, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
DICTIONARY: [RF02, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
DICTIONARY: [RF03, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
DICTIONARY: [RF04, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
DICTIONARY: [RF05, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
DICTIONARY: [RF06, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
(...)
```
## samples

The `samples` function  extract the samples' name from BAM/CRAM/SAM/VCF/BCF files as a collection of Strings.

| parameter | type | description |
|-----------|------|-------------|
| defaultName | string | if no sample is available, use this mock name. Otherwise there will be no output |


```nextflow
htsfiles_ch.
	filter{!(it.toString().endsWith(".fai") || it.toString().endsWith(".fa"))}.
    map{[samples(it,["defaultName":"NO_SAMPLE"]),it]}.
    flatMap{row->row[0].collect{sn->[sn,row[1]]}}.
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
SAMPLE: [NO_SAMPLE, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.interval_list]
SAMPLE: [NO_SAMPLE, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.interval_list.gz]
```

# build

The `build` function returns an object of type `Build` that is associated to a dictionary, or null if the dictionary cannot be found. (buils are currently defined in `./plugins/nf-htsjdk/build/resources/main/META-INF/builds.xml`

The Build interface is: 
```
public interface Build {
    /** return the build identifier */
    public String getId();
    }
```



```nextflow

htsfiles_ch.
	map{[build(it),it]}.
	map{[it[0]==null?"NO_BUILD":it[0].getId(),it[1]]}.
	view{"BUILD: ${it}"}

```

output

```
BUILD: [GRCh38, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa]
BUILD: [GRCh37, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/HG00101/exome_alignment/HG00101.mapped.illumina.mosaik.GBR.exome.20110411.bam]
BUILD: [GRCh38, https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.chr1.vcf.bgz]
BUILD: [GRCh38, https://www.encodeproject.org/files/ENCFF656TNN/@@download/ENCFF656TNN.bam]
BUILD: [NO_BUILD, https://github.com/lindenb/jvarkit/raw/master/src/test/resources/toy.cram]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/S1.rota.bam]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.noRG.sam]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.bcf]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/vcf_no_genotype.vcf]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.vcf.gz]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.fa.fai]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.interval_list]
BUILD: [rotavirus, /home/lindenb/src/nf-htsjdk/data/rotavirus_rf.interval_list.gz]
```

## Author

Pierre Lindenbaum PhD
