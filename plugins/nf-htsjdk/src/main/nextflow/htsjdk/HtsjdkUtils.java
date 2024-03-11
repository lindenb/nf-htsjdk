package nextflow.htsjdk;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

public class HtsjdkUtils {
	
	static boolean pathHasSuffix(Path path,final String suffix) {
		return path.getFileName().toString().endsWith(suffix);
	}
	
	static List<String> extractSamples(final Path path) throws IOException  {
		if(FileExtensions.VCF_LIST.stream().anyMatch(EXT->pathHasSuffix(path,EXT))) {
			return extractVcfHeader(path).getGenotypeSamples();
			}
		if(pathHasSuffix(path, FileExtensions.BAM) || pathHasSuffix(path, FileExtensions.CRAM)) {
			SamReaderFactory srf = SamReaderFactory.makeDefault();
			try(SamReader sr=srf.open(path)) {
				return new ArrayList<>(sr.getFileHeader().
						getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(SM->SM!=null).
						collect(Collectors.toSet())
						);
					}
				}
		throw new IOException("cannot extract samples from "+path);
		}
	
	static List<String> extractContigs(final Path path) throws IOException  {
		if(pathHasSuffix(path,".gz") ) {
			final Path tbi = path.resolveSibling(path.getFileName().toString() + FileExtensions.TABIX_INDEX);
			if(Files.exists(tbi)) {
				try(TabixReader tbr = new TabixReader(path.toString())) {
					return new ArrayList<>(tbr.getChromosomes());
					}
				}
			}
		final Path tribble_index  = Tribble.indexPath(path);
		if(Files.exists(tribble_index)) {
			try(InputStream in = Files.newInputStream(tribble_index)) {
				final Index idx = IndexFactory.loadIndex(path.toString(),in);
				return idx.getSequenceNames();
				}
			}
		if(pathHasSuffix(path, FileExtensions.BAM) || pathHasSuffix(path, FileExtensions.CRAM)) {
			SamReaderFactory srf = SamReaderFactory.makeDefault();
			try(SamReader sr=srf.open(path)) {
				if(!sr.hasIndex()) throw new IOException("BAM/CRAM file "+path+" is not indexed");
				SAMSequenceDictionary dict = sr.getFileHeader().getSequenceDictionary();
				List<String> L = new ArrayList<>();
				BAMIndex bai = sr.indexing().getIndex();
				for( int i=0;i< dict.getReferenceLength();i++) {
					BAMIndexMetaData meta= bai.getMetaData(i);
					if(meta==null ||meta.getAlignedRecordCount()==0) continue;
					L.add(dict.getSequence(i).getSequenceName());
					}
				return L;
				}
			}
		throw new IOException("cannot extract chromosomes from "+path);
		}

	
	static SAMSequenceDictionary extractDictionary(final Path path) throws IOException  {
		final SAMSequenceDictionary dict;
		if(FileExtensions.VCF_LIST.stream().anyMatch(EXT->pathHasSuffix(path,EXT))) {
			final VCFHeader header = extractVcfHeader(path);
			if(header==null) throw new IOException("Cannot extract header from VCF file "+path);
			dict=header.getSequenceDictionary();
			if(dict==null)  throw new IOException("there is no dictionary (lines starting with '##contig=') in header of VCF file "+path);
			}
		// fai not implemented in fasta
		else if(pathHasSuffix(path,FileExtensions.FASTA_INDEX)) {
			return dictionaryFromFai(path);
			}
		else if(FileExtensions.FASTA.stream().anyMatch(EXT->pathHasSuffix(path.getFileName(),EXT))) {
			final Path fai = path.resolveSibling(path.getFileName().toString() + FileExtensions.FASTA_INDEX);
			final Path dictPath = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(path);
			if(Files.exists(fai) && !Files.exists(dictPath)) return dictionaryFromFai(fai);
			dict = SAMSequenceDictionaryExtractor.extractDictionary(path);
			}
		else
			{
			dict =	SAMSequenceDictionaryExtractor.extractDictionary(path);
			}
		
		
		if(dict==null) throw new SAMException("Cannot extract dictionary from "+path);
		if(dict.isEmpty()) throw new SAMException("Cannot empty dictionary in "+path);
		return dict;
		}
	
	static SAMSequenceDictionary dictionaryFromFai(Path path) throws IOException {
		if(!pathHasSuffix(path,FileExtensions.FASTA_INDEX)) {
			throw new SAMException("not a valid extension for a fasta index file "+path);
			}
		try(BufferedReader br = Files.newBufferedReader(path,Charset.defaultCharset())) {
			return new SAMSequenceDictionary(
				br.lines().map(T->T.split("[\t]")).
					map(T->new SAMSequenceRecord(T[0],Integer.parseInt(T[1]))).
					collect(Collectors.toList())
					);
			}
		}
	
	static VCFHeader extractVcfHeader(Path path) throws IOException {
		if(FileExtensions.VCF_LIST.stream().anyMatch(EXT->path.getFileName().toString().endsWith(EXT))) {
			throw new SAMException("not a valid extension ("+String.join("|", FileExtensions.VCF_LIST)+") for a VCF file "+path);
			}
		if(path.getFileName().toString().endsWith(FileExtensions.BCF)) {
			final BCF2Codec codec = new BCF2Codec();
			  try (InputStream fis =  new BlockCompressedInputStream(Files.newInputStream(path))){
				 try(PositionalBufferedStream pbs = new PositionalBufferedStream(fis)) {
				  return (VCFHeader)codec.readHeader(pbs).getHeaderValue();
				 }
		        }
			}
		try(VCFReader r=new VCFFileReader(path,false)) {
			 return r.getHeader();
		}
	}
}
