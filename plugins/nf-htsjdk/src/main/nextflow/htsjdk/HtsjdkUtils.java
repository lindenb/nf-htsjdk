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
import htsjdk.samtools.util.IOUtil;

/** TODO

 def source = targetObj instanceof List ? findSource((List)targetObj) : targetObj

    def findSource( List tuple ) {

        if( elem >= 0 )
            return tuple.get(elem)

        // find the elem-th item having Path or File type
        int pos = elem != null ? -elem : 1
        int count = 0
        for( int i=0; i<tuple.size(); i++ ) {
            def it = tuple[i]
            if(  it instanceof Path || it instanceof File ) {
                if( ++count == pos ) {
                    elem = i
                    return tuple.get(i)
                }
            }
        }

        // -- not found, if default was null just return the first item
        if( elem == null ) {
            elem = 0
            return tuple.get(0)
        }

        // -- otherwise if a `elem` value was specified but not found, raise an exception
        throw new IllegalArgumentException("Cannot find splittable file (elem=$elem)")
}

*/

public class HtsjdkUtils {
    static interface HtsSource {
        String getPath():
        String getFilename();
        boolean isRemote();
        default boolean isLocal() { return !isRemote();}
        default boolean hasSuffix(final String suff) {
            return getFilename().endsWith(suff);
            }
        default boolean hasSuffix(Collection<String> suffixes) {
            return suffixes.stream().anyMatch(S->hasSuffix(S));;
            }
        default boolean isFasta() {
            return hasSuffix(FileExtensions.FASTA_LIST);
            }
        default boolean isVcf() {
            return hasSuffix(FileExtensions.VCF_LIST);
            }
        default boolean isBam() {
            return hasSuffix(FileExtensions.BAM) || hasSuffix(FileExtensions.CRAM);
            }
        HtsSource resolveSibling(String fn);
        Collection<String> getSamples();
        Collection<String> getContigs();
        }
    static abstract class AbstractHtsSource<T> implements HtsSource {
        protected final T source;
        AbstractHtsSource(final T  source) {
            this.source = source;
            }
        @Overide String getPath() {return this.source.toString();}}
        @Override int hashCode() { return this.source.hashCode();}
        @Override public String toString() { return this.getPath();}
        }

    static class HtsPath extends AbstractHtsSource<Path> {
        HtsPath(Path source) {super(source);}
        HtsPath(File source) {super(source.toPath());}
        @Overide String getFilename() {
            return this.source.getFilename().toString();
            }
        @Override HtsSource resolveSibling(String fn) {
            return new HtsPath(source.toPath().resolve(fn));
            }
        @Override final boolean isRemote() {  return false;}
        @Override Path toPath() { return super.source;}
        }

    static class HtsUrl extends AbstractHtsSource<URL> {
        HtsUrl(URL source) {super(source);}
        @Overide String getFilename() {
            String str= this.source.getPath();
            int i=str.lastIndexOf('/');
            if(i>=0) str=str.substring(i+1);
            return str;
            }
        @Override HtsSource resolveSibling(String fn) {
            return new HtsUrl(source.toURI().resolve(fn).toURL());
            }
        @Override final boolean isRemote() {  return true;}
        @Override Path toPath() { throw new IllegalStateException("cannot get Path for URL");}
        }



    static Optional<HtsSource> toHtsFile(final Object path) {
	    if(path==null) {
	        return Optional.empty();
	        }
	    else if(path instanceof URL) {
	        return Optional.of(new HtsUrl(URL.class.cast(path));
	        }
	    else if(path instanceof String) {
	        String str = (String)path;
	        if(IOUtils.isUrl(str)) {
	            return Optional.of(new HtsUrl(new URL(path));
	            }
	        return Optional.empty();
	        }
	    else if(path instanceof File) {
	        return Optional.of(new HtsPath(File.class.cast(path));
	        }
	    else if(path instanceof Path) {
	        return Optional.of(new HtsPath(Path.class.cast(path));
	        }
	   return Optional.empty();
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

	static List<String> extractContigs(final HtsFile htsFile) throws IOException  {
		if(htsFile.hasSuffix(".gz") ) {
			try(TabixReader tbr = new TabixReader(htsFile.getPath())) {
				return new ArrayList<>(tbr.getChromosomes());
			    }
			}
	    if(htsFile.isLocal()) {
		    final Path tribble_index  = Tribble.indexPath(htsFile.toPath());
		    if(Files.exists(tribble_index)) {
    			try(InputStream in = Files.newInputStream(tribble_index)) {
	    			final Index idx = IndexFactory.loadIndex(path.toString(),in);
		    		return idx.getSequenceNames();
			    	}
			 }
		if(htsFile.isBam()) {
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

	static SAMSequenceDictionary extractDictionary(final HtsFile htsFile) throws IOException  {
		final SAMSequenceDictionary dict;
		if(htsFile.isVcf()) {
			final VCFHeader header = extractVcfHeader(path);
			if(header==null) throw new IOException("Cannot extract header from VCF file "+path);
			dict=header.getSequenceDictionary();
			if(dict==null)  throw new IOException("there is no dictionary (lines starting with '##contig=') in header of VCF file "+path);
			}
		// fai not implemented in fasta
		else if(htsFile.hasSuffix(FileExtensions.FASTA_INDEX)) {
			return dictionaryFromFai(path);
			}
		else if(htsFile.isFasta()) {
			final Path fai = path.resolveSibling(path.getFileName().toString() + FileExtensions.FASTA_INDEX);
			final Path dictPath = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(path);
			if(Files.exists(fai) && !Files.exists(dictPath)) return dictionaryFromFai(fai);
			dict = SAMSequenceDictionaryExtractor.extractDictionary(path);
			}
		else if(htsFile.isLocal()) {
			dict =	SAMSequenceDictionaryExtractor.extractDictionary(htsFile.toPath());
			}
	    if(dict==null) throw new SAMException("Cannot extract dictionary from "+htsFile);
		if(dict.isEmpty()) throw new SAMException("Cannot empty dictionary in "+htsFile);
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
