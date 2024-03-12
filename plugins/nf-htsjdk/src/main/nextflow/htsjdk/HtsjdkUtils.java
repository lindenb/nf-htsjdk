package nextflow.htsjdk;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;



public class HtsjdkUtils {
	
	static final String KEY_ELEMENT="element";
	
    public static interface HtsSource {
        String getPath();
        String getFilename();
        default String getBasename() {
        	String str= getFilename();
        	int i = str.lastIndexOf('.');
        	return i==-1?"":str.substring(0,i);
        	}
        boolean isRemote();
        default boolean isLocal() { return !isRemote();}
        default boolean hasSuffix(final String suff) {
            return getFilename().endsWith(suff);
            }
        default boolean hasSuffix(Collection<String> suffixes) {
            return suffixes.stream().anyMatch(S->hasSuffix(S));
            }
        default boolean isFasta() {
            return hasSuffix(FileExtensions.FASTA);
            }
        default boolean isVcf() {
            return hasSuffix(FileExtensions.VCF_LIST);
            }
        default boolean isBam() {
            return hasSuffix(FileExtensions.BAM) || hasSuffix(FileExtensions.CRAM);
            }
        Path asPath();
        HtsSource resolveSibling(String fn);
        
        default SamReaderFactory samReaderFactory() {
        	return SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
        	}
        default SamReader openSamReader() throws IOException {
        	return samReaderFactory().open(asSamInputResource());
        	}
        public InputStream openInputStream() throws IOException;

        public VCFIterator openVcfIterator() throws IOException ;
        SamInputResource asSamInputResource();
        public default BufferedReader openBufferedReader() throws IOException {
        	return new BufferedReader(new InputStreamReader(mayBeGzippedInputStream(openInputStream())));
        	}
        
        
        
    	public default SAMFileHeader extractSamFileHeader() throws IOException {
    		if(!this.isBam()) {
    			throw new SAMException("not a valid extension for BAM "+getPath());
    			}
    		try(InputStream is = this.openInputStream()) {
    			try(SamReader sr=this.samReaderFactory().open(SamInputResource.of(is))) {
    				return sr.getFileHeader();
    				}
    			}
    		}
    	
    	public default VCFHeader extractVcfHeader() throws IOException {
    		if(!this.isVcf()) {
    			throw new SAMException("not a valid extension for VCF "+getPath());
    			}
    		
    		try(VCFIterator r=this.openVcfIterator()) {
    			 return r.getHeader();
    		}
    	}

    	
    	
    	default SAMSequenceDictionary dictionaryFromFai() throws IOException {
    		if(!this.hasSuffix(FileExtensions.FASTA_INDEX)) {
    			throw new SAMException("not a valid extension for a fasta index file "+getPath());
    			}
    		try(BufferedReader br = this.openBufferedReader()) {
    			return new SAMSequenceDictionary(
    				br.lines().map(T->T.split("[\t]")).
    					map(T->new SAMSequenceRecord(T[0],Integer.parseInt(T[1]))).
    					collect(Collectors.toList())
    					);
    			}
    		}

    	
    	public default SAMSequenceDictionary extractDictionary() throws IOException  {
    		SAMSequenceDictionary dict = null;
    		if(this.isLocal() && !this.hasSuffix(FileExtensions.FASTA_INDEX)) {
    			dict =	SAMSequenceDictionaryExtractor.extractDictionary(this.asPath());
    			}
    		else if(this.isVcf()) {
    			final VCFHeader header = extractVcfHeader();
    			if(header==null) throw new IOException("Cannot extract header from VCF file "+getPath());
    			dict=header.getSequenceDictionary();
    			if(dict==null)  throw new IOException("there is no dictionary (lines starting with '##"+VCFConstants.CONTIG_HEADER_KEY+"') in header of VCF file "+getPath());
    			}
    		// fai not implemented in fasta
    		else if(this.hasSuffix(FileExtensions.FASTA_INDEX)) {
    			dict= this.dictionaryFromFai();
    			}
    		else if(this.isFasta()) {
    			try {
    				final HtsSource dictPath = this.resolveSibling(this.getBasename()+ FileExtensions.DICT);
    				  try (BufferedLineReader bufferedLineReader = new BufferedLineReader( dictPath.openInputStream())) {
    					  final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
    					  final SAMFileHeader header = codec.decode(bufferedLineReader, dictPath.getPath());
    					  dict = header.getSequenceDictionary();
    				  	  }
    				  }
    			catch (final IOException e) {
    				final HtsSource fai = this.resolveSibling(this.getFilename() + FileExtensions.FASTA_INDEX);
        			dict= fai.dictionaryFromFai();
    				}
    			}
    		else if(this.isBam()) {
    			dict = extractSamFileHeader().getSequenceDictionary();
    			}
    	    if(dict==null) throw new SAMException("Cannot extract dictionary from "+getPath());
    		if(dict.isEmpty()) throw new SAMException("Cannot empty dictionary in "+getPath());
    		return dict;
    		}
    	
    	
    	public default Collection<String> extractSamples(String rgAttribute) throws IOException  {
    		if(this.isVcf()) {
    			return extractVcfHeader().getGenotypeSamples();
    			}
    		else if(this.isBam()) {
    			final String rgAtt=StringUtil.isBlank(rgAttribute)?SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG:rgAttribute;
				return extractSamFileHeader().
					getReadGroups().
					stream().
					map(RG->RG.getAttribute(rgAtt)).
					filter(SM->!StringUtil.isBlank(SM)).
					collect(Collectors.toSet()
					);
    			}
    		throw new IOException("cannot extract samples from "+getPath());
    		}
    	
    	public default List<String> extractContigs() throws IOException  {
    		if(this.hasSuffix(".gz") ) {
    			try(TabixReader tbr = new TabixReader(this.getPath())) {
    				return new ArrayList<>(tbr.getChromosomes());
    			    }
    			}
    	    if(this.isLocal()) {
    		    final Path tribble_index  = Tribble.indexPath(this.asPath());
    		    if(Files.exists(tribble_index)) {
        			try(InputStream in = mayBeGzippedInputStream(Files.newInputStream(tribble_index))) {
    	    			final Index idx = IndexFactory.loadIndex(this.getPath(),in);
    		    		return idx.getSequenceNames();
    			    	}
    		    	}
    			 }
    		if(this.isBam()) {
    			try(SamReader sr= this.openSamReader()) {
    				if(!sr.hasIndex()) throw new IOException("BAM/CRAM file "+getPath()+" is not indexed");
    				final SAMSequenceDictionary dict = sr.getFileHeader().getSequenceDictionary();
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
    		throw new IOException("cannot extract chromosomes from "+getPath());
    		}


        }
    
    public static abstract class AbstractHtsSource<T> implements HtsSource {
        protected final T source;
        AbstractHtsSource(final T  source) {
            this.source = source;
            }
        @Override public String getPath() {return this.source.toString();}
        @Override public int hashCode() { return this.source.hashCode();}
        @Override public String toString() { return this.getPath();}
       
        }

    public static class HtsPath extends AbstractHtsSource<Path> {
        HtsPath(Path source) {super(source);}
        HtsPath(File source) {super(source.toPath());}
        @Override public String getFilename() {
            return super.source.getFileName().toString();
            }
        @Override public HtsSource resolveSibling(String fn) {
            return new HtsPath(super.source.resolve(fn));
            }
        @Override public final boolean isRemote() { return false;}
        @Override public  Path asPath() { return super.source;}
        @Override public SamInputResource asSamInputResource() {
        	return SamInputResource.of(asPath());
        	}
        @Override
        public VCFIterator openVcfIterator() throws IOException {
        	return new VCFIteratorBuilder().open(asPath());
        	}
        public InputStream openInputStream() throws IOException {
        	return Files.newInputStream(asPath());
        	}
        }

    public  static class HtsUrl extends AbstractHtsSource<URL> {
        HtsUrl(URL source) {super(source);}
        @Override public String getFilename() {
            String str= this.source.getPath();
            int i=str.lastIndexOf('/');
            if(i>=0) str=str.substring(i+1);
            return str;
            }
        @Override public HtsSource resolveSibling(String fn) {
            try {
				return new HtsUrl(source.toURI().resolve(fn).toURL());
			} catch (MalformedURLException | URISyntaxException e) {
				throw new IllegalArgumentException("cannot resolve "+fn+" to URL", e);
			}
            }
        @Override public final boolean isRemote() { return true;}
        @Override public Path asPath() { throw new IllegalStateException("cannot get Path for URL");}
        
        @Override public SamInputResource asSamInputResource() {
        	return SamInputResource.of(super.source);
        	}
        @Override
        public VCFIterator openVcfIterator() throws IOException {
        	return new VCFIteratorBuilder().open(getPath());
        	}
        @Override
        public InputStream openInputStream() throws IOException {
        	return ParsingUtils.openInputStream(getPath());
        	}

    	}



    static Optional<HtsSource> toHtsSource(final Object path) {
	    if(path==null) {
	        return Optional.empty();
	        }
	    else if(path instanceof URL) {
	        return Optional.of(new HtsUrl(URL.class.cast(path)));
	        }
	    else if(path instanceof String) {
	        if(IOUtil.isUrl(String.class.cast(path))) {
	            try {
					return toHtsSource(new URL(String.class.cast(path)));
				} catch (MalformedURLException e) {
					throw new IllegalArgumentException("cannot convert "+path+" to URL", e);
				}
	            }
	        return Optional.empty();
	        }
	    else if(path instanceof File) {
	        return Optional.of(new HtsPath(File.class.cast(path)));
	        }
	    else if(path instanceof Path) {
	        return Optional.of(new HtsPath(Path.class.cast(path)));
	        }
	   return Optional.empty();
       }

    static HtsSource findHtsSource(final Object o, Object key, final Collection<String> extensions)  {
    	if(o==null) {
    		throw new IllegalArgumentException("Object cannot be null");
    		}
    	else if(o instanceof Map) {
    		final Map<?,?> hash = Map.class.cast(o);
    		if(key!=null) {
    			if(!hash.containsKey(key)) {
    				throw new IllegalArgumentException("Cannot find object key="+key+" in map : "+o);
    				}
    			Object value = hash.get(key);
    			HtsSource source = toHtsSource(value).orElse(null);
    			if(source==null) {
    				throw new IllegalArgumentException("Cannot convert "+value+" to HTS file");
    				}
    			if(extensions!=null && !source.hasSuffix(extensions)) {
    				throw new IllegalArgumentException("HTS file "+source+" was found but has no valid suffixes : "+String.join(",", extensions));
    				}
    			return source;
    			}
    		else
    			{
    			List<HtsSource> L = hash.values().
    					stream().
    					map(F->HtsjdkUtils.toHtsSource(F).orElse(null)).
    					filter(F->F!=null && F.hasSuffix(extensions)).
    					collect(Collectors.toList());
    			if(L.isEmpty()){
    				throw new IllegalArgumentException("Cannot find a valid source in "+o+" with suffixes : "+String.join(",", extensions));
    				}
    			// return first
    			return L.get(0);
    			}
    		}
    	else if(o instanceof List) {
    		final List<?> list = List.class.cast(o);
    		int index=-1;
    		if(key!=null) {
    			if(key instanceof Integer) {
    				index = Integer.class.cast(key).intValue();
    				}
    			else
    				{
    				throw new IllegalArgumentException("Index for list provided  "+key+" but it's a "+key.getClass()+" not an integer");
    				}
    			if(index< 0 || index >= list.size()) {
    				throw new IndexOutOfBoundsException("Index provided idx="+index+" but list has size:"+list.size());
    				}
    			Object value = list.get(index);
    			HtsSource source = toHtsSource(value).orElse(null);
    			if(source==null) {
    				throw new IllegalArgumentException("Cannot convert "+value+" to HTS file");
    				}
    			if(extensions!=null && !source.hasSuffix(extensions)) {
    				throw new IllegalArgumentException("HTS file "+source+" was found but has no valid suffixes : "+String.join(",", extensions));
    				}
    			return source;
    			}
    		else
    			{
    			List<HtsSource> L = list.stream().map(F->HtsjdkUtils.toHtsSource(F).orElse(null)).
					filter(F->F!=null && (extensions==null || F.hasSuffix(extensions))).
					collect(Collectors.toList());
    			if(L.isEmpty()){
    				throw new IllegalArgumentException("Cannot find a valid source in "+o+" with suffixes : "+String.join(",", extensions));
    				}
    			// return first
    			return L.get(0);
    			}
    		}
    	else
    		{
    		HtsSource hts = toHtsSource(o).orElse(null);
    		if(hts==null) throw new IllegalArgumentException("Cannot convert "+o+" to a HTS file");
    		if(!hts.hasSuffix(extensions)) throw new IllegalArgumentException("got hts file"+o+" but with wrong suffixes : "+String.join(",", extensions));
			return hts;
    		}
    	}
	
	
	private static InputStream mayBeGzippedInputStream(InputStream in) throws IOException {
	    // wrap the input stream into a BufferedInputStream to reset/read a BCFHeader or a GZIP
	    // buffer must be large enough to contain the BCF header and/or GZIP signature
	    BufferedInputStream  bufferedinput = new BufferedInputStream(in, Math.max(BCF2Codec.SIZEOF_BCF_HEADER, IOUtil.GZIP_HEADER_READ_LENGTH));
	    // test for gzipped inputstream
	    if(IOUtil.isGZIPInputStream(bufferedinput)) {
	        // this is a gzipped input stream, wrap it into GZIPInputStream
	        // and re-wrap it into BufferedInputStream so we can test for the BCF header
	        bufferedinput = new BufferedInputStream(new GZIPInputStream(bufferedinput), BCF2Codec.SIZEOF_BCF_HEADER);
	    	}
	    return bufferedinput;
		}
}
