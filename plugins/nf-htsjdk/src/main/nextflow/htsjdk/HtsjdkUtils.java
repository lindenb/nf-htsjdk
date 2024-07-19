/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package nextflow.htsjdk;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;


import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCF2Type;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;



public class HtsjdkUtils {
    private static final String BUILD_RESOURCE = "/META-INF/builds.xml";
	
	/** search for compatible contig in a SAMSequenceDictionary */
	private static interface ContigMatcher {
		abstract public  boolean test(boolean resolveChromosome,final SAMSequenceDictionary dict);
		}
	
	/** matcher for contig name and contig length */
	private static class ContigLengthMatcher implements ContigMatcher {
		private final String contigName;
		private final int contigLen;
		ContigLengthMatcher(String contigName,int contigLen) {
			this.contigName = contigName;
			this.contigLen = contigLen;
			}
		@Override
		public boolean test(boolean resolveChromosome,final SAMSequenceDictionary dict) {
			if(!resolveChromosome) {
				final SAMSequenceRecord ssr= dict.getSequence(this.contigName);
				if(ssr==null) return false;
				return ssr.getSequenceLength()==this.contigLen;
				}
			
			return dict.getSequences().
					stream().
					anyMatch(SRR->contigLen==(SRR.getSequenceLength()) && contigsMatch(contigName,SRR.getContig()));
			
			}
		private boolean contigsMatch(String s1,String s2) {
			return simpleChromName(s1).equals(simpleChromName(s2));
			}
		
		private String simpleChromName(String s) {
			s=s.toLowerCase();
			if(s.startsWith("chr")) s=s.substring(3);
			if(s.equals("m")) return "mt";
			return s;
			}
		}
	private static class ContigMD5Matcher implements ContigMatcher {
		private final String md5;
		ContigMD5Matcher(String md5) {
			this.md5 = md5;
			}
		@Override
		public boolean test(boolean resolveChromosome,final SAMSequenceDictionary dict) {
			return dict.getSequences().
					stream().
					anyMatch(SRR->md5.equals(SRR.getMd5()));
			}
		}

	
	
	public static interface Build {
	    public boolean match(boolean resolveChromosome,SAMSequenceDictionary dict);
	    public String getOrganism();
	    public default String getVersion() { return ".";}
	    public String getId();
	    /** alias for getId */
	    public default String getName() {
	    	return this.getId();
	    	}
	    };
	
	private static class BuildImpl implements Build  {
		String id=".";
		String version=".";
		String organism=".";
		final List<ContigMatcher> predicates = new ArrayList<>();
		@Override  
		public boolean match(boolean resolveChromosome,final SAMSequenceDictionary dict) {
	    	return predicates.stream().allMatch(P->P.test(resolveChromosome,dict));
	    	}
	    @Override public String getOrganism() {return organism;}
	    @Override public String getVersion() {return version;}
	    @Override public String getId() {return id;}
	    };
		
	    
	private static List<Build> BUILDS = null;
	
	/** warpper around HTS files (bam, cram, vcf, fasta, etc...) */
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
            return hasSuffix(FileExtensions.VCF_LIST) || hasSuffix(".vcf.bgz") /* gnomad */;
            }
        default boolean isBamCramSam() {
            return isBam() || isCram() || isSam();
            }
        default boolean isBam() {
            return hasSuffix(FileExtensions.BAM);
            }
        default boolean isCram() {
            return hasSuffix(FileExtensions.CRAM);
            }
        default boolean isSam() {
            return hasSuffix(FileExtensions.SAM);
            }
        default boolean isFai() {
            return hasSuffix(FileExtensions.FASTA_INDEX);
            }
        default boolean isDict() {
            return hasSuffix(FileExtensions.DICT);
            }
        default boolean isIntervalList() {
        	 return hasSuffix(FileExtensions.INTERVAL_LIST) || hasSuffix(FileExtensions.COMPRESSED_INTERVAL_LIST);
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
    		if(isIntervalList()) {
    			try(BufferedReader r=this.openBufferedReader()) {
    				return IntervalList.fromReader(r).getHeader();
    				}
    			}
    		if(!isBamCramSam()) {
    			throw new SAMException("not a valid extension for BAM "+getPath());
    			}
    		try(InputStream is = this.openInputStream()) {
    			try(SamReader sr=this.samReaderFactory().open(SamInputResource.of(is))) {
    				return sr.getFileHeader();
    				}
    			}
    		}
    	public default List<SAMReadGroupRecord> extractReadGroups() throws IOException {
    		final SAMFileHeader header= extractSamFileHeader();
    		return header==null?Collections.emptyList():header.getReadGroups();
    		}
    	
    	public default VCFHeader extractVcfHeader() throws IOException {
    		if(!this.isVcf()) {
    			throw new SAMException("not a valid extension for VCF "+getPath());
    			}
    		
    	
    		try(InputStream in=this.openInputStream()) {
    			 return decodeVCFHeader(in);
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
    		if(this.isVcf()) {
    			final VCFHeader header = extractVcfHeader();
    			if(header==null) throw new IOException("Cannot extract header from VCF file "+getPath());
    			dict=header.getSequenceDictionary();
    			if(dict==null)  throw new IOException("there is no dictionary (lines starting with '##"+VCFConstants.CONTIG_HEADER_KEY+"') in header of VCF file "+getPath());
    			}
    		else if(this.isLocal() && !this.hasSuffix(FileExtensions.FASTA_INDEX)) {
    			dict =	SAMSequenceDictionaryExtractor.extractDictionary(this.asPath());
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
    		else if(isBamCramSam()) {
    			dict = extractSamFileHeader().getSequenceDictionary();
    			}
    	    if(dict==null) throw new SAMException("Cannot extract dictionary from "+getPath());
    		if(dict.isEmpty()) throw new SAMException("Cannot empty dictionary in "+getPath());
    		return dict;
    		}

    	public default Collection<String> extractSamples() throws IOException  {
    	    return extractSamples(SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG);
    	    }
    	
    	
    	public default Collection<String> extractSamples(final String rgAttribute) throws IOException  {
    		if(this.isVcf()) {
    			return extractVcfHeader().getGenotypeSamples();
    			}
    		else if(this.isBamCramSam() || this.isIntervalList()) {
    			final String rgAtt=StringUtil.isBlank(rgAttribute)?SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG:rgAttribute;
				return extractReadGroups().
					stream().
					map(RG->RG.getAttribute(rgAtt)).
					filter(SM->!StringUtil.isBlank(SM)).
					collect(Collectors.toSet()
					);
    			}
    		throw new IOException("cannot extract samples from "+getPath());
    		}
    	
    	public default List<String> extractMappedContigs() throws IOException  {
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
    		if(this.isBam() || this.isCram()) {
    			try(SamReader sr= this.openSamReader()) {
    				if(!sr.hasIndex()) throw new IOException("BAM/CRAM file "+getPath()+" is not indexed");
    				final SAMSequenceDictionary dict = sr.getFileHeader().getSequenceDictionary();
    				final List<String> L = new ArrayList<>();
    				final BAMIndex bai = sr.indexing().getIndex();
    				for( int i=0;i< dict.size();i++) {
    				    if(hasSuffix(".cram")) {
    				        try(CloseableIterator<SAMRecord> iter= sr.queryAlignmentStart(dict.getSequence(i).getSequenceName(),1)) {
    				            if(iter.hasNext())  L.add(dict.getSequence(i).getSequenceName());
    				            }
    				        }
    				    else {
    					    final BAMIndexMetaData meta= bai.getMetaData(i);
        					if(meta==null ||meta.getAlignedRecordCount()==0) {
        					        continue;
    	    				        }
    		    			L.add(dict.getSequence(i).getSequenceName());
    		    			}
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

    private static Build parseBuild(XMLEventReader xr) throws XMLStreamException {
    	final BuildImpl b = new BuildImpl();
    	while(xr.hasNext()) {
			XMLEvent evt = xr.nextEvent();
			if(evt.isStartElement() ) {
				StartElement e = evt.asStartElement();
				String name = e.getName().getLocalPart();
				if(name.equals("version")) {
					b.version = xr.getElementText().trim();
					}
				else if(name.equals("id")) {
					b.id = xr.getElementText().trim();
					}
				else if(name.equals("organism")) {
					b.organism = xr.getElementText().trim();
					}
				else if(name.equals("md5")) {
					final String md5 =  xr.getElementText().trim();
					b.predicates.add(new ContigMD5Matcher(md5));
					}
				else if(name.equals("contig")) {
					Attribute a =e.getAttributeByName(new QName("name"));
					if(a==null) throw new XMLStreamException("@name missing ",e.getLocation());
					final String chrom = a.getValue();
					a =e.getAttributeByName(new QName("length"));
					if(a==null) throw new XMLStreamException("@length missing ",e.getLocation());
					final int len = Integer.parseInt(a.getValue());
					b.predicates.add(new ContigLengthMatcher(chrom, len));
					}
				else
					{
					throw new XMLStreamException("undefined tag <"+name+">",e.getLocation());
					}
				}
			else if(evt.isEndElement()) {
				final EndElement e = evt.asEndElement();
				if(e.getName().getLocalPart().equals("build")) {
					if(b.predicates.isEmpty()) throw new XMLStreamException("no predicate for build",e.getLocation());
					return b;
					}
				}
			}
    	throw new IllegalStateException("error in XML");
    	}
    
    private static synchronized List<Build> getDefaultBuilds() {
    	if(BUILDS==null) {
    		synchronized (HtsjdkUtils.class) {
    			if(BUILDS==null) {
    				BUILDS= new ArrayList<>();
    				try(InputStream in = HtsjdkUtils.class.getResourceAsStream(BUILD_RESOURCE)) {
    					
    					if(in!=null) {
    						XMLInputFactory xif = XMLInputFactory.newFactory();
    						XMLEventReader xr = xif.createXMLEventReader(in);
    						while(xr.hasNext()) {
    							XMLEvent evt = xr.nextEvent();
    							if(evt.isStartElement() && evt.asStartElement().getName().getLocalPart().equals("build")) {
    								BUILDS.add(parseBuild(xr));
    								}
    							}
    						}
    					}
    				catch(final IOException|XMLStreamException err) {
    					err.printStackTrace();
	    				}
	    			}
				}
	    	}
    	return BUILDS;
    	}
    

   

   @SuppressWarnings("rawtypes")
   static List<Build> decodeBuilds(final Object o1) {
	if(o1==null) {
		return getDefaultBuilds();
		}
	if(!(o1 instanceof List)) throw new IllegalArgumentException("in htsjdk config. Expected htsjdk.builds as List but got a "+o1.getClass());
	List L1 = List.class.cast(o1);
	if(L1.isEmpty()) {
		return Collections.emptyList();
		}
	final List<Build> builds = new ArrayList<>(L1.size());
	
	for(int i=0;i< L1.size();i++) {
		final Object o2 = L1.get(i);
		if(o2==null)  throw new IllegalArgumentException("config htsjdk : builds["+i+"] is null");
		if(!(o2 instanceof Map)) throw new IllegalArgumentException("config htsjdk : builds["+i+"] expected a Map but got "+o2.getClass());
		final Map hash2 = Map.class.cast(o2);
		if(!hash2.containsKey("name"))  throw new IllegalArgumentException("config htsjdk : builds["+i+"].name is undefined");
		Object buildName = hash2.get("name");
		if(!(buildName instanceof String))  throw new IllegalArgumentException("config htsjdk : builds["+i+"].name is not a string "+buildName.getClass());
		final BuildImpl buildImpl = new BuildImpl();
		buildImpl.id = String.valueOf(buildName);
		if(!hash2.containsKey("chromosomes"))  throw new IllegalArgumentException("config htsjdk : builds["+i+"].chromosomes is undefined");
		final Object o3 = hash2.get("chromosomes");
		if(!(o3 instanceof List)) throw new IllegalArgumentException("in htsjdk config. Expected builds["+i+"].chromosomes as List but got a "+o3.getClass());
		List L3 = List.class.cast(o3);
		if(L3.isEmpty()) {
			throw new IllegalArgumentException("in htsjdk config. Expected builds["+i+"].chromosomes is empty");
			}
		for(int j=0;j< L3.size();j++) {
			final Object o4 = L3.get(j);
			if(o4==null)  throw new IllegalArgumentException("config htsjdk : builds["+i+"].chromosomes["+j+"] is null");
			if(!(o4 instanceof Map)) throw new IllegalArgumentException("config htsjdk : builds["+i+"].chromosomes["+j+"] expected a Map but got "+o2.getClass());
			final Map hash4 = Map.class.cast(o4);
			if(hash4.containsKey("md5"))  {
				final String md5 =  String.class.cast(hash4.get("md5"));
				buildImpl.predicates.add(new ContigMD5Matcher(md5));
				continue;
				}

			
			
			if(!hash4.containsKey("name"))  throw new IllegalArgumentException("config htsjdk : builds["+i+"].chromosomes["+j+"].name is undefined");
			final String contigName = String.class.cast(hash4.get("name"));
			if(hash4.containsKey("length")) {
				final int length = Integer.class.cast(hash4.get("length"));
				buildImpl.predicates.add(new ContigLengthMatcher(contigName, length));
				}
			else
				{
				throw new IllegalArgumentException("config htsjdk : builds["+i+"].chromosomes["+j+"].length is undefined");
				}
			}
		}
	return builds;
   	}
   
   
    /** bild an HTS source from an object */
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
	        return Optional.of(new HtsPath(Paths.get(String.class.cast(path))));
	        }
	    else if(path instanceof File) {
	        return Optional.of(new HtsPath(File.class.cast(path)));
	        }
	    else if(path instanceof Path) {
	        return Optional.of(new HtsPath(Path.class.cast(path)));
	        }
	   return Optional.empty();
       }


    /**
     * convert Object 'source ' to 'HtsSource'
     */
    static HtsSource findHtsSource(final Object source, final Predicate<HtsSource> acceptHtsSource)  {
    	if(source==null) {
    		throw new IllegalArgumentException("source cannot be null");
    		}
    	final HtsSource htsSource = toHtsSource(source).orElse(null);
		if(htsSource==null) {
			throw new IllegalArgumentException("Cannot convert "+source+" to HTS file");
			}
		if(!acceptHtsSource.test(htsSource)) {
			throw new IllegalArgumentException("HTS file "+source+" was found but has no valid suffixes");
			}
    	
		return htsSource;
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

	
	

	/** extract VCFHeader from stream. Once decoded, the stream is invalid (uncompressed, bytes skipped, etc..) */
	private static VCFHeader decodeVCFHeader(InputStream in0) throws IOException {
		BufferedInputStream in = new BufferedInputStream(mayBeGzippedInputStream(in0),BCFVersion.MAGIC_HEADER_START.length);
        in.mark(2);
		final VCFCodec headerParser = new VCFCodec();

		// read BCF version
        final byte[] magicBytes = new byte[BCFVersion.MAGIC_HEADER_START.length];
        in.read(magicBytes);
        
        if ( !Arrays.equals(magicBytes, BCFVersion.MAGIC_HEADER_START) ) {
        	in.reset();
        	/* looks like a regular VCF */
 	       try( final PositionalBufferedStream bps = new PositionalBufferedStream(in)) {
		        final LineIterator lineIterator = new LineIteratorImpl(new SynchronousLineReader(bps));
		        return (VCFHeader)headerParser.readActualHeader(lineIterator);
		       	}
        	}
        else
	        {
	        // we're a BCF file
	        final int majorByte = in.read();
	        final int minorByte = in.read();
	        new BCFVersion( majorByte, minorByte );
	    	
	        final int headerSizeInBytes = BCF2Type.INT32.read(in);
	
	        if ( headerSizeInBytes <= 0 ) throw new IOException("BCF2 header has invalid length: " + headerSizeInBytes + " must be >= 0");
	        final byte[] headerBytes = new byte[headerSizeInBytes];
	        readFully(in, headerBytes);
	        
	        try( final PositionalBufferedStream bps = new PositionalBufferedStream(new ByteArrayInputStream(headerBytes))) {
		        final LineIterator lineIterator = new LineIteratorImpl(new SynchronousLineReader(bps));
		
		        return (VCFHeader)headerParser.readActualHeader(lineIterator);
		       	}
			}
		}
		
		/** Reads len bytes in a loop */
	    private static void readFully(final InputStream in, final byte buf[], int off, int len ) throws IOException {
	        int toRead = len;
	        while ( toRead > 0 ) {
	          int ret = in.read( buf, off, toRead );
	          if ( ret < 0 ) {
	            throw new IOException( "Premeture EOF from inputStream");
	          }
	          toRead -= ret;
	          off += ret;
	        }
	      }
      /** Reads len bytes in a loop.
       * @param in The InputStream to read from
       * @param buf The buffer to fill
       * @throws IOException if it could not read requested number of bytes 
       * for any reason (including EOF)
       */
	    private static void readFully(final InputStream in, final byte buf[]) throws IOException {
	      	readFully(in, buf,0,buf.length);
	      }
	   
}
