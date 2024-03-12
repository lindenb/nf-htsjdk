package nextflow.htsjdk


import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Channel
import nextflow.Session
import nextflow.extension.CH
import nextflow.extension.DataflowHelper
import nextflow.plugin.extension.Factory
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.Operator
import nextflow.plugin.extension.PluginExtensionPoint
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.SAMException 
import java.util.Arrays;
import htsjdk.samtools.util.FileExtensions;
import nextflow.htsjdk.HtsjdkUtils;





/**
 *
 * @author : Pierre Lindenbaum univ-nantes.fr
 *
 */
@Slf4j
@CompileStatic
class HtsjdkExtension extends PluginExtensionPoint {
    /*
     * A session hold information about current execution of the script
     */
    private Session session

    /*
     * A Custom config extracted from nextflow.config under htsjdk tag
     * nextflow.config
     * ---------------
     * docker{
     *   enabled = true
     * }
     * ...
     * htsjdk{
     *    prefix = 'Mrs'
     * }
     */
     private HtsjdkConfig config

    /*
     * nf-core initializes the plugin once loaded and session is ready
     * @param session
     */
    @Override
    protected void init(Session session) {
        this.session = session
        this.config = new HtsjdkConfig(session.config.navigate('htsjdk') as Map)
    }

	
	private Object bind2(object, List keyValues) {
		if(object instanceof List) {
			def list = []
			for(int i=0;i< keyValues.size();i+=2) {
				list.add(keyValues.get(i+1));
				}
			return List.class.cast(object).plus(list);
			}
		else if(object instanceof Map) {
			def hash = [:]
			for(int i=0;i< keyValues.size();i+=2) {
				hash.put(keyValues.get(i), keyValues.get(i+1));
				}
			return Map.class.cast(object).plus(hash)
			}
		else {
			return bind2([object],keyValues);
			}
		}



    @Operator
    DataflowWriteChannel build(DataflowReadChannel source, Map params = null) {
        if(params==null) params=[:]	
		//validate params
		for(Object k: params.keySet()) {
			if(HtsjdkUtils.KEY_ELEMENT.equals(k)) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}


	final target = CH.createBy(source)
        final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get(HtsjdkUtils.KEY_ELEMENT),null);
			def dict  = htsfile.extractDictionary();
			def build = HtsjdkUtils.findBuild(dict)
			target.bind(bind2(it,[
					"chrom_count",dict.size(),
					"dict_md5",dict.md5(),
					"organism",(build==null?null:build.getOrganism()),
					"build",(build==null?null:build.getId())
					]))	
			}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    	}

	

    @Operator
    DataflowWriteChannel dictionary(DataflowReadChannel source, Map params = null) {
        if(params==null) params=[:]	
		//validate params
		for(Object k: params.keySet()) {
			if(HtsjdkUtils.KEY_ELEMENT.equals(k)) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}


		final target = CH.createBy(source)
        final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get(HtsjdkUtils.KEY_ELEMENT),null);
			def dict  = htsfile.extractDictionary();
			dict.getSequences().each{
				V->target.bind(bind2(it,[
					"tid",V.getSequenceIndex(),
					"chrom",V.getSequenceName(),
					"chromLength", V.getSequenceLength(),
					"altNames", new java.util.ArrayList(V.getAlternativeSequenceNames())
					]))
				}
			}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    	}
		
	@Operator
	DataflowWriteChannel samples(DataflowReadChannel source, Map params = null) {
		def att = null;
		if(params==null) params=[:]
		//validate params
		for(Object k: params.keySet()) {
			if(HtsjdkUtils.KEY_ELEMENT.equals(k)) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}


		final target = CH.createBy(source)
		final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get(HtsjdkUtils.KEY_ELEMENT),null);
			def samples  = htsfile.extractSamples();
			samples.each{
				V->target.bind(bind2(it,[
					"sample",V
					]))
				}
			}
		final done = { target.bind(Channel.STOP) }
		DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
		return target
		}
		
}
