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
					"dict_length",dict.getReferenceLength(),
					"organism",(build==null?null:build.getOrganism()),
					"build",(build==null?null:build.getId())
					]))	
			}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    	}

		@Operator
		DataflowWriteChannel chromosomes(DataflowReadChannel source, Map params = null) {
			if(params==null) params=[:]
			//validate params
			for(Object k: params.keySet()) {
				throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
				}
			
			final target = CH.createBy(source)
			final next = {
				def htsfile = HtsjdkUtils.findHtsSource(it,null,null);
				final SAMSequenceDictionary dict  = htsfile.extractDictionary();
				
				dict.getSequences().each{
					V->target.bind(
						[
						V.getSequenceName(),
						it
						]
						);
					}
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
			if(k.equals("header")) continue;
			if(k.equals("elem")) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}
		final Object withHeaderObject = params.getOrDefault("header", true);
		if(!(withHeaderObject instanceof Boolean)) throw new IllegalArgumentException("\"header\" is not a valid key.");
		final boolean withHeader = Boolean.class.cast(withHeaderObject);
		final Object elem = params.getOrDefault("elem", null);
		

		final target = CH.createBy(source)
        final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get("elem"),null);
			final SAMSequenceDictionary dict  = htsfile.extractDictionary();
			if(withHeader) {
				/*
				dict.getSequences().each{
					V->target.bind(bind2([
						"tid",V.getSequenceIndex(),
						"chrom",V.getSequenceName(),
						"chromLength", V.getSequenceLength(),
						"altNames", new java.util.ArrayList(V.getAlternativeSequenceNames())
						],it));
						}
				*/
				}
			else
				{
				dict.getSequences().each{
					V->target.bind(
						[
						V.getSequenceIndex(),
						V.getSequenceName(),
						V.getSequenceLength(),
						new java.util.ArrayList(V.getAlternativeSequenceNames()),
						it
						]
						);
					}
				}
			}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    	}
		
	@Operator
	DataflowWriteChannel samples(DataflowReadChannel source, Map params = null) {
		if(params==null) params=[:]
		def enableEmpty = (params.getOrDefault(HtsjdkUtils.KEY_ENABLE_EMPTY,false) as boolean);
		//validate params
		for(Object k: params.keySet()) {
			if(HtsjdkUtils.KEY_ELEMENT.equals(k)) continue;
			if(HtsjdkUtils.KEY_ENABLE_EMPTY.equals(k)) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}


		final target = CH.createBy(source)
		final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get(HtsjdkUtils.KEY_ELEMENT),null);
			def samples  = htsfile.extractSamples();
			if(enableEmpty  && samples.isEmpty()) {
				target.bind(bind2(it,[
					"sample",HtsjdkUtils.NO_SAMPLE
					]))
				}
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
	
	@Operator
	DataflowWriteChannel mappedContigs(DataflowReadChannel source, Map params = null) {
		if(params==null) params=[:]
		def enableEmpty = (params.getOrDefault(HtsjdkUtils.KEY_ENABLE_EMPTY,false) as boolean);
		//validate params
		for(Object k: params.keySet()) {
			if(HtsjdkUtils.KEY_ELEMENT.equals(k)) continue;
			if(HtsjdkUtils.KEY_ENABLE_EMPTY.equals(k)) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}


		final target = CH.createBy(source)
		final next = {
			def htsfile = HtsjdkUtils.findHtsSource(it, params.get(HtsjdkUtils.KEY_ELEMENT),null);
			def samples  = htsfile.extractMappedContigs();
			if(enableEmpty  && samples.isEmpty()) {
				target.bind(bind2(it,[
					"contig",HtsjdkUtils.NO_CONTIG
					]))
				}
			samples.each{
				V->target.bind(bind2(it,[
					"contig",V
					]))
				}
			}
		final done = { target.bind(Channel.STOP) }
		DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
		return target
		}	
}
