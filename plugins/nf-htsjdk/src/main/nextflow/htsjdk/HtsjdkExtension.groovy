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
import htsjdk.samtools.util.StringUtil
import nextflow.htsjdk.HtsjdkUtils;
import nextflow.htsjdk.HtsjdkUtils.Build


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

	
	private HtsjdkConfig getConfig() {
		return this.config;
		}
	
	private Object bind2(Object userData,Object row) {
		if(row instanceof List) {
			List L=[];
			
			if(userData instanceof List) {
				L.addAll(List.class.cast(userData));
				}
			else
				{
				L.add(userData);
				}
			L.addAll(List.class.cast(row));
			return L;
			}
		else if(row instanceof Map && userData instanceof Map) {
			Map hash=[:];
			hash.putAll(Map.class.cast(row));
			hash.putAll(Map.class.cast(userData));
			return hash;
			}
		else {
			return bind2(userData,[row]);
			}
		}

    @Function
    Object dictionary(Object source, Map params = null) {
        if(params==null) params=[:]	
		//validate params
		for(Object k: params.keySet()) {
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}		
		final HtsjdkUtils.HtsSource htsfile = HtsjdkUtils.findHtsSource(source ,{HTS->HTS.isBamCramSam() || HTS.isVcf() || HTS.isDict()|| HTS.isFai()| HTS.isFasta() || HTS.isIntervalList()});
		return htsfile.extractDictionary();
    	}
	
    @Function
	Object build(Object source, Map params = null) {
		if(params==null) params=[:]
		final SAMSequenceDictionary dict= SAMSequenceDictionary.class.cast( dictionary(source,[:]) );
		if(dict==null) return null;
		//validate params
		for(Object k: params.keySet()) {
			if(k.equals("resolveContig")) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}
		final boolean resolveContigName = params.containsKey("resolveContig")
			? (params.get("resolveContig") as boolean)
			: getConfig().isResolveContigName()
			;
		
		return this.findBuild(resolveContigName,dict);
		}

	@Function
	Collection readGroups(Object source, Map params = null) {
		if(params==null) params=[:]
		//validate params
		for(Object k: params.keySet()) {
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}
						
		final HtsjdkUtils.HtsSource htsfile = HtsjdkUtils.findHtsSource(source,{HTS->HTS.isBamCramSam() || HTS.isIntervalList() });
		return htsfile.extractReadGroups();
		}
		
	@Function
	Collection samples(Object source, Map params = null) {
		if(params==null) params=[:]	
		//validate params
		for(Object k: params.keySet()) {
			if(k.equals("defaultName")) continue;
			throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
			}
		
		final Object defaultName = params.getOrDefault("defaultName", null);
				
		final HtsjdkUtils.HtsSource htsfile = HtsjdkUtils.findHtsSource(source,{HTS->HTS.isBamCramSam() || HTS.isVcf() || HTS.isIntervalList()});
		def samples = htsfile.extractSamples();
		if(samples.isEmpty() && defaultName!=null) {
			samples = Collections.singletonList(defaultName);
			}
		return samples;
		}
		
	private Build findBuild(boolean resolveContig,final SAMSequenceDictionary dict) {
		return getConfig().getBuilds().stream().filter(B->B.match(resolveContig,dict)).findFirst().orElse(null);
		}
	}
