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

    /*
     * {@code reverse} is a `producer` method and will be available to the script because:
     *
     * - it's public
     * - it returns a DataflowWriteChannel
     * - it's marked with the @Factory annotation
     *
     * The method can require arguments but it's not mandatory, it depends of the business logic of the method.
     *
     */
    @Factory
    DataflowWriteChannel reverse(String message) {
        final channel = CH.create()
        session.addIgniter((action) -> reverseImpl(channel, message))
        return channel
    }

    private void reverseImpl(DataflowWriteChannel channel, String message) {
        channel.bind(message.reverse());
        channel.bind(Channel.STOP)
    }

    /*
    * {@code goodbye} is a *consumer* method as it receives values from a channel to perform some logic.
    *
    * Consumer methods are introspected by nextflow-core and include into the DSL if the method:
    *
    * - it's public
    * - it returns a DataflowWriteChannel
    * - it has only one arguments of DataflowReadChannel class
    * - it's marked with the @Operator annotation 
    *
    * a consumer method needs to proportionate 2 closures:
    * - a closure to consume items (one by one)
    * - a finalizer closure
    *
    * in this case `goodbye` will consume a message and will store it as an upper case
    */
    @Operator
    DataflowWriteChannel goodbye(DataflowReadChannel source) {
        final target = CH.createBy(source)
        final next = { target.bind("Goodbye $it".toString());}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }


    private Object findPathWithExtension(target, suffixes, Integer index) {
 	 tuple = (target instanceof List?target:[target])
                if(element_idx != null && element_idx>=0) {
                        return tuple[element_idx]
                        }
                else
                        {
                        int i = 0;
                        while(i< args.size()) {
                                def arg = args[i]
                                if((arg instanceof Path) || (arg instanceof File)) {
					}
                                }
			}
	  }


    @Operator
    DataflowWriteChannel fastaDict(DataflowReadChannel source, Map params = null) {
	Integer element_idx = null;
        if(params==null) params=[:]	
	//validate params
	for(Object k: params.keySet()) {
		//throw new IllegalArgumentException("\""+k+"\" is not a valid key.");
		}


	final target = CH.createBy(source)
        final next = {
		tuple = (it instanceof List?it:[it])
		if(element_idx != null && element_idx>=0) {
			fasta = tuple[element_idx]
			}
		else
			{
			int i = 0;
			while(i< args.size()) {
				def arg = args[i]
				if((arg instanceof Path) || (arg instanceof 
				}		

		final def htsfile = (java.nio.file.Path)it;
		final def dict  = SAMSequenceDictionaryExtractor.extractDictionary(htsfile)
		if(dict==null) {
			throw new SAMException("Cannot extract dictionary from \""+ htsfile + "\". Fasta files must be indexed TODO");
			}
		dict.getSequences().each{V->target.bind(["contig":V.getSequenceName(),"length":V.getSequenceLength(),"file":htsfile])}
		}
        final done = { target.bind(Channel.STOP) }
        DataflowHelper.subscribeImpl(source, [onNext: next, onComplete: done])
        return target
    }


    /*
     * Generate a random string
     *
     * Using @Function annotation we allow this function can be imported from the pipeline script
     */
    @Function
    String randomString(int length=9){
        new Random().with {(1..length).collect {(('a'..'z')).join(null)[ nextInt((('a'..'z')).join(null).length())]}.join(null)}
    }

}
