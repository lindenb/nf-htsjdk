package nextflow.htsjdk

import java.nio.file.Files
import java.util.jar.Manifest

import nextflow.Channel
import nextflow.plugin.Plugins
import nextflow.plugin.TestPluginDescriptorFinder
import nextflow.plugin.TestPluginManager
import nextflow.plugin.extension.PluginExtensionProvider
import org.pf4j.PluginDescriptorFinder
import spock.lang.Shared
import spock.lang.Timeout
import test.Dsl2Spec

import java.nio.file.Path


/**
 * Unit test for Htsjdk DSL
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 */
@Timeout(10)
class HtsjdkDslTest extends Dsl2Spec{

    @Shared String pluginsMode

    def setup() {
        // reset previous instances
        PluginExtensionProvider.reset()
        // this need to be set *before* the plugin manager class is created
        pluginsMode = System.getProperty('pf4j.mode')
        System.setProperty('pf4j.mode', 'dev')
        // the plugin root should
        def root = Path.of('.').toAbsolutePath().normalize()
        def manager = new TestPluginManager(root){
            @Override
            protected PluginDescriptorFinder createPluginDescriptorFinder() {
                return new TestPluginDescriptorFinder(){

                    @Override
                    protected Manifest readManifestFromDirectory(Path pluginPath) {
                        if( !Files.isDirectory(pluginPath) )
                            return null

                        final manifestPath = pluginPath.resolve('build/resources/main/META-INF/MANIFEST.MF')
                        if( !Files.exists(manifestPath) )
                            return null

                        final input = Files.newInputStream(manifestPath)
                        return new Manifest(input)
                    }
                }
            }
        }
        Plugins.init(root, 'dev', manager)
    }

    def cleanup() {
        Plugins.stop()
        PluginExtensionProvider.reset()
        pluginsMode ? System.setProperty('pf4j.mode',pluginsMode) : System.clearProperty('pf4j.mode')
    }


    def 'dictionary with dict' () {
        when:
        def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.dict')
                .dictionary()
				.take(2)
        '''
        and:
        	def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
        then:
		    result.val[0] == "RF01"
		    result.val[0] == "RF02"
	        result.val == Channel.STOP
        
    }

	def 'dictionary with dict and length' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.dict')
                .dictionary(withLength:true)
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[1] == 3302
			result.val[1] == 2687
			result.val == Channel.STOP
		
	}
	
	def 'dictionary with dict and tid' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.dict')
                .dictionary(withTid:true)
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[1] == 0
			result.val[1] == 1
			result.val == Channel.STOP
		
	}
	
	
	def 'dictionary with BAM' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/S1.rota.bam')
                .dictionary()
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "RF01"
			result.val[0] == "RF02"
			result.val == Channel.STOP
		
	}
	
	def 'dictionary with VCF' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.vcf.gz')
                .dictionary()
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "RF01"
			result.val[0] == "RF02"
			result.val == Channel.STOP
		
	}
	
	def 'dictionary with BCF' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.bcf')
                .dictionary()
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "RF01"
			result.val[0] == "RF02"
			result.val == Channel.STOP
		
	}
	
	def 'dictionary with dict with header=true' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.dict')
                .dictionary(header:true)
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0].chrom == "RF01"
			result.val[0].chrom == "RF02"
			result.val == Channel.STOP
		
	}

	def 'samples with vcf' () {
		when:
		def SCRIPT = '''
            include {samples} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.vcf.gz')
                .samples()
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "S1"
			result.val[0] == "S2"
			result.val == Channel.STOP
		
	}
	
	def 'samples with bcf' () {
		when:
		def SCRIPT = '''
            include {samples} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.bcf')
                .samples()
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "S1"
			result.val[0] == "S2"
			result.val == Channel.STOP
		
	}
}
