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
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
				.take(2)
        '''
        and:
        	def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
        then:
		    result.val[0] == "RF01"
		    result.val[0] == "RF02"
	        result.val == Channel.STOP
        
    }

	def 'dictionary with interval_list' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.interval_list')
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "RF01"
			result.val[0] == "RF02"
			result.val == Channel.STOP
		}
	
		def 'dictionary with compressed interval_list' () {
			when:
			def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.interval_list.gz')
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
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
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getLengthOnReference(),row[1]]}}
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == 3302
			result.val[0] == 2687
			result.val == Channel.STOP
		
	}
	
	def 'dictionary with dict and tid' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.dict')
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getSequenceIndex(),row[1]]}}
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == 0
			result.val[0] == 1
			result.val == Channel.STOP
		
	}
	
	
	def 'dictionary with BAM' () {
		when:
		def SCRIPT = '''
            include {dictionary} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/S1.rota.bam')
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
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
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
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
                .map{[dictionary(it),it]}
                .flatMap{row->row[0].getSequences().collect{dict->[dict.getContig(),row[1]]}}
				.take(2)

        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "RF01"
			result.val[0] == "RF02"
			result.val == Channel.STOP
		
	}
	

	def 'samples with vcf' () {
		when:
		def SCRIPT = '''
            include {samples} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.vcf.gz')
                .map{[samples(it),it]}
                .flatMap{row->row[0].collect{sn->[sn,row[1]]}}
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
                .map{[samples(it),it]}
                .flatMap{row->row[0].collect{sn->[sn,row[1]]}}
				.take(2)
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "S1"
			result.val[0] == "S2"
			result.val == Channel.STOP
		
	}

	def 'samples with sam' () {
		when:
		def SCRIPT = '''
            include {samples} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/S1.rota.bam')
                .map{[samples(it),it]}
                .flatMap{row->row[0].collect{sn->[sn,row[1]]}}
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "S1"
			result.val == Channel.STOP
		
	}
	
	def 'readGroups with sam' () {
		when:
		def SCRIPT = '''
            include {readGroups} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/S1.rota.bam')
                .map{[readGroups(it),it]}
                .flatMap{row->row[0].collect{sn->[sn.getLibrary(),row[1]]}}
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "L1"
			result.val == Channel.STOP
		
	}

	def 'build with sam' () {
		when:
		def SCRIPT = '''
            include {build} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/S1.rota.bam')
                .map{[build(it),it]}
                .filter{it[0]!=null}
				.map{[it[0].getId(),it[1]]}
        '''
		and:
			def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
		then:
			result.val[0] == "rotavirus"
			result.val == Channel.STOP
		
	}
}
