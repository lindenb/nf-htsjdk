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


    def 'test with chromosomes' () {
        when:
        def SCRIPT = '''
            include {chromosomes} from 'plugin/nf-htsjdk'
            channel
                .fromPath('../../data/rotavirus_rf.fa')
                .chromosomes()
				.take(1)
        '''
        and:
        	def result = new MockScriptRunner([:]).setScript(SCRIPT).execute()
        then:
		    result.val[0] == "RF01"
	        result.val == Channel.STOP
        
    }

}
