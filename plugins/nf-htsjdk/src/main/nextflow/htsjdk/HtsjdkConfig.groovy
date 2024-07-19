package nextflow.htsjdk

import groovy.transform.PackageScope


/**
 * This class allows model an specific configuration, extracting values from a map and converting
 *
 * In this plugin, the user can configure how the messages are prefixed with a String, i.e.
 * due a nextflow.config
 *
 * htsjdk {
 *     prefix = '>>'
 * }
 *
 * when the plugin reverse a String it will append '>>' at the beginning instead the default 'Mr.'
 *
 * We anotate this class as @PackageScope to restrict the access of their methods only to class in the
 * same package
 *
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@PackageScope
class HtsjdkConfig {

    final private List<HtsjdkUtils.Build> builds

    HtsjdkConfig(Map map){
        def config = map ?: Collections.emptyMap()
        this.builds = config.containsKey("builds") 
		? HtsjdkUtils.decodeBuilds(config.get("builds"))
		: HtsjdkUtils.getDefaultBuilds();
    }

    List<HtsjdkUtils.Build> getBuilds() { 
		return this.builds;
		}
}
