
config ?= compileClasspath


all:  check jar assemble

launch: 
	./launch.sh run  -plugins  nf-htsjdk@0.1.0 data/test01.nf 

jar: ./plugins/nf-htsjdk/build/libs/nf-htsjdk-0.1.0.jar

./plugins/nf-htsjdk/build/libs/nf-htsjdk-0.1.0.jar : ./plugins/nf-htsjdk/src/main/nextflow/htsjdk/HtsjdkExtension.groovy ./plugins/nf-htsjdk/src/main/nextflow/htsjdk/HtsjdkUtils.java
	./gradlew jar

ifdef module 
mm = :${module}:
else 
mm = 
endif 

clean:
	rm -rvf .nextflow*
	rm -rvf work
	rm -rvf build
	rm -rvf plugins/*/build
	./gradlew clean

compile:
	./gradlew :nextflow:exportClasspath compileGroovy
	@echo "DONE `date`"


check:
	./gradlew check --warning-mode all


#
# Show dependencies try `make deps config=runtime`, `make deps config=google`
#
deps:
	./gradlew -q ${mm}dependencies --configuration ${config}

deps-all:
	./gradlew -q dependencyInsight --configuration ${config} --dependency ${module}

#
# Refresh SNAPSHOTs dependencies
#
refresh:
	./gradlew --refresh-dependencies 

#
# Run all tests or selected ones
#
test:
ifndef class
	./gradlew ${mm}test --warning-mode all
else
	./gradlew ${mm}test --tests ${class} --warning-mode all
endif

assemble:
	./gradlew assemble

#
# generate build zips under build/plugins
# you can install the plugin copying manually these files to $HOME/.nextflow/plugins
#
buildPlugins:
	./gradlew copyPluginZip

#
# Upload JAR artifacts to Maven Central
#
upload:
	./gradlew upload


upload-plugins:
	./gradlew plugins:upload

publish-index:
	./gradlew plugins:publishIndex


run:
	$(MAKE)
	$(MAKE) buildPlugins
	./launch.sh data/main.nf
