
config ?= compileClasspath


all:  check jar assemble

jar: ./plugins/nf-htsjdk/build/libs/nf-htsjdk-0.1.0.jar

./plugins/nf-htsjdk/build/libs/nf-htsjdk-0.1.0.jar : ./plugins/nf-htsjdk/src/main/nextflow/htsjdk/HtsjdkExtension.groovy
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
	./gradlew check


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
	./gradlew ${mm}test
else
	./gradlew ${mm}test --tests ${class}
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
