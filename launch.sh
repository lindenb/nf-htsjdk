export NXF_PLUGINS_DEV=$PWD/plugins

 ../nextflow/launch.sh "$@" || \
${HOME}/src/nextflow/launch.sh "$@" 
