#! /bin/bash

# ==================== MakeModelGrid ====================

# computes the model grid. This is expensive, so we only do it if the filters,
# templates, or wavelength range has changed.

source ./simulation.params
source utils.sh

errlog "Generating the model grid."

# if the output file already exists and it is not older than the $ModelGridSpec and 
# $TemplateLibrary files, then use the existing grid. 

if test -e "${ModelGrid}.photoz"; then

    errlog "Note: The file $ModelGrid already exists. I'll use it instead of re-generating it."
    exit
fi

if [ "$STRATOS_HOME" != "" ]; then

    # StratOS has been detected. Split the file into shards, submit the command to the 
    # cluster, then combine the output files. TODO: do the splitting automatically 
    # with the StratOS partitioner and recombine, automatically.

    errlog "Copying filters and templates to the DFS: $DFS_DIR"

    cp -r $FilterDir $DFS_DIR/

    cp -r $TemplateDir $DFS_DIR/

    shard_base="${ModelGridSpec}.part-"

    split -d -n l/$StratOS_Shards $ModelGridSpec $shard_base

    errlog "Copying $ModelGridSpec to the DFS: $DFS_DIR"

    cp ${shard_base}* $DFS_DIR/

    # make sure that the files are really written to disk (not cached)

    ansible all -a "sync"

    stratos-run $STRATOS_HOME/bin/make-grid $DFS_DIR/$TemplateLibrary \
                                            $DFS_DIR/$FittingDataFilterList \
                                            "%c% %c%.photoz" \
                                            "$DFS_DIR/$shard_base*"

    errlog "Writing the model grid in $(ls $DFS_DIR/*.photoz)"

    rm ${shard_base}*

else
    errlog "Running: make-grid $TemplateLibrary $FilterList $ModelGridSpec ${ModelGrid}.photoz"   
    make-grid $TemplateLibrary $FilterList $ModelGridSpec ${ModelGrid}.photoz
fi

if [[ $? == 0 ]]; then
    # TODO: verify that the output file exists and is non-empty.
    errlog "Finished generating the model grid."
else
    errlog "Grid generation failed."
    exit 1
fi
