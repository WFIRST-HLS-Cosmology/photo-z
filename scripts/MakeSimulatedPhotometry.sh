#! /bin/bash

# This script makes a new simulated photometry input catalog.
# It calculates photometry for set of objects with a specific set of parameters.

source ./simulation.params
source utils.sh

errlog "\nGenerating simulated photometry model."

# Integrate filters, perform redshift adjustment, and modify SEDs to account for dust

if [ "$STRATOS_HOME" != "" ]; then

    # StratOS has been detected. Split the file into shards, submit the command to the 
    # cluster, then combine the output files. TODO: do the splitting automatically 
    # with the StratOS partitioner and recombine, automatically.

    errlog "StratOS has been detected."

    if [ ! -e "$DFS_DIR" ]; then

         errlog "Creating $DFS_DIR"

         mkdir -p $DFS_DIR
    fi

    errlog "Copying filters and templates to the DFS: $DFS_DIR"

    cp -r $FilterDir $DFS_DIR/

    cp -r $TemplateDir $DFS_DIR/

    shard_base="${InputModels}.part-"

    split -d -n l/$StratOS_Shards $InputModels $shard_base

    errlog "Copying $InputModels to the DFS: $DFS_DIR"

    cp ${shard_base}* $DFS_DIR/

    # make sure that the files are really written to disk

    ansible all -a "sync"

    stratos-run $STRATOS_HOME/bin/make-grid $DFS_DIR/$InputSeds \
                                            $DFS_DIR/$SimulatedDataFilterList \
                                            "%c% %c%.output" \
                                            "$DFS_DIR/$shard_base*"

    errlog "Writing $DFS_DIR/$OutputModels"

    cat $DFS_DIR/${shard_base}*.output > $DFS_DIR/$OutputModels

   # rm $DFS_DIR/${shard_base}* ${shard_base}*
 
    errlog "Copying the photometry file to $DFS_DIR/$OutputModels to simulation directory."

    cp $DFS_DIR/$OutputModels .

else
    errlog "make-grid $TemplateLibrary $FilterList $InputModels $OutputModels"

    make-grid $TemplateLibrary $FilterList $InputModels $OutputModels
fi

cat $OutputModels | add_noise_to_photometry.py 26.1 27.4 27.5 26.8 26.1 24.9 25.4 25.6 25.4 24.7 25.3 > $NoisyPhotometry

