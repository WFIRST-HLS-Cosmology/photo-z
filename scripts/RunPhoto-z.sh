#! /bin/bash

source ./simulation.params
source utils.sh

errlog "\nPreparing to run the Photo-z fitting code."


if [ "$STRATOS_HOME" != "" ]; then

    # StratOS has been detected. Split the file into shards, submit the command to the 
    # cluster, then combine the output files. TODO: do the splitting automatically 
    # with the StratOS partitioner and recombine, automatically.

    ModelGridGlobPattern="${ModelGridSpec}.part-*.photoz"

    shard_base="${NoisyPhotometry}.part-"

    split -d -n l/$StratOS_Shards $NoisyPhotometry $shard_base

    errlog "Copying $NoisyPhotometry to the DFS: $DFS_DIR"

    cp ${shard_base}* $DFS_DIR/

    errlog "Copying $TemplateProbabilities to the DFS: $DFS_DIR"

    cp $TemplateProbabilities $DFS_DIR/

    # make sure that the files are really written to disk (not cached)

    ansible all -a "sync"

    sleep 10 # poor coding here: this is intended to give FUSE-DFS + HDFS time to update. 

    #TODO: possibly check to make sure that all nodes are able to read all of the shards.

    sleep 3

    stratos-run $STRATOS_HOME/bin/photo-z "$DFS_DIR/$ModelGridGlobPattern" \
                                          "%c% $DFS_DIR/$TemplateProbabilities %c%.output" \
                                          "$DFS_DIR/$shard_base*"

    errlog "Saving estimated redshifts in $DFS_DIR/$OutputRedshiftCatalog"

    cat $DFS_DIR/${shard_base}*.output > $DFS_DIR/$OutputRedshiftCatalog

    rm $DFS_DIR/${shard_base}* ${shard_base}*

    errlog "Copying $OutputRedshiftCatalog to simulation directory."

    cp $DFS_DIR/$OutputRedshiftCatalog .

else
    # run locally

    nlines=($(wc -l $NoisyPhotometry))

    if ((nlines > 100000)); then

        nslices=$(($nlines / 50000))

        split -n l/$nslices $NoisyPhotometry "$NoisyPhotometry-part-"

        parts=$(ls $NoisyPhotometry-part-*)

        first_file=1

        for file in $parts; do

            echo "photo-z ${ModelGrid}.photoz $file $TemplateProbabilities $file-output"

            photo-z "${ModelGrid}.photoz" $file $TemplateProbabilities $file-output

            # only include the column headers once
            if (($first_file == 1)); then

                head -n 1 $file-output > $OutputRedshiftCatalog

                first_file=0

            fi

            tail -n +2 $file-output >> $OutputRedshiftCatalog

        done

        #cat $NoisyPhotometry-part-*-output > $OutputRedshiftCatalog

        rm $NoisyPhotometry-part-*

    else

        echo "photo-z ${ModelGrid}.photoz $NoisyPhotometry $TemplateProbabilities $OutputRedshiftCatalog"

        photo-z "${ModelGrid}.photoz" $NoisyPhotometry $TemplateProbabilities $OutputRedshiftCatalog

    fi

fi

make_catalog.py --in_catalog $InputModels --z_catalog $OutputRedshiftCatalog --out_catalog tmpcatalog 

cat tmpcatalog | grep -v ' -9 ' | grep -v 'inf'  > $FinalCatalog
