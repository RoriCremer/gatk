#!/usr/bin/env bas


PROJECT_ID="spec-ops-sand-box"
#"broad-dsp-spec-ops"
ORIGIN_DATASET="chr20"
DESTINATION_DATASET="chr20_dummy_million"


VET_COUNTER=0
       while [  $VET_COUNTER -lt 24 ]; do
           echo $VET_COUNTER
           bq query  --allow_large_results --destination_table=$DESTINATION_DATASET.vet  --append_table --project_id=$PROJECT_ID \
           "SELECT position, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, CONCAT(sample, '_d') as sample, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL FROM $ORIGIN_DATASET.vet"
           let VET_COUNTER=VET_COUNTER+1
       done


PET_COUNTER=0
       while [  $PET_COUNTER -lt 24 ]; do
           echo $PET_COUNTER
           bq query  --allow_large_results --destination_table=$DESTINATION_DATASET.pet --append_table --project_id=$PROJECT_ID \
           "SELECT position, CONCAT(sample, '_d$COUNTER') as sample, state FROM $ORIGIN_DATASET.pet"
           let PET_COUNTER=PET_COUNTER+1
       done


COUNTER=0
       while [  $COUNTER -lt 24 ]; do
           echo $COUNTER
           bq query  --allow_large_results --destination_table=$DESTINATION_DATASET.sample_list  --append_table  --project_id=$PROJECT_ID \
           "SELECT CONCAT(sample, '_d$COUNTER') as sample FROM metadata.sample_list"
           let COUNTER=COUNTER+1
       done
