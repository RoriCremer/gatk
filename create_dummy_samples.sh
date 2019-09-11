#!/usr/bin/env bas


COUNTER=0
       while [  $COUNTER -lt 24 ]; do
           echo $COUNTER
           bq query  --allow_large_results --destination_table=joint_genotyping_chr20_dummy_data_scale_testing.vet  --append_table --project_id=broad-dsp-spec-ops "SELECT position, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, CONCAT(sample, '_d') as sample, call_GT, call_AD, call_DP, call_GQ, call_PGT, call_PID, call_PL FROM `broad-dsp-spec-ops.joint_genotyping_chr20_dalio_40000_july_updated.vet"
           let COUNTER=COUNTER+1
       done


COUNTER=0
       while [  $COUNTER -lt 24 ]; do
           echo $COUNTER
           bq query  --allow_large_results --destination_table=joint_genotyping_chr20_dummy_data_scale_testing.pet_without_gq60 --append_table --project_id=broad-dsp-spec-ops "SELECT position, CONCAT(sample, '_d$COUNTER') as sample, state FROM joint_genotyping_chr20_dalio_40000_july_updated.pet_without_gq60"
           let COUNTER=COUNTER+1
       done


COUNTER=0
       while [  $COUNTER -lt 24 ]; do
           echo $COUNTER
           bq query  --allow_large_results --destination_table=joint_genotyping_chr20_dummy_data_scale_testing.sample_list  --append_table  --project_id=broad-dsp-spec-ops "SELECT CONCAT(sample, '_d$COUNTER') as sample FROM joint_genotyping_chr20_dalio_40000_july_updated.sample_list"
           let COUNTER=COUNTER+1
       done




