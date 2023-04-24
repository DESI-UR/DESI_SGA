#!/bin/bash

#- Script to run `prospect_pages`

#- To be run with DESI data (eg at NERSC computing center)
# Requires to have DESI_SPECTRO_REDUX and DESI_ROOT environment variables set.
# If not done so, see https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/NERSC and manually load software stack

#- Usage:
#  ./current_prospect_code    => run test for bright and dark

OUTPUT_ROOT=/global/homes/m/mjkeller/DESI_SGA/TF/fuji_test_outputs
TARGETLIST=/global/homes/m/mjkeller/DESI_SGA/TF/tfuji_IDs.txt

# -------------------------
# Mode: Scan directory tree
# -------------------------

if [[ $1 == sv1 ]] || [[ $1 == '' ]]; then

    if [[ $2 == bright ]] || [[ $2 == '' ]]; then
        echo "------ Bright ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv1_VIbright
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv1 bright \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi


    if [[ $2 == dark ]] || [[ $2 == '' ]]; then
        echo "------ Dark ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv1_VIdark
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv1 dark \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi

fi

if [[ $1 == sv2 ]] || [[ $1 == '' ]]; then

    if [[ $2 == bright ]] || [[ $2 == '' ]]; then
        echo "------ Bright ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv2_VIbright
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv2 bright \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi


    if [[ $2 == dark ]] || [[ $2 == '' ]]; then
        echo "------ Dark ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv2_VIdark
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv2 dark \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi

fi

if [[ $1 == sv3 ]] || [[ $1 == '' ]]; then

    if [[ $2 == bright ]] || [[ $2 == '' ]]; then
        echo "------ Bright ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv3_VIbright
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv3 bright \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi


    if [[ $2 == dark ]] || [[ $2 == '' ]]; then
        echo "------ Dark ------"
        DATADIR=${DESI_SPECTRO_REDUX}/fuji/healpix
        OUTPUTDIR=${OUTPUT_ROOT}/sv3_VIdark
        mkdir ${OUTPUTDIR}
        prospect_pages --datadir ${DATADIR} \
                        --dirtree_type healpix \
                        --survey_program sv3 dark \
                        --target_list ${TARGETLIST} \
                        --with_zcatalog \
                        -o ${OUTPUTDIR} 

    fi

fi