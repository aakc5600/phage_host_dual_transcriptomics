#!/usr/bin/env nextflow

include { AUTO_PROCESSING } from '../../modules/downstreamProcessing'
include { ANNOTATE_PHA } from '../../modules/downstreamProcessing'
include { ANNOTATE_PHO } from '../../modules/downstreamProcessing'
include { ANNOTATE_PHY } from '../../modules/downstreamProcessing'

workflow PROCESSING {
    take:
    auto
    gbkPath
    phageGenome
    pha
    pho
    phy
    pha_db
    pho_db
    phy_model
    threads
    countTable
    metaPath
    dualGFF
    sampleDict
    timeUntilMiddle
    timeUntilLate
    pwd

    main:
    if ( auto ) {
        /* pha_db = "-d $pha_db"
        pho_db = "-d $pho_db"
        phy_model = "-m $phy_model" */
        if ( gbkPath != "false" ) {
            processing = AUTO_PROCESSING(countTable, metaPath, dualGFF, gbkPath, sampleDict, timeUntilMiddle, timeUntilLate, pwd)
        }
        else if ( phy ) {
            pha = ANNOTATE_PHA(phageGenome, pha_db, threads)
            pho = ANNOTATE_PHO(pha.out, pho_db, threads)
            phy = ANNOTATE_PHY(pho.out, phy_model)
            processing = AUTO_PROCESSING(countTable, metaPath, dualGFF, phy.out, sampleDict, timeUntilMiddle, timeUntilLate, pwd)
        }
        else if ( pho ) {
            pha = ANNOTATE_PHA(phageGenome, pha_db, threads)
            pho = ANNOTATE_PHO(pha.out, pho_db, threads)
            processing = AUTO_PROCESSING(countTable, metaPath, dualGFF, pho.out, sampleDict, timeUntilMiddle, timeUntilLate, pwd)
        }
        else if ( pha ) {
            pha = ANNOTATE_PHA(phageGenome, pha_db, threads)
            processing = AUTO_PROCESSING(countTable, metaPath, dualGFF, pha.out, sampleDict, timeUntilMiddle, timeUntilLate, pwd)
        }
        else {
            processing = AUTO_PROCESSING(countTable, metaPath, dualGFF, "false", sampleDict, timeUntilMiddle, timeUntilLate, pwd)
        }
    }
    else {
        processing = channel.empty()
    }

    emit:
    analysis = processing
    pha = pha
    pho = pho
    phy = phy
}
