/*
 * Copyright (C) 2017 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */


#include <htslib/vcf.h>
#include <htslib/sam.h>
#include "CuTest.h"
#include "sonLib.h"
#include "stRPHmm.h"
#include "../externalTools/sonLib/C/impl/sonLibListPrivate.h"

void printSequenceStats(stList *profileSequences) {
    /*
     * Print stats about the set of profile sequences.
     */
    int64_t totalLength=0;
    for(int64_t i=0; i<stList_length(profileSequences); i++) {
        stProfileSeq *profileSeq = stList_get(profileSequences, i);
        totalLength += profileSeq->length;
    }
    st_logDebug("Got %" PRIi64 " profile sequences, with total length: %" PRIi64 ", average length: %f\n",
            stList_length(profileSequences), totalLength, ((float)totalLength)/stList_length(profileSequences));
}

double getExpectedNumberOfMatches(uint64_t *haplotypeString, int64_t start, int64_t length, stProfileSeq *profileSeq) {
    /*
     * Returns the expected number of positions in the profile sequence that are identical to the given haplotype string.
     */
    double totalExpectedMatches = 0.0;

    for(int64_t i=0; i<profileSeq->length; i++) {
        // Get base in the haplotype sequence
        int64_t j = i + profileSeq->refStart - start;
        if(j >= 0 && j < length) {
            uint64_t hapBase = haplotypeString[j];
            assert(hapBase < ALPHABET_SIZE);

            // Expectation of a match
            totalExpectedMatches += getProb(&(profileSeq->profileProbs[i * ALPHABET_SIZE]), hapBase);
        }
    }
    return totalExpectedMatches;
}

double getExpectedIdentity(uint64_t *haplotypeString, int64_t start, int64_t length, stSet *profileSeqs) {
    /*
     * Returns the expected fraction of positions in the profile sequences that match their corresponding position in the
     * given haplotype string.
     */
    double totalExpectedNumberOfMatches = 0.0;
    int64_t totalLength = 0;
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        totalExpectedNumberOfMatches += getExpectedNumberOfMatches(haplotypeString, start, length, pSeq);
        totalLength += pSeq->length;
    }
    stSet_destructIterator(it);
    return totalExpectedNumberOfMatches/totalLength;
}

double getIdentityBetweenHaplotypes(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        }
    }
    return ((double)totalMatches)/length;
}

double getIdentityBetweenHaplotypesExcludingIndels(uint64_t *hap1String, uint64_t *hap2String, int64_t length) {
    /*
     * Returns the fraction of positions in two haplotypes that are identical.
     */
    int64_t totalMatches = 0;
    int64_t numGaps = 0;
    for(int64_t i=0; i<length; i++) {
        if(hap1String[i] == hap2String[i]) {
            totalMatches++;
        } else if (hap1String[i] == ALPHABET_SIZE-1 || hap2String[i] == ALPHABET_SIZE-1) {
            numGaps++;
        }
    }
    return ((double)totalMatches)/(length - numGaps);
}

void getExpectedMatcheBetweenProfileSeqs(stProfileSeq *pSeq1, stProfileSeq *pSeq2, int64_t *totalAlignedPositions, double *totalExpectedMatches) {
    /*
     * Calculates the number of base overlaps and expected base matches between two profile sequences.
     */

    for(int64_t i=0; i<pSeq1->length; i++) {
        // Establish if the coordinate is in both sequences
        int64_t j = i + pSeq1->refStart - pSeq2->refStart;
        if(j >= 0 && j < pSeq2->length) {
            (*totalAlignedPositions)++;

            // Calculate expectation of match
            for(int64_t k=0; k<ALPHABET_SIZE; k++) {
                double e1 = getProb(&(pSeq1->profileProbs[i * ALPHABET_SIZE]), k);
                double e2 = getProb(&(pSeq2->profileProbs[j * ALPHABET_SIZE]), k);
                assert(e1 * e2 <= 1.0);
                *totalExpectedMatches += e1 * e2;
            }
        }
    }
}

void printAvgIdentityBetweenProfileSequences(stList *profileSequences, int64_t maxSeqs) {
    /*
     * Prints the average base identity between pairwise base overlaps between the given profile sequences
     */

    double totalExpectedMatches = 0.0;
    int64_t totalAlignedPositions = 0;

    for(int64_t i=0; i<stList_length(profileSequences) && i<maxSeqs; i++) {
        for(int64_t j=i+1; j<stList_length(profileSequences) && j<maxSeqs; j++) {
            getExpectedMatcheBetweenProfileSeqs(stList_get(profileSequences, i), stList_get(profileSequences, j),
                                                &totalAlignedPositions, &totalExpectedMatches);

        }
    }

    st_logDebug("Avg. pairwise identity between profile sequences: %f measured at %" PRIi64 " overlapping sites\n",
            totalExpectedMatches/totalAlignedPositions, totalAlignedPositions);
}

double *getHaplotypeBaseComposition(uint64_t *hapString, int64_t length) {
    /*
     * Get the count of each alphabet character in the haplotype sequence, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for(int64_t i=0; i<length; i++) {
        baseCounts[hapString[i]] += 1;
    }
    return baseCounts;
}

double *getExpectedProfileSequenceBaseComposition(stSet *profileSeqs) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        for(int64_t i=0; i<pSeq->length; i++) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[i*ALPHABET_SIZE]), j);
            }
        }
    }
    stSet_destructIterator(it);
    return baseCounts;
}

double *getProfileSequenceBaseCompositionAtPosition(stSet *profileSeqs, int64_t pos) {
    /*
     * Get the expected count of each alphabet character in the profile sequences, returned
     * as an array.
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    stSetIterator *it = stSet_getIterator(profileSeqs);
    stProfileSeq *pSeq;
    while((pSeq = stSet_getNext(it)) != NULL) {
        if (pos > pSeq->refStart && pos < pSeq->refStart+pSeq->length) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(pSeq->profileProbs[(pos - pSeq->refStart)*ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;
}

double *getColumnBaseComposition(stRPColumn *column, int64_t pos) {
    /*
     * Get the observed counts for each base seen at a particular position in a column
     */
    double *baseCounts = st_calloc(ALPHABET_SIZE, sizeof(double));
    for (int64_t i=0; i<column->depth; i++) {
        stProfileSeq *seq = column->seqHeaders[i];

        if (pos >= seq->refStart && pos < seq->length+seq->refStart) {
            for(int64_t j=0; j<ALPHABET_SIZE; j++) {
                baseCounts[j] += getProb(&(seq->profileProbs[(pos - seq->refStart) * ALPHABET_SIZE]), j);
            }
        }
    }
    return baseCounts;

}

void printBaseComposition(double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("Base %" PRIi64 " count: %f fraction: %f\n", i, baseCounts[i], baseCounts[i]/totalCount);
    }
}

void printBaseComposition2(double *baseCounts) {
    /*
     * Print the counts/fraction of each alphabet character in a slightly more compressed form.
     */
    double totalCount = 0;
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        totalCount += baseCounts[i];
    }
    st_logDebug("\t\t0 (A)\t1 (C)\t2 (G)\t3 (T)\t4 (-) \n");
    st_logDebug("    Counts:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.1f", baseCounts[i]);
    }
    st_logDebug("\n    Fraction:");
    for(int64_t i=0; i<ALPHABET_SIZE; i++) {
        st_logDebug("\t%0.3f", baseCounts[i]/totalCount);
    }
    st_logDebug( "\n");
}


void printColumnAtPosition(stRPHmm *hmm, int64_t pos) {
    /*
     * Print out the columns of the hmm at a specific position
     */
    stRPColumn *column = hmm->firstColumn;
    while(1) {
        if (pos >= column->refStart && pos < column->refStart+column->length) {
            double *columnBaseCounts = getColumnBaseComposition(column, pos);
            printBaseComposition2(columnBaseCounts);
            free(columnBaseCounts);
        }
        if (column->nColumn == NULL) {
            break;
        }
        column = column->nColumn->nColumn;
    }
}

void printPartitionComposition(int64_t evalPos, stRPHmm *hmm, stSet *reads1, stSet *reads2) {
    double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, evalPos);
    double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, evalPos);
    printColumnAtPosition(hmm, evalPos);
    st_logDebug("\tPartition 1: \n");
    printBaseComposition2(read1BaseCounts);
    st_logDebug("\tPartition 2: \n");
    printBaseComposition2(read2BaseCounts);
    free(read1BaseCounts);
    free(read2BaseCounts);
}

/*
 * Stores information about relevant test results.
 */
struct genotypeResults {
    // Variants in reference
    int64_t negatives;
    int64_t positives;

    // Variants in evaluated vcf
    int64_t truePositives;
    int64_t falsePositives;
    int64_t trueNegatives;
    int64_t falseNegatives;

    // Types of errors
    int64_t error_trueVariantWrong;
    int64_t error_badPartition;
    int64_t error_missedIndels;
};

/*
 * Prints information contained in genotypeResults struct.
 */
void printGenotypeResults(struct genotypeResults *results) {
    // Sensitivity
    st_logInfo("\nSensitivity: %f \n(= fraction of true positives compared to reference, \t%"
                       PRIi64 " out of %"PRIi64 ")\n",
            (float)results->truePositives/results->positives,
               results->truePositives, results->positives) ;
    st_logInfo("\t \t(Number of false negatives: %" PRIi64 ")\n", results->falseNegatives);
    st_logInfo("Sensitivity ignoring false negatives not supported by reads:\t %f \n",
            (float)results->truePositives/(results->positives-results->error_trueVariantWrong)) ;

    // Specificity
    st_logInfo("\nSpecificity: %f \n(= fraction of true negatives compared to reference, \t%"
                       PRIi64 " out of % "PRIi64 ")\n",
            (float)results->trueNegatives/results->negatives,
               results->trueNegatives, results->negatives);
    st_logInfo("\t \t(Number of false positives: %" PRIi64 ")\n", results->falsePositives);

    // More detailed numbers about errors
    st_logInfo("\nFalse negatives:\n");
    st_logInfo("Partition bad: %" PRIi64 " \t\t\t\t(%f)\n",
               results->error_badPartition, (float)results->error_badPartition/results->falseNegatives);

    st_logInfo("Involve indels in the ref vcf: %" PRIi64 " \t\t(%f)\n",
               results->error_missedIndels, (float)results->error_missedIndels/results->falseNegatives);

    st_logInfo("True variant not supported by reads: %" PRIi64 " \t(%f)\n",
               results->error_trueVariantWrong,
               (float)results->error_trueVariantWrong/results->falseNegatives);
}

/*
     * Test to compare a vcf to a truth vcf containing known variants for the region.
     *
     * Test depends on the format of the vcf files written in vcfWriter.c
     * (Currently they don't follow a quite standard format)
     *
     * Information about some of the results saved in the genotypeResults struct
     *
     */
void compareVCFs(FILE *fh, stList *hmms,
                 char *vcf_toEval, char *vcf_ref, double threshold,
                 stBaseMapper *baseMapper, struct genotypeResults *results) {

    st_logInfo("VCF reference: %s \n", vcf_ref);
    st_logInfo("VCF being evaluated: %s \n", vcf_toEval);

    vcfFile *inRef = vcf_open(vcf_ref,"r"); //open vcf file
    if (inRef == NULL) {
        st_logCritical("ERROR: cannot open reference vcf, %s\n", vcf_ref);
        return;
    }
    bcf_hdr_t *hdrRef = bcf_hdr_read(inRef); //read header
    bcf1_t *refRecord = bcf_init1(); //initialize for reading

    vcfFile *inEval = vcf_open(vcf_toEval,"r"); //open vcf file
    if (inEval == NULL) {
        st_logCritical("ERROR: cannot open vcf to evaluate, %s\n", vcf_toEval);
        return;
    }
    bcf_hdr_t *hdrEval = bcf_hdr_read(inEval); //read header
    bcf1_t *evalRecord = bcf_init1(); //initialize for reading
    int64_t referencePos = 0;

    st_logInfo("> Comparing vcf files \n");

    // Start by looking at the first hmm
    int64_t hmmIndex = 0;
    stRPHmm *hmm = stList_get(hmms, hmmIndex);

    // Inefficient, but recalculate the info relevant to the hmm to get bipartitions
    stRPHmm_forwardBackward(hmm);
    stList *path = stRPHmm_forwardTraceBack(hmm);
    stGenomeFragment *gF = stGenomeFragment_construct(hmm, path);
    stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
    stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);

    // Iterate through the vcf being checked until getting to the start of the specified interval
    // Don't bother analyzing these records
    int64_t refStart = 0;
    int64_t evalPos =  0;
    bcf1_t *unpackedRecord;
    int recordCount = 0;

    while(bcf_read(inRef, hdrRef, refRecord) == 0) {
        // To take care of the case where a false positive may have been skipped
        // over if the previous eval location was a false negative
        bool maybeFalsePositive = false;
        if (referencePos < evalPos) {
            maybeFalsePositive = true;
        }
        // Unpack reference record
        bcf1_t *unpackedRecordRef = refRecord;
        bcf_unpack(unpackedRecordRef, BCF_UN_INFO);
        referencePos = unpackedRecordRef->pos+1;

        if (maybeFalsePositive && evalPos < referencePos) {
            results->falsePositives++;
            char *evalRefChar = unpackedRecord->d.als;
            char *altEvalChar = unpackedRecord->d.allele[1];
            st_logDebug("\nFALSE POSITIVE \n\t pos: %" PRIi64 " ref:%s alt: %s \n",
                        evalPos, evalRefChar, altEvalChar);
            printColumnAtPosition(hmm, evalPos);
        }

        // Skip to the first known location of variation in file being evaluated
        if (results->positives == 0) refStart = referencePos;

        // Make sure to only look at records in the specified interval
        if (referencePos < hmm->refStart) continue;

        // If the position is beyond the end of this hmm, get the next one
        while ((hmm->refStart + hmm->refLength) < referencePos) {
            hmmIndex++;
            if (hmmIndex < stList_length(hmms)) {
                hmm = stList_get(hmms, hmmIndex);
                path = stRPHmm_forwardTraceBack(hmm);
                gF = stGenomeFragment_construct(hmm, path);
                reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, 1);
                reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, 0);
            } else {
                break;
            }
        }

        results->positives++;
        char *refChar = unpackedRecordRef->d.als;
        char *refAltChar = unpackedRecordRef->d.allele[1];


        // Iterate through vcf until getting to the position of the variant
        // from the reference vcf currently being looked at
        while (evalPos < referencePos) {

            if (bcf_read(inEval, hdrEval, evalRecord) != 0) {
                st_logInfo("Error: bcf_read\n");
                break;
            }

            unpackedRecord = evalRecord;                // unpack record
            bcf_unpack(unpackedRecord, BCF_UN_INFO);
            evalPos = unpackedRecord->pos+1;
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];
            if (evalPos < refStart) continue;           // skip this record
            recordCount++;

            // Check for false positives - variations found not in reference
            if (evalPos < referencePos) {
                results->falsePositives++;
                st_logDebug("\nFALSE POSITIVE \n\t pos: %" PRIi64 " ref:%s alt: %s \n",
                            evalPos+1, evalRefChar, evalAltChar);
                printColumnAtPosition(hmm, evalPos);
            } else {
                break;
            }
        }
        // At locus of known variation
        if (evalPos == referencePos) {
            char *evalRefChar = unpackedRecord->d.als;
            char *evalAltChar = unpackedRecord->d.allele[1];

            if (strcmp(refChar, evalRefChar) == 0 && strcmp(evalAltChar, refAltChar) == 0) {
                results->truePositives++;
            } else if (strcmp(refChar, evalAltChar) == 0 && strcmp(evalRefChar, refAltChar) == 0) {
                results->truePositives++;
            } else {
                results->falsePositives++;
                st_logDebug("INCORRECT POSITIVE \n");
                st_logDebug("\t pos: %" PRIi64 "\n\tref: %s\talt: ", referencePos, refChar);
                for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
                    if (i != 1) st_logDebug( ",");
                    st_logDebug("%s", unpackedRecordRef->d.allele[i]);
                }
                st_logDebug("\n\toutput alleles: %s, %s\n", evalRefChar, evalAltChar);
            }

        } else {
            // Missed the variant
            // False negative - no variation was found, but truth vcf has one
            results->falseNegatives++;

            double *read1BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads1, referencePos);
            double *read2BaseCounts = getProfileSequenceBaseCompositionAtPosition(reads2, referencePos);

            // Print out info about the record and its composition
            st_logDebug("\npos: %" PRIi64 "\n\tref: %s\talt: ", referencePos, refChar);
            for (int i = 1; i < unpackedRecordRef->n_allele; i++) {
                if (i != 1) st_logDebug(",");
                st_logDebug("%s", unpackedRecordRef->d.allele[i]);
            }

            char h1AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString1[referencePos-gF->refStart]);
            char h2AlphChar = stBaseMapper_getCharForValue(baseMapper, gF->haplotypeString2[referencePos-gF->refStart]);
            st_logDebug("\n\toutput alleles: %c, %c\n", h1AlphChar, h2AlphChar);
            printColumnAtPosition(hmm, referencePos);


            // Check if record was an insertion or deletion
            if (strlen(refChar) > 1 || strlen(refAltChar) > 1) {
                results->error_missedIndels++;
                size_t indelLen = strlen(refChar);
                if (strlen(refAltChar) > indelLen) indelLen = strlen(refAltChar);

                st_logDebug("INSERTION / DELETION!!\n");
                for (int64_t j = 1; j < indelLen; j++) {
                    st_logDebug("\tNext pos: %" PRIi64 "\n", referencePos+j);
                    printColumnAtPosition(hmm, referencePos+j);
                }
            } else {
                // Quantify the type of false negative it was
                double totalCount = 0;
                for (int64_t i = 0; i < ALPHABET_SIZE; i++) {
                    totalCount += read1BaseCounts[i];
                    totalCount += read2BaseCounts[i];
                }
                int refBase = stBaseMapper_getValueForChar(baseMapper, *refChar);
                int altBase = stBaseMapper_getValueForChar(baseMapper, *refAltChar);
                float fractionRefBase = (read1BaseCounts[refBase] + read2BaseCounts[refBase]) / totalCount;
                float fractionAltBase = (read1BaseCounts[altBase] + read2BaseCounts[altBase]) / totalCount;


                if (fractionRefBase < threshold || fractionAltBase < threshold) {
                    results->error_trueVariantWrong++;
                    st_logDebug("   TRUE VARIANT WRONG\n");
                } else {
                    results->error_badPartition++;
                    st_logDebug("\tPartition 1: \n");
                    printBaseComposition2(read1BaseCounts);
                    st_logDebug("\tPartition 2: \n");
                    printBaseComposition2(read2BaseCounts);
                    st_logDebug( "*** BAD PARTITION\n");
                }
                st_logDebug("\tfraction of reference base seen in reads: %f\n", fractionRefBase);
                st_logDebug("\tfraction of alternate base seen in reads: %f\n", fractionAltBase);
            }
            st_logDebug("\tposterior prob: %f\n", gF->genotypeProbs[referencePos-gF->refStart]);

            free(read1BaseCounts);
            free(read2BaseCounts);
        }
    }

    // Remaining positions after the last variant in the reference are not currently being looked through
    // False positives in this region could therefore be missed
    // (In addition to false positives after the first variant)
    results->negatives = referencePos - refStart - results->positives;
    results->trueNegatives = results->negatives - results->falsePositives;


    // cleanup
    vcf_close(inRef);
    vcf_close(inEval);
    bcf_hdr_destroy(hdrRef);
    bcf_hdr_destroy(hdrEval);
    bcf_destroy(refRecord);
    bcf_destroy(evalRecord);
    stSet_destruct(reads1);
    stSet_destruct(reads2);
    stGenomeFragment_destruct(gF);
    stList_destruct(path);
}

/*
 * Performs test to report how well genotypes for the sequence are determined.
 * Reports statistics such as sensitivity and specificity.
 */
struct genotypeResults *genotypingTest(char *paramsFile, char *bamFile, char *vcfOutFile, char *vcfOutFileDiff, char *referenceFile, char *vcfReference, double threshold, int64_t iterationsOfParameterLearning) {

    st_logInfo("> Parsing parameters\n");
    stBaseMapper *baseMapper = stBaseMapper_construct();
    stRPHmmParameters *params = parseParameters(paramsFile, baseMapper);
    // Print a report of the parsed parameters
    if (st_getLogLevel() >= debug) stRPHmmParameters_printParameters(params, stderr);

    st_logInfo("> Creating profile sequences\n");
    stList *profileSequences = stList_construct3(0, (void (*)(void *))stProfileSeq_destruct);;
    parseReads(profileSequences, bamFile, baseMapper);

    // Print some stats about the input sequences
    printSequenceStats(profileSequences);
    printAvgIdentityBetweenProfileSequences(profileSequences, 100);

    st_logInfo("> Creating reference prior probabilities\n");
    stHash *referenceNamesToReferencePriors = createReferencePriorProbabilities(referenceFile, profileSequences, baseMapper, params);

    st_logInfo("> Learning parameters for HMM model (%" PRIi64 ")\n", iterationsOfParameterLearning);
    stRPHmmParameters_learnParameters(params, profileSequences, referenceNamesToReferencePriors, iterationsOfParameterLearning);

    st_logInfo("> Building hmms\n");
    stList *hmms = getRPHmms(profileSequences, referenceNamesToReferencePriors, params);
    stList *l = stList_construct3(0, (void (*)(void *))stRPHmm_destruct2);
    while(stList_length(hmms) > 0) {
        stList_appendAll(l, stRPHMM_splitWherePhasingIsUncertain(stList_pop(hmms)));
    }
    hmms = l;
    st_logInfo("\tGot %" PRIi64 " hmms\n", stList_length(hmms));

    st_logInfo("> Writing vcf files\n");
    vcfFile *vcfOutFP = vcf_open(vcfOutFile, "w");
    bcf_hdr_t *bcf_hdr = writeVcfHeader(vcfOutFP, l, referenceFile);

    vcfFile *vcfOutFP_diff = vcf_open(vcfOutFileDiff, "w");
    bcf_hdr_t *bcf_hdr_diff = writeVcfHeader(vcfOutFP_diff, l, referenceFile);

    struct genotypeResults *results = st_calloc(1, sizeof(struct genotypeResults));
    stGenomeFragment *gF;

    // Prep for BAM outputs
    stSet *read1Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);
    stSet *read2Ids = stSet_construct3(stHash_stringKey, stHash_stringEqualKey, NULL);

    for(int64_t i=0; i<stList_length(hmms); i++) {
        stRPHmm *hmm = stList_get(hmms, i);

        // Run the forward-backward algorithm
        stRPHmm_forwardBackward(hmm);

        // Now compute a high probability path through the hmm
        stList *path = stRPHmm_forwardTraceBack(hmm);

        // Compute the genome fragment
        gF = stGenomeFragment_construct(hmm, path);

        // Get the reads which mapped to each path
        stSet *reads1 = stRPHmm_partitionSequencesByStatePath(hmm, path, true);
        stSet *reads2 = stRPHmm_partitionSequencesByStatePath(hmm, path, false);


        addProfileSeqIdsToSet(reads1, read1Ids);
        addProfileSeqIdsToSet(reads2, read2Ids);


        if (st_getLogLevel() >= debug && stList_length(hmms) <= 5) {
            // Print stats about the HMM
            st_logDebug("Hmm # %" PRIi64 "\n", i);
            stRPHmm_print(hmm, stderr, 0, 0);


            st_logDebug("\nThere are %" PRIi64 " reads covered by the hmm, bipartitioned into sets of %" PRIi64 " and %" PRIi64 " reads\n",
                        stList_length(hmm->profileSeqs), stSet_size(reads1), stSet_size(reads2));

            // Print the similarity between the two imputed haplotypes sequences
            st_logDebug("The haplotypes have an %f identity\n", getIdentityBetweenHaplotypes(gF->haplotypeString1, gF->haplotypeString2, gF->length));

            // Print the base composition of the haplotype sequences
            double *hap1BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString1, gF->length);
            st_logDebug("The base composition of haplotype 1:\n");
            printBaseComposition(hap1BaseCounts);
            free(hap1BaseCounts);

            double *hap2BaseCounts = getHaplotypeBaseComposition(gF->haplotypeString2, gF->length);
            st_logDebug("The base composition of haplotype 2:\n");
            printBaseComposition(hap2BaseCounts);
            free(hap2BaseCounts);

            // Print the base composition of the reads
            double *reads1BaseCounts =getExpectedProfileSequenceBaseComposition(reads1);
            st_logDebug("The base composition of reads1 set:\n");
            printBaseComposition(reads1BaseCounts);
            free(reads1BaseCounts);

            double *reads2BaseCounts =getExpectedProfileSequenceBaseComposition(reads2);
            st_logDebug("The base composition of reads2 set:\n");
            printBaseComposition(reads2BaseCounts);
            free(reads2BaseCounts);

            // Print some summary stats about the differences between haplotype sequences and the bipartitioned reads
            st_logDebug("hap1 vs. reads1 identity: %f\n",
                        getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads1));
            st_logDebug("hap1 vs. reads2 identity: %f\n",
                        getExpectedIdentity(gF->haplotypeString1, gF->refStart, gF->length, reads2));

            st_logDebug("hap2 vs. reads2 identity: %f\n",
                        getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads2));
            st_logDebug("hap2 vs. reads1 identity: %f\n",
                        getExpectedIdentity(gF->haplotypeString2, gF->refStart, gF->length, reads1));

            stSet_destruct(reads1);
            stSet_destruct(reads2);

        }


        // Write vcf
        writeVcfFragment(vcfOutFP, bcf_hdr, gF, referenceFile, baseMapper, false);
        writeVcfFragment(vcfOutFP_diff, bcf_hdr_diff, gF, referenceFile, baseMapper, true);

        stList_destruct(path);
    }

    vcf_close(vcfOutFP);
    vcf_close(vcfOutFP_diff);
    compareVCFs(stderr, hmms, vcfOutFileDiff, vcfReference, threshold,
                baseMapper, results);

    char *samOutBase = "test";
    st_logInfo("\n> Writing out SAM files for each partition into files: %s.1.sam and %s.1.sam\n", samOutBase, samOutBase);
    writeSplitSams(bamFile, samOutBase, read1Ids, read2Ids);


    // cleanup
    stRPHmmParameters_destruct(params);
    stBaseMapper_destruct(baseMapper);
    stList_destruct(hmms);
    stList_destruct(profileSequences);
    stGenomeFragment_destruct(gF);
    stSet_destruct(read1Ids);
    stSet_destruct(read2Ids);

    bcf_hdr_destroy(bcf_hdr);
    bcf_hdr_destroy(bcf_hdr_diff);

    return results;
}

int genotypingTest2(FILE *fh, char *paramsFile, char *bamFile, char *vcfOutFile,
        char *referenceFile, char *vcfReference, int64_t refStart, int64_t refEnd,
        bool printVerbose, int64_t iterationsOfParameterLearning) {

    // Run margin phase
    char *logString = printVerbose ? "--logLevel DEBUG" : "";
    char *command = stString_print("./marginPhase %s %s %s --params %s --vcfFile %s "
            "--iterationsOfParameterLearning %" PRIi64 "",
            bamFile, referenceFile, logString,
            paramsFile, vcfOutFile, iterationsOfParameterLearning);
    fprintf(fh, "> Running margin phase with command: %s\n", command);
    return st_system(command);

    // TODO : Do VCF comparison using VCF eval
}

void test_5kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_5kb.vcf";
    char *vcfOutFileDiff = "test_5kb_diff.vcf";

    char *bamFile = "../tests/NA12878.pb.chr3.5kb.bam";
    // TODO: create vcf specifically for this 5 kb region

    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    double falseNegativeThreshold = 0.20;


    st_logInfo("Testing haplotype inference on %s\n", bamFile);
    genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, 3);

    //int i = genotypingTest2(stderr, paramsFile, bamFile, vcfOutFile, referenceFile, vcfReference, 190000, 195000, 1, 3);

    //CuAssertTrue(testCase, i == 0);

}

void test_100kbGenotyping(CuTest *testCase) {

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    double falseNegativeThreshold = 0.20;
    int64_t iterationsOfParameterLearning = 1;

    char *bamFile = "../tests/NA12878.pb.chr3.100kb.2.bam";
    //char *bamFile = "../tests/NA12878.ihs.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    //char *vcfReference = "../tests/HG001.GRCh37.chr3.100kb.vcf";

    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile,
                                                     vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);
    printGenotypeResults(results);

    free(results);

}

void test_multiple100kbGenotyping(CuTest *testCase) {

    st_setLogLevelFromString("info");

    char *paramsFile = "../tests/params.json";
    char *referenceFile = "../tests/hg19.chr3.9mb.fa";
    char *vcfOutFile = "test_100kb.vcf";
    char *vcfOutFileDiff = "test_100kb_diff.vcf";
    double falseNegativeThreshold = 0.20;
    int64_t iterationsOfParameterLearning = 0;

    char *bamFile = "../tests/NA12878.ihs.chr3.100kb.0.bam";
    char *vcfReference = "../tests/NA12878.PG.chr3.100kb.0.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results1 = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.1.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.1.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results2 = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.2.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.2.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results3 = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.3.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.3.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results4 = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);

    bamFile = "../tests/NA12878.ihs.chr3.100kb.4.bam";
    vcfReference = "../tests/NA12878.PG.chr3.100kb.4.vcf";
    st_logInfo("Testing haplotype inference on %s\n", bamFile);

    struct genotypeResults *results5 = genotypingTest(paramsFile, bamFile, vcfOutFile, vcfOutFileDiff, referenceFile, vcfReference, falseNegativeThreshold, iterationsOfParameterLearning);

    st_logInfo("***** Test 1 *****\n");
    printGenotypeResults(results1);
    st_logInfo("***** Test 2 *****\n");
    printGenotypeResults(results2);
    st_logInfo("***** Test 3 *****\n");
    printGenotypeResults(results3);
    st_logInfo("***** Test 4 *****\n");
    printGenotypeResults(results4);
    st_logInfo("***** Test 5 *****\n");
    printGenotypeResults(results5);

    free(results1);
    free(results2);
    free(results3);
    free(results4);
    free(results5);
}



CuSuite *marginPhaseTestSuite(void) {

    CuSuite* suite = CuSuiteNew();

    st_setLogLevelFromString("debug");
    //st_setLogLevelFromString("info");

//    SUITE_ADD_TEST(suite, test_5kbGenotyping);
    SUITE_ADD_TEST(suite, test_100kbGenotyping);
//    SUITE_ADD_TEST(suite, test_multiple100kbGenotyping);

    return suite;
}
