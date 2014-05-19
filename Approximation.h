/* 
 * File:   Approximation.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 10:21 AM
 */

#ifndef APPROXIMATION_H
#define	APPROXIMATION_H

#include <iostream>
#include <stdio.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstring>

#include "base.h"
#include "ICipher.h"

// Maximum order of terms.
#define MAX_ORDER 7

class Approximation {
private:
    ICipher  * cip;
    
    // Key & input parameters prepared for cipher input.
    uchar * finput;
    
    // Input & output blocks defined on ulongs to speed up computation.
    // During approximated function evaluation these are used as buffers.
    ULONG * ulongInp;
    ULONG * ulongOut;
    
    // Cache polynomial coefficients for low order.
    // Indexing is as follows coefficients[order][coefficient index][polyout].
    //
    // Size of this array can be precomputed from the cipher.
    // Let iw = input width of the cipher in bytes (message+key).
    // Let ow = output width of the cipher in bytes.
    //
    // Then second dimension is number of all possible terms with specified order ord.
    // In particular it is Binomial(8*iw, ord).
    // Third dimension is ceil(ow/sizeof(ulong)).
    // In order to optimize memory storage/access of/to this structure
    // 2nd and 3rd dimension are merged to one.
    std::vector<ULONG> coefficients[5];
    
    // Limit on the term order for storage.
    // Only terms of order/degree less than or equals to this limit will
    // be stored and evaluated.
    uint orderLimit;
    
    // Flag telling whether to dump coefficients to the file or not.
    bool dumpCoefsToFile;
    
    // Byte width of the input cipher.
    // Key block size + message block size.
    ULONG byteWidth;
    
    // Size of the output block of the cipher in ulong types.
    ULONG outputWidthUlong;
    
    // Size of the input block (message+key) of the cipher in ulong types.
    ULONG inputWidthUlong;
    
    // Binomial sums for computing coefficient indexes for x1x2x3.
    ULONG * cubeBinomialSums;
    
public:
    Approximation();
    virtual ~Approximation();
    
    /**
     * Entry point.
     */
    void work();
    
    /**
     * Initialization of the internal pre-computed tables.
     */
    void init();
    
    /**
     * Returns the number of particular combination, assuming N=cipher input width.
     * Uses precomputed values to optimize computation - cubeBinomialSums. 
     * @param ULONG
     */
    ULONG getCubeIdx(ULONG x1, ULONG x2, ULONG x3);
    
    /**
     * Evaluates function determined by coefficients 
     * @param input
     * @param key
     * @param output
     * @return 
     */
    int evaluateCoefficients(const unsigned char * input, unsigned char * output);
    
    /**
     * Tests polynomial approximation of the cipher.
     * Pre-computed coefficients are used.
     * @return 
     */
    int testPolynomialApproximation();
    
    ICipher * getCipher(){ return cip; }
    void      setCipher(ICipher * cip);
    
    void genMessages();
};

#endif	/* APPROXIMATION_H */

