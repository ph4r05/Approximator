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
    
    // Array of binomial sums for computing combination index in the lexicographic ordering.
    // The first array element corresponds to the order 3, next to the order 4 and so on...
    ULONG ** binomialSums;
    
    // Number of threads to use for parallelized computation.
    uint threadCount;
    
public:
    Approximation(uint orderLimit);
    virtual ~Approximation();
    
    /**
     * Entry point.
     */
    void work();
    
    /**
     * Computes coefficients for polynomial approximation.
     * Init has to be called before this function.
     * Memory for coefficient is allocated in this step.
     */
    void computeCoefficients();
    
    /**
     * Initialization of the internal pre-computed tables.
     */
    void init();
    
    /**
     * Returns the number of particular combination, assuming N=cipher input width.
     * Uses precomputed values to optimize computation - cubeBinomialSums. 
     * @param ULONG
     */
    ULONG getCubeIdx(ULONG x1, ULONG x2, ULONG x3) const;
    
    /**
     * General function for computing an combination index in an lexicographic ordering.
     * Works only up to defined degree.
     * 
     * @param order       number of elements in the combination.
     * @param xs          ULONG representation of the combination (output of the comb. generator).
     * @param xsOffset    offset for the combination in xs to take into account (for recursive computation).
     * @param Noffset     minus offset to the N (for recursion).
     * @param combOffset  minus offset for combination elements (recursion).
     * @return 
     */
    ULONG getCombinationIdx(uint order, const ULONG * xs, uint xsOffset=0, ULONG Noffset=0, ULONG combOffset=0) const;
    
    /**
     * Evaluates function determined by coefficients 
     * @param input
     * @param output
     * @param iBuff  input  ULONG buffer to use (has to be allocated to exact size)
     * @param oBuff  output ULONG buffer to use (has to be allocated to exact size)
     * @return 
     */
    int evaluateCoefficients(const unsigned char * input, unsigned char * output, ULONG * iBuff, ULONG * oBuff) const;
    
    /**
     * Tests polynomial approximation of the cipher.
     * Pre-computed coefficients are used.
     * @return 
     */
    int testPolynomialApproximation(unsigned long numSamples);
    
    /**
     * Performs self test on determined coefficient values.
     * For low degree input we have to obtain exactly the same results from 
     * cipher and from the approximation function.
     * 
     * @return 
     */
    int selftestApproximation(unsigned long numSamples);
    
    /**
     * Simple test for combination indexing correctness.
     */
    int selftestIndexing();
    
    /**
     * Computes maximal number of terms the polynomial with given number
     * of variables and with the defined maximal order can have.
     * @return 
     */
    ULONG numberOfTerms(ULONG variables, ULONG maxOrder);
    
    /**
     Procedure for solving equation for keys with using GB.
     */
    void solveKeyGrobner(uint samples);
    
    ICipher * getCipher() const { return cip; }
    void      setCipher(ICipher * cip);
    
    uint getOrderLimit() const { return orderLimit; }
    
    uint getThreadCount() const { return threadCount; }
    void setThreadCount(uint threadCount) { this->threadCount = threadCount; }


    void genMessages();
};

#endif	/* APPROXIMATION_H */

