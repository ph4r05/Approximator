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
#define MAX_ORDER 6

//
// FGb
//
#include "faugere/fgb.h"
// NTL library.
#include <NTL/vec_vec_GF2.h>

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
    // Indexing is as follows coefficients[order][coefficient_index][polyout].
    //    Where order is the order of terms represented in the structure.
    //    coefficient_index is the order number of the combination of variables
    //      in particular term, terms follow lexicographic order.
    //    polyout: Output block for polynomials. Each bit in this block
    //      tells whether in the particular polynomial given term is present or not.
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
    std::vector<ULONG> coefficients[MAX_ORDER];
    
    // Limit on the term order for storage.
    // Only terms of order/degree less than or equals to this limit will
    // be stored and evaluated.
    uint orderLimit;
    
    // Flag telling whether to dump coefficients to the file or not.
    bool dumpCoefsToFile;
    
    // Byte width of the input cipher.
    // Key block size + message block size.
    ULONG byteWidth;
    
    // Number of bits needed to store the 
    uint logBitInputWidth;
    
    // Size of the output block of the cipher in ulong types.
    ULONG outputWidthUlong;
    
    // Size of the input block (message+key) of the cipher in ulong types.
    ULONG inputWidthUlong;
    
    // Array of binomial sums for computing combination index in the lexicographic ordering.
    // The first array element corresponds to the order 1, next to the order 2 and so on...
    ULONG ** binomialSums;
    
    // Number of threads to use for parallelized computation.
    uint threadCount;
    
    // Variable names for FGb.
    char ** varNames;
    
    // Number of key-bits set to zero.
    uint keybitsToZero;
    
    // Bitmap of polynomials to take into consideration (for example during GB solving).
    ULONG * poly2take;
    // Hamming weight of the poly2take.
    uint numPolyActive;
    // Log file for FGb library.
    FILE * fgbFile;
    
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
     * Determines combination from its index.
     * O(order*numvariables)
     * 
     * @param input
     * @param output
     * @param iBuff
     * @param oBuff
     * @return 
     */
    int getCombinationFromIdx(uint order, ULONG * xs, ULONG idx);
    
    /**
     * Computes combination ULONG number. 
     * Generated numbers are not continuous (not space effective)!
     * Format:
     *  4bits to store order (number of elements to choose)
     *  1. element
     *  2. element
     *  ...
     * 
     * Benefit: Easily reversible.
     * 
     * @param order
     * @param xs
     * @return 
     */
    ULONG getCombinationULong(uint order, const ULONG * xs);
    
    /**
     * Determines combination from its ULONG index.
     * 
     * @param xs    Buffer the combination will be stored in. Has to be big enough to store the whole combination.
     * @param comb  Combination ULONG number.
     * @return 
     */
    int getCombinationFromULong(ULONG * xs, ULONG combUlong);
    
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
     * Partially evaluates approximated function on a given input.
     * Result is again a function, but with smaller number of variables. 
     * 
     * @param newVariables          Number of the variables in a new partially evaluated function.
     * @param variablesValueMask    Variables bit mask for evaluation (we have values for them.) 
     *                              Hamming_weight(variablesValueMask) = bitWIdth - newVariables.
     * @param iBuff                 Input to the evaluation. Values for variables to evaluate.
     *                              Only masked values are taken into consideration.
     * @param coeffEval             Storage provided by the caller to store the new function in.
     * @return 
     */
    int partialEvaluation(uint numVariables, ULONG * variablesValueMask, ULONG * iBuff, std::vector<ULONG> * coeffEval);
    
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
    void solveKeyGrobner(uint samples, bool dumpInputBase=false);
    
    void solveGb(uint numVariables, Dpol * basis, uint numPoly);
    
    /**
     * Generated FGb polynomial representation.
     * Allocates a new memory.
     * @param coefs     Coefficient storage for the polynomials.
     * @param maxOrder  Maximal order of the terms stored in coefs.
     * @param polyIdx   Which polynomial to represent.
     * @param numTerms [OPTIONAL] If non-null, it will contain number of terms in the polynomial.
     */
    Dpol_INT polynomial2FGb(uint numVariables, std::vector<ULONG> * coefs, uint maxOrder, uint polyIdx, ULONG * numTerms);
    
    /**
     * Dumps FGb polynomial to the standard output.
     * @param numVariables
     * @param poly
     */
    void dumpFGbPoly(uint numVariables, Dpol poly);
    
    /**
     * Dumps polynomial basis.
     */
    void dumpBasis(uint numVariables, Dpol * basis, uint numPoly);
    
    /**
     * Initializes FGb library.
     * @param numVariables
     */
    void initFGb(uint numVariables);
    
    /**
     * Deinitializes FGb library.
     */
    void deinitFGb();
    
    /**
     * Reads 8-bit buffer to the 64 bit buffer.
     * iBuff has to be big enough to fit the input buffer.
     * 
     * @param input     input buffer to read.
     * @param size      size of the input buffer in bytes to read.
     * @param iBuff     destination buffer.
     */
    void readUcharToUlong(const uchar * input, uint size, ULONG * iBuff) const;
    
    /**
     * Reads 64 bit buffer to the 8 bit buffer.
     * 
     * @param output    output buffer to write.
     * @param size      size of the input buffer to read in bytes.
     * @param iBuff     input buffer to read (from LSB).
     */
    void readUlongToUchar(uchar * output, uint size, const ULONG * iBuff) const;
    
    ICipher * getCipher() const { return cip; }
    void      setCipher(ICipher * cip);
    
    uint getOrderLimit() const { return orderLimit; }
    
    uint getThreadCount() const { return threadCount; }
    void setThreadCount(uint threadCount) { this->threadCount = threadCount; }

    uint getKeybitsToZero() const { return keybitsToZero; }
    void setKeybitsToZero(uint keybitsToZero) { this->keybitsToZero = keybitsToZero; }

    void setPoly2Take(const std::vector<std::string> & map);
    bool isPoly2Take(uint polyIdx) const;
    
    void genMessages();
};

#endif	/* APPROXIMATION_H */

