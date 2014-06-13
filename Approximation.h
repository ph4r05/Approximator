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

// Maximal order used during cube attack.
#define MAX_SUPERPOLY_ORDER 2

//
// FGb
//
#include "faugere/fgb.h"
#include "CombinatorialIndexer.h"
#include "FGbHelper.h"
// NTL library.
#include <NTL/vec_vec_GF2.h>
#include <NTL/vec_GF2.h>

class Approximation {
private:
    ICipher  * cip;
    
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
    
    // Combinatorial indexer - helps with computing index of combinations and vice versa.
    CombinatorialIndexer combIndexer;
    
    // FGb helper class.
    FGbHelper fgb;
    
    // Number of threads to use for parallelized computation.
    uint threadCount;
    
    // Number of key-bits set to zero.
    uint keybitsToZero;
    
    // Bitmap of polynomials to take into consideration (for example during GB solving).
    ULONG * poly2take;
    // Hamming weight of the poly2take.
    uint numPolyActive;
    
public:
    Approximation(uint orderLimit);
    virtual ~Approximation();
    
    /**
     * Computes coefficients for polynomial approximation.
     * Init has to be called before this function.
     * Memory for coefficient is allocated in this step.
     */
    void computeCoefficients(std::vector<ULONG> * coefficients);
    
    /**
     * Initialization of the internal pre-computed tables.
     */
    void init();
    
    /**
     * Evaluates function determined by coefficients 
     * @param input
     * @param output
     * @param iBuff  input  ULONG buffer to use (has to be allocated to exact size)
     * @param oBuff  output ULONG buffer to use (has to be allocated to exact size)
     * @return 
     */
    int evaluateCoefficients(const std::vector<ULONG> * coefficients,
        const unsigned char * input, unsigned char * output, ULONG * iBuff, ULONG * oBuff) const;
    
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
    int partialEvaluation(const std::vector<ULONG> * coefficients,
        uint numVariables, ULONG * variablesValueMask, ULONG * iBuff, std::vector<ULONG> * coeffEval) const;
    
    /**
     * Tests polynomial approximation of the cipher.
     * Pre-computed coefficients are used.
     * @return 
     */
    int testPolynomialApproximation(unsigned long numSamples) const;
    
    /**
     * Performs self test on determined coefficient values.
     * For low degree input we have to obtain exactly the same results from 
     * cipher and from the approximation function.
     * 
     * @return 
     */
    int selftestApproximation(unsigned long numSamples) const;
    
    /**
     * Simple test for combination indexing correctness.
     */
    int selftestIndexing() const;
    
    /**
     * Computes maximal number of terms the polynomial with given number
     * of variables and with the defined maximal order can have.
     * @return 
     */
    ULONG numberOfTerms(ULONG variables, ULONG maxOrder) const;
    
    /**
     Procedure for solving equation for keys with using GB.
     * FGb has to be initialized before and deinitialized after calling this method!
     */
    void solveKeyGrobner(uint samples, bool dumpInputBase=false, bool selfTest=false, int basisReduction=0) const;
    
    int solveGb(uint numVariables, Dpol* basis, uint numPoly, uchar * solvedKey) const;
    
    /**
     * Function computes subCube of a given term specified by termWeight and 
     * term bitmask (with 1 on positions corresponding to a variable present in 
     * term to cube).
     * 
     * Finput is function input that will be used for cube process in target
     * function evaluation, except bits specified in termMask, those will be
     * always set accordingly to the cubing. Term bits has to be set to 0!
     * 
     * Function allows parallelization, does not modify any internal state, 
     * can be configured to compute only a sub-cube (each x-th term in the cube).
     * 
     * Result is XOR of all possible terms constructible from given term.
     * 
     * Result is written to the subcube argument. 
     * 
     * @param termWeight Term to cube, number of variables, Hamming weight of its bitmask.
     * @param termMask   Term to cube, bitmask
     * @param finput     Input for the function evaluation (contains keys).
     * @param subcube    Output parameter.
     * @param step       Parallelization. (i*step) + offset
     * @param offset     Parallelization. (i*step) + offset
     * @param incldueTerm If true, also the term specified in the arguments will
     *                   be part of the result.
     * @return 
     */
    int subCubeTerm(uint termWeight, const ULONG * termMask, const uchar * finput, ULONG * subcube,
        uint step, uint offset, bool includeTerm, bool precompKey) const;
    
    /**
     * Threaded variant of subCubeTerm.
     * Difference: Key is assumed to be constant, saved in finput...
     * 
     * @param termWeight
     * @param termMask
     * @param finput
     * @param subcube
     * @param step
     * @param offset
     * @param includeTerm
     * @return 
     */
    int subCubeTermThreaded(uint termWeight, const ULONG * termMask, const uchar * finput, ULONG * subcube, bool includeTerm) const;
    
    int cubeAttack(uint wPlain, uint wKey, uint numRelations) const;
    
    /**
     * Initializes FGb library.
     * @param numVariables
     */
    void initFGb(uint numVariables) const;
    
    /**
     * Deinitializes FGb library.
     */
    void deinitFGb() const;
    
    /**
     * Reset FGb internal memory.
     * @param output
     * @param size
     * @param iBuff
     */
    void resetFGb() const;
    
    /**
     * Returns number of variables for current cipher and key-bits-to-zero setting.
     * @return 
     */
    uint getNumVariables() const;
    
    std::vector<ULONG> * getCoefficients() { return coefficients; }

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

typedef struct CubeRelations_t_ {
    // Public variables (plaintext bits)
    std::vector<ULONG> termMask;
    // Bitmask of output polynomials having superpoly in this polynomial.
    std::vector<ULONG> isSuperpoly;
    // Hamming weight of the vector above.
    uint numSuperpolys;
    // Superpolys for each output polynomial (vectorized representation).
    std::vector<ULONG> superpolys[MAX_SUPERPOLY_ORDER];
} CubeRelations_t;

#endif	/* APPROXIMATION_H */

