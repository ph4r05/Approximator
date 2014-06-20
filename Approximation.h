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
#define MAX_SUPERPOLY_ORDER 4

//
// FGb
//
#include "faugere/fgb.h"
#include "CombinatorialIndexer.h"
#include "FGbHelper.h"
// NTL library.
#include <NTL/vec_vec_GF2.h>
#include <NTL/vec_GF2.h>
// Boost serialization
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>

class CubeRelations_t;
class CubeRelations_vector;

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
    
    // Signal blocking variables.
    sigset_t pendingSignals, blockingMask;
    
    // Verbosivity level
    uint verboseLvl;
    
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
     * Function
     * 
     * @param termWeight Term to cube, number of variables, Hamming weight of its bitmask.
     * @param termMask   Term to cube, bitmask
     * @param finput     Input for the function evaluation (contains keys).
     * @param subcube    Output buffer.
     * @param step       Parallelization. (i*step) + offset
     * @param offset     Parallelization. (i*step) + offset
     * @param doSubCubes If true, also the cubes of size N-1 are computed.
     * @return 
     */
    int subCubeTerm(uint termWeight, const ULONG * termMask, const uchar * finput, ULONG * subcube,
        uint step, uint offset, uint subCubes, bool precompKey) const;
    
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
     * @param doSubCubes
     * @return 
     */
    int subCubeTermThreaded(uint termWeight, const ULONG * termMask, const uchar * finput, ULONG * subcube, uint subCubes) const;
    
    /**
     * Returns cache file name for the cube attack with specified settings.
     * @param wPlain
     * @param wKey
     * @return 
     */
    std::string getCubeCacheName(uint wPlain, uint wKey) const;
    
    /**
     * Loads cube attack archive to the cube relations vector.
     * @param fname
     * @param 
     * @return 
     */
    int readCubeArchive(const char * fname, CubeRelations_vector & vct) const;
    
    /**
     * Writes complete cube relations vector to the given file. 
     * Blocks sigterm, sigquit, sigint signals to avoid archive corruption 
     * during save procedure. 
     * 
     * @param fname
     * @param vct
     * @return 
     */
    int writeCubeArchive(const char * fname, CubeRelations_vector & vct) const;
    
    /**
     * Write given relation to the archive.
     * 
     * @param fname
     * @param vct
     * @param wKey
     * @param termMask
     * @param keyCubes
     * @param isSuperpoly
     * @return 
     */
    int writeRelationToArchive(const char * fname, CubeRelations_vector & vct,
        uint wKey, ULONG * termMask, std::vector<ULONG> * keyCubes, ULONG * isSuperpoly) const;
    
    /**
     * Dumps function defined by the order vectors to the stream in plaintext.
     * 
     * @param c
     * @param coefficients
     * @param maxOrder
     */
    ULONG dumpCoefficients(std::ostream & c, const std::vector<ULONG> * coefficients, uint maxOrder, uint numVariables, uint polyIdx, uint fmt=1) const;
    
    /**
     * Dumps function representation of a multiple function in a polynomial form.
     * 
     * @param c             Output stream to write data to.
     * @param coefficients  Storage structure for polynomial representation.
     * @param maxOrder      Maximum order of a polynomial terms to dump.
     * @param numVariables  Number of variables in one polynomial (describes structure).
     * @param numPoly       Number of polynomials to dump.
     * @param nonNullOnly   if true only non-null polynomials will be dumped.
     * @param fmt   formatting: 1=polynomial textual, 2=ASCII binary, 3=binary
     * @return 
     */
    ULONG dumpOutputFunctions(std::ostream & c, const std::vector<ULONG> * coefficients, uint maxOrder, 
        uint numVariables, uint numPoly, bool nonNullOnly=true, uint fmt=1) const;
    
    /**
     * Generate array of bit positions from left to right. 
     * @param termBitPositions
     * @param bitWidth
     * @param termMask
     * @return 
     */
    uint genBitPosMask(uint * termBitPositions, uint bitWidth, const ULONG* termMask, uint termWeight) const;
    
    /**
     * On provided plaintext part performs cube computation on key variables.
     * Plaintext remains fixed during the computation (stored in termMask),
     * key cube starts from order <b>startOrder</b> and stops in order
     * <b>stopOrder</b> inclusively. 
     * 
     * keyCubes and isSuperpoly has to be already initialized. No memory reset
     * is performed so this can be called to finish a computation, e.g., 
     * to compute quadratic cube on already computed linear cube.
     * Results are stored to these variables and these variables are used 
     * for computation (e.g., for computing quadratic key cube, results for
     * linear key cube are used). 
     * 
     * isSuperpoly is standard ULONG array of size outputWidthUlong,
     * keyCubes is array of vectors, each vector for each order of computation.
     * 
     * No sub-cubes are computed (no optimization).
     * 
     * @param wPlain
     * @param wKey
     * @param startOrder
     * @param stopOrder
     * @param termMask
     * @param keyCubes
     * @param isSuperpoly
     * @return 
     */
    int keyCube(uint wPlain, uint wKey, uint startOrder, uint stopOrder,
        ULONG * termMask, std::vector<ULONG> * keyCubes, ULONG * isSuperpoly) const;
    
    /**
     * Cube attack.
     * 
     * @param wPlain
     * @param wKey
     * @param numRelations
     * @param subCubesLimit
     * @param saveRelations
     * @param dumpAllRelations  If true all relations (also null ones) will be dumped.
                                If true, cube attack will not work, just for randomness test of the 
                                coefficient distribution.
     * @return 
     */
    int cubeAttack(uint wPlain, uint wKey, uint numRelations, uint subCubesLimit, 
        bool saveRelations=true, bool dumpAllRelations=false) const;
    
    /**
     * Online phase of the cube attack. Key is fixed, goal is to determine it.
     * We have to determine b_t for each relation:
     *
     * \Sum_{v \in C_t} f(v,x) = b_t
     * a_1x_1 + a_2x_2 + \dots + a_nx_n + c = b_t       (this is for one relation)
     * 
     * From this we get system of n variables and more than n equations we want 
     * to solve with GB or Gaussian elimination.
     * 
     * @param keyRelationsVector
     * @param input Input block with prepared MAC key.
     * @param solution Recovered key will be stored here
     * @return 
     */
    long cubeOnlineAttack(CubeRelations_vector & keyRelationsVector, const uchar * input, uchar * solvedKey) const;
    
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
    
    uint getVerboseLvl() const { return verboseLvl; }
    void setVerboseLvl(uint verboseLvl) { this->verboseLvl = verboseLvl; }
    
    void genMessages();
};

class CubeRelations_t {
public:
    // Public variables (plaintext bits)
    std::vector<ULONG> termMask;
    // Bitmask of output polynomials having superpoly in this polynomial.
    std::vector<ULONG> isSuperpoly;
    // Hamming weight of the vector above.
    uint numSuperpolys;
    // Size of the key cube.
    uint wkey;
    // Superpolys for each output polynomial (vectorized representation).
    std::vector<ULONG> superpolys[MAX_SUPERPOLY_ORDER];
    
    CubeRelations_t() {}
    virtual ~CubeRelations_t() {}
    
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(termMask);
        ar & BOOST_SERIALIZATION_NVP(isSuperpoly);
        ar & BOOST_SERIALIZATION_NVP(numSuperpolys);
        ar & BOOST_SERIALIZATION_NVP(wkey);
        ar & BOOST_SERIALIZATION_NVP(superpolys);
    }
};

class CubeRelations_vector {
public:
    // Public variables (plaintext bits).
    std::vector<CubeRelations_t> stor;
    // Total number of relations stored (for all polynomials together).
    uint totalRelations;
    
    CubeRelations_vector() : totalRelations(0) {}
    virtual ~CubeRelations_vector() {}
    
    std::vector<CubeRelations_t> & get() { return stor; }
    uint getTotal() const { return totalRelations;    }
    void setTotal(uint t) { this->totalRelations = t; }
    uint * getTotalPtr()  { return &totalRelations; }
    
    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(totalRelations);
        ar & BOOST_SERIALIZATION_NVP(stor);
    }
};

#endif	/* APPROXIMATION_H */

