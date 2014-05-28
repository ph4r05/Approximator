/* 
 * File:   FGbHelper.h
 * Author: ph4r05
 *
 * Created on May 28, 2014, 11:51 AM
 */

#ifndef FGBHELPER_H
#define	FGBHELPER_H

#include "base.h"
#include <string>
#include "faugere/fgb.h"

class FGbHelper {
private:
    // Variable names for FGb.
    char ** varNames;
    
    // Log file for FGb library.
    FILE * fgbFile;
    
    // Limit on the term order for storage.
    // Only terms of order/degree less than or equals to this limit will
    // be stored and evaluated.
    uint orderLimit;
    
    // Byte width of the input cipher.
    // Key block size + message block size.
    ULONG byteWidth;
    
    // Size of the output block of the cipher in ulong types.
    ULONG outputWidthUlong;   
public:
    FGbHelper();
    virtual ~FGbHelper();
    
    /**
     * Initialize helper object.
     * 
     * @param byteWidth     Function input byte width.
     * @param orderLimit    Order limit on terms.
     * @param outputBits    Number of bits on output function.
     */
    void init(ULONG byteWidth, uint orderLimit, uint outputBits);
    
    /**
     * Deinitialize internal state.
     */
    void deinit();
    
    /**
     * Returns number of terms in given polynomial.
     * @param poly
     * @return 
     */
    I32 getNumberOfTerms(Dpol poly) const;
    
    /**
     * Determines if polynomial is 0.
     * @param poly
     * @return 
     */
    bool isPolyNull(Dpol poly) const;
    
    /**
     * Determines if polynomial is 1*1.
     * @param poly
     * @return 
     */
    bool isPoly1(Dpol poly, I32 numVariables) const;
    
    /**
     * Calls underlying FGb routine to export polynomial to the 
     * 
     * @param numVariables  Number of variables the polynomial may contain.
     * @param numTerms      Number of terms in the polynomial, use getNumberOfTerms().
     * @param exps          Exponents for each term and each variable in the term.
     *                      Size is numVariables*numTerms. Terms, then variables.
     * @param coefs         Array of size numTerms, stores coefficients for terms.
     * @param polynomial
     * @return 
     */
    I32 exportPolynomial(I32 numVariables, I32 numTerms, I32* exps, I32 * coefs, Dpol polynomial) const;
    
    /**
     * Generated FGb polynomial representation.
     * Allocates a new memory.
     * @param coefs     Coefficient storage for the polynomials.
     * @param maxOrder  Maximal order of the terms stored in coefs.
     * @param polyIdx   Which polynomial to represent.
     * @param numTerms [OPTIONAL] If non-null, it will contain number of terms in the polynomial.
     * @param hash [OPTIONAL] If non-null, hash of the polynomial will be computed and set here.
     */
    Dpol_INT polynomial2FGb(uint numVariables, std::vector<ULONG> * coefs, uint maxOrder, uint polyIdx, ULONG * numTerms = NULL, ULONG * hash = NULL) const;
    
    /**
     * Dumps FGb polynomial to the standard output.
     * @param numVariables
     * @param poly
     */
    void dumpFGbPoly(uint numVariables, Dpol poly) const;
    
    /**
     * Dumps polynomial basis.
     */
    void dumpBasis(uint numVariables, Dpol * basis, uint numPoly) const;
    
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
     * Computes Gb with standard settings.
     * @param n_input
     * @param inputBasis
     * @param outputBasis
     * @param t0
     * @return 
     */
    int computeFGb(int n_input, Dpol * inputBasis, Dpol * outputBasis, double * t0) const;
    
private:

};

#endif	/* FGBHELPER_H */

