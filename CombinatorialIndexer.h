/* 
 * File:   CombinatorialIndexer.h
 * Author: ph4r05
 *
 * Created on May 28, 2014, 11:08 AM
 */

#ifndef COMBINATORIALINDEXER_H
#define	COMBINATORIALINDEXER_H

#include "base.h"
class CombinatorialIndexer {
private:
    // Limit on the term order for storage.
    // Only terms of order/degree less than or equals to this limit will
    // be stored and evaluated.
    uint down;    
    
    // Byte width of the input cipher.
    // Key block size + message block size.
    uint up;
    
    // Number of bits needed to store the 
    uint logBitInputWidth;
    
    // Array of binomial sums for computing combination index in the lexicographic ordering.
    // The first array element corresponds to the order 1, next to the order 2 and so on...
    ULONG ** binomialSums;
    
public:
    CombinatorialIndexer();
    virtual ~CombinatorialIndexer();
    
    /**
     * Initialization of the internal pre-computed tables.
     * @param up = number of total elements available to choose.
     * @param down = size of a set of element to be chosen.
     */
    void init(uint up, uint down);
    
    /**
     * Frees internal memory, after de-init this class has to be initialized again
     * before normal use.
     */
    void deinit();
    
    /**
     * Returns the number of particular combination, assuming N=cipher input width.
     * Uses precomputed values to optimize computation - cubeBinomialSums. 
     * @param ULONG
     */
    inline ULONG getCubeIdx(ULONG x1, ULONG x2, ULONG x3) const;
    
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
    int getCombinationFromIdx(uint order, ULONG * xs, ULONG idx) const;
    
    /**
     * Computes combination ULONG number. 
     * Generated numbers are not continuous (not space effective)!
     * Format:
     *  4bits to store order (number of elements to choose)
     *  1. element
     *  2. element
     *  ...
     * 
     * Warning: This format can store only small order!
     * Benefit: Easily reversible.
     * 
     * @param order
     * @param xs
     * @return 
     */
    ULONG getCombinationULong(uint order, const ULONG * xs) const;
    
    /**
     * Determines combination from its ULONG index.
     * 
     * @param xs    Buffer the combination will be stored in. Has to be big enough to store the whole combination.
     * @param comb  Combination ULONG number.
     * @return 
     */
    int getCombinationFromULong(ULONG * xs, ULONG combUlong) const;
private:

};

#endif	/* COMBINATORIALINDEXER_H */

