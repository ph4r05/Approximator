/* 
 * File:   CombinatorialIndexer.cpp
 * Author: ph4r05
 * 
 * Created on May 28, 2014, 11:08 AM
 */

#include "CombinatorialIndexer.h"
#include "CombinatiorialGenerator.h"
#include <string>
#include <cassert>
#include <cmath>

using namespace std;

CombinatorialIndexer::CombinatorialIndexer() : down(0), up(0), logBitInputWidth(0), binomialSums(NULL){
}

CombinatorialIndexer::~CombinatorialIndexer() {
    deinit();
}

void CombinatorialIndexer::init(uint up, uint down) {
    // Deinitialize first (may have been initialized before)
    deinit();
    
    this->down = down;
    this->up = up;
    this->logBitInputWidth = ceil(log2(up));
    
    // Pre-compute binomial sums for combination index computation.
    if (down >= 1){
        binomialSums = new ULONG * [down];
        for(uint order = 1; order <= down; order++){
            const uint idx = order-1;
            binomialSums[idx]    = new ULONG[up];
            binomialSums[idx][0] = 0ul;

            ULONG res = 0;
            for(uint i = 1; i<up; i++){
                res += CombinatiorialGenerator::binomial(up-i, order-1);
                binomialSums[idx][i] = res;
            }
        }
    }
}

void CombinatorialIndexer::deinit() {
    if (binomialSums!=NULL){
        for(uint order = 1; order <= down; order++){
            delete[] binomialSums[order-1];
            binomialSums[order-1]=NULL;
        }
        
        delete[] binomialSums;
        binomialSums = NULL;
    }
}


inline ULONG CombinatorialIndexer::getCubeIdx(ULONG x1, ULONG x2, ULONG x3) const {
    return down < 3 ? 0 : binomialSums[2][x1] + CombinatiorialGenerator::getQuadIdx(up-1-x1, x2-x1-1, x3-x1-1);
}

ULONG CombinatorialIndexer::getCombinationIdx(uint order, const ULONG* xs, uint xsOffset, ULONG Noffset, ULONG combOffset) const {
    assert(order<=down);
    
    // Small order.
    if (order==0) return 0;
    if (order==1) return xs[0+xsOffset] - combOffset;
    if (order==2) return CombinatiorialGenerator::getQuadIdx(up-Noffset, xs[0+xsOffset]-combOffset, xs[1+xsOffset]-combOffset);
    
    // Order 3.
    // Take recursive parameters into consideration.
    const ULONG x1 = xs[0+xsOffset] - combOffset;
    
    // This is simple optimization to remove 1 recursion step.
    if (order==3) {
        const ULONG sum = binomialSums[2][x1+Noffset] - binomialSums[2][Noffset];
        return sum + CombinatiorialGenerator::getQuadIdx(
                up-1-x1-Noffset, 
                xs[1+xsOffset]-combOffset-x1-1, 
                xs[2+xsOffset]-combOffset-x1-1);
    }
    
    // Order 3 and higher.
    // Sum of the previous combinations is shifted due to Noffset.
    // N is reduced thus the binomial sum is shifted (originating point is smaller).
    return (binomialSums[order-1][x1+Noffset] - binomialSums[order-1][Noffset]) + 
            getCombinationIdx(
            order-1, 
            xs, 
            xsOffset+1, 
            Noffset+1+x1, 
            combOffset+1+x1
           );
}

int CombinatorialIndexer::getCombinationFromIdx(uint order, ULONG* xs, ULONG idx) const {
    int nOffset=0;
    for(int x=order-1; x>=1; x--){
        int i=up-1-nOffset;
        for(; i>=0; i--){
            const ULONG biSum = (binomialSums[x][i+nOffset] - binomialSums[x][nOffset]);
            // If current index is bigger than the current sum, the previous sum
            // helps us to determine index of the order.
            if (idx >= biSum){
                break;
            }
        }

        idx-=(binomialSums[x][i+nOffset] - binomialSums[x][nOffset]);
        
        xs[order-x-1] = i+nOffset;
        nOffset = xs[order-x-1]+1;
    }
    
    // Final element is already determined.
    xs[order-1] = idx + nOffset;
    return 1;
}

ULONG CombinatorialIndexer::getCombinationULong(uint order, const ULONG* xs) const {
    assert(order <= 0xf);
    ULONG res = order & 0xf;
    
    uint offset = 4;  // order.
    for(uint x=0; x<order; x++){
        res |= xs[x] << offset;
        offset += logBitInputWidth;
    }
    
    return res;
}

int CombinatorialIndexer::getCombinationFromULong(ULONG* xs, ULONG combUlong) const {
    uint order = combUlong & 0x7;
    
    combUlong = combUlong >> 4; // remove order.
    for(uint x=0; x<order; x++){
        xs[x] = combUlong & ((ULONG1 << logBitInputWidth)-1);
        combUlong = combUlong >> logBitInputWidth;
    }
    
    return 1;
}


