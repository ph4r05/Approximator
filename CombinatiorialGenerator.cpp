/* 
 * File:   CombinatiorialGenerator.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 2:05 PM
 */

#include <string.h>
#include <cassert>
#include <iostream>
#include "CombinatiorialGenerator.h"
using namespace std;

CombinatiorialGenerator::CombinatiorialGenerator(ULONG up, ULONG down) {
    assert(up >= down);
    
    this->up = up;
    this->down = down;
    this->byteWidth = OWN_CEIL((double)up / 8.0);
    this->byteUlongWidth = OWN_CEIL((double)up / (8.0 * SIZEOF_ULONG));
    this->totalNum = binomial(up, down);
    
    this->curCombinationValid = false;
    this->curUlongCombinationValid = false;
    this->counter = 0;
    
    this->curState = down>0 ? new ULONG[this->down] : NULL;
    this->curCombination = new uchar[this->byteWidth];
    this->curUlongCombination = new ULONG[this->byteUlongWidth];
    this->reset();
}

CombinatiorialGenerator::~CombinatiorialGenerator() {
    if (this->curState!=NULL){
        delete[] this->curState;
        this->curState = NULL;
    }
    
    if (this->curCombination != NULL){
        delete[] this->curCombination;
        this->curCombination = NULL;
    }
    
    if (this->curUlongCombination != NULL){
        delete[] this->curUlongCombination;
        this->curUlongCombination = NULL;
    }
}

ULONG CombinatiorialGenerator::binomial(ULONG n, ULONG k) {
    ULONG r = 1, d = n - k;
    if (k==0) return ULONG1;
    if (k==1) return n;
    if (k==2) return (n*(n-1))/2ul;
    if (k==3) return (n*(n-1)*(n-2))/6ul;    
    
    /* choose the smaller of k and n - k */
    if (d > k) { k = d; d = n - k; }

    while (n > k) {
        if (r >= UINT_MAX / n) return 0; /* overflown */
        r *= n--;

        /* divide (n - k)! as soon as we can to delay overflows */
        while (d > 1 && !(r % d)) r /= d--;
    }
    return r;
}

void CombinatiorialGenerator::reset() {
    // Reset internal state, set to the first 
    if (down>0){
        memset(curState, 0, sizeof(ULONG)*down);
    }
    memset(curCombination, 0, sizeof(uchar)*byteWidth);
    memset(curUlongCombination, 0, sizeof(ULONG)*byteUlongWidth);
        
    // Reset counter
    counter = 0;
    started=false;
    curCombinationValid=false;
    curUlongCombinationValid=false;
}

const uchar * CombinatiorialGenerator::getCurCombination() {
    // Generating combinations only on demand.
    // Using cached version of the bit representation of the current combination.
    if (curCombinationValid) return curCombination;
    
    // Set bit representation from the current state.
    memset(curCombination, 0, sizeof(uchar)*byteWidth);
    for(unsigned i = 0; i<down; i++){
        //cout << "; c_"<<i<<"=" << curState[i] << " ";
        curCombination[ curState[i]/8 ] |= ULONG1 << (curState[i]%8);
    }
    curCombinationValid=true;
    return curCombination;
}

const ULONG* CombinatiorialGenerator::getCurUlongCombination() {
    // Generating combinations only on demand.
    // Using cached version of the bit representation of the current combination.
    if (curUlongCombinationValid) return curUlongCombination;
    
    // Set bit representation from the current state.
    memset(curUlongCombination, 0, SIZEOF_ULONG*byteUlongWidth);
    for(unsigned i = 0; i<down; i++){
        curUlongCombination[ curState[i]/(8*SIZEOF_ULONG) ] |= ULONG1 << (curState[i]%(8*SIZEOF_ULONG));
    }
    curUlongCombinationValid=true;
    return curUlongCombination;
}

void CombinatiorialGenerator::firstCombination() {
    // Initialize current state properly.
    for(unsigned j = 0; j<down; j++){
        curState[j] = j;   
    }
}

bool CombinatiorialGenerator::next() {
    // Very special case, down=0, at the end right now.
    if (down<0){
        return false;
    }
    
    // If the generator was not started, calling
    // the first next() moves it to the first combination.
    if (!started){
        firstCombination();
        started=true;
        counter=0;
        return true;
    } else if(down==0) {
        return false;
    }
    
    //
    // Move internal state to the next combination.
    //
    bool inext = true;
    
    // Base case: move the last element of the combination.
    // E.g., [1,2,3,4] -> [1,2,3,5].
    curState[down-1]+=1;
    
    // If the last counter overflowed, switch is needed.
    // E.g., [1,2,3,128] for up=128. 128 is not legal, thus overflowed.
    if (curState[down-1] > up-1){
        
        // Find nearest "digit" that does not overflow and can be incremented.
        long x=down-1;
        while(curState[x] >= up-down+x && x>=0){
            x-=1;        
        }

        // Terminating algorithm? All combinations were generated.
        if (x<0){
            inext=false;
        }

        // If some shift happened, do the shifting.
        if (inext && x!=down-1){
            // Increment the non-overflowing digit.
            curState[x]+=1;
            // Since we have combinations here, all digits
            // to the right from the newly shifted digit have
            // to have the least possible combination.
            // E.g., [1,2,126,127] -> [1,3,4,5]
            for(unsigned y=1; y<=(down-1-x); y++){
                curState[x+y]=curState[x]+y;
            }
        }
    }
    
    // If is already terminated, nothing to do next.
    if (!inext) return inext;
    
    // Move happened.
    // Invalidate current bit representation of the combination.
    curCombinationValid=false;
    curUlongCombinationValid = false;
    counter+=1;    
    return true;
}

ULONG CombinatiorialGenerator::getQuadIdx(ULONG N, ULONG x1, ULONG x2) {
    // Index of this combination is sum of all previous combinations (N-2 + N-3 + ...)
    // plus the ordering number of the x1x2 combination from the beginning of the x1 
    // starting combinations.
    //
    // The result is the same as SUM_{i=0}^{x1-1} Binomial(N-1-i, 1) + (x2-x1);
    // 
    const ULONG n = N-1;
    ULONG idx = ((2*n*x1-x1*x1+x1)/2) + (x2-x1) - 1;
    return idx;
}

ULONG CombinatiorialGenerator::getCubeIdx(ULONG N, ULONG x1, ULONG x2, ULONG x3) {
    // Index of this combination is sum of all previous combinations
    // plus the ordering number of the x1x2x3 combination from the beginning of the x1
    // starting combinations.
    //
    // All previous combinations not starting with x1: SUM_{i=1}^{x1} Binomial(N-i, 2)
    // because we are looking for combinations of positions for x2x3, while
    // the space is decreasing by 1 since x1 is moving also.
    // 
    // Example, assume N=128.
    // if x1=1, then we want to compute all previous combinations, i.e., when x1=0 is fixed
    // and x2,x3 are floating. x1 already takes the first bit, thus there are 127 
    // remaining bits. There are Binomial(127, 2) combinations for x2x3, thus
    // if x1=1 there were Binomial(127,2) combination before x1=1,x2=2,x3=3
    // (first combination with x1=1 in the lexicographic ordering).
    //
    // result = SUM_{i=0}^{x1-1} Binomial(N-1-i, 2) + getQuadIdx(N-x1-1, x2, x3)
    
    ULONG res = 0;
    for(uint i = 1; i<=x1; i++){
        res += binomial(N-i, 2);
    }
    
    res += getQuadIdx(N-1-x1, x2-x1-1, x3-x1-1);
    return res;
}
