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
    this->totalNum = binomial(up, down);
    
    this->curCombinationValid = false;
    this->counter = 0;
    this->curState = new ULONG[down];
    this->curCombination = new uchar[byteWidth];
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
}

ULONG CombinatiorialGenerator::binomial(ULONG n, ULONG k) {
    ULONG r = 1, d = n - k;

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
    memset(curState, 0, sizeof(ULONG) * down);
    memset(curCombination, 0, byteWidth);
        
    // Reset counter
    counter = 0;
    started=false;
    
    // Knuth lexicographic ordering initialization.
//    for(unsigned j = 1; j<=down; j++){
//        c[j] = j-1;
//    }
//    
//    c[down+1] = up;
//    c[down+2] = 0;
//    this->j = down-1;
//    this->x = 0;
}

const uchar * CombinatiorialGenerator::getCurCombination() {
    // Generating combinations only on demand.
    // Using cached version of the bit representation of the current combination.
    if (curCombinationValid) return curCombination;
    
    // Set bit representation from the current state.
    memset(curCombination, 0, byteWidth);
    for(unsigned i = 0; i<down; i++){
        //cout << "; c_"<<i<<"=" << curState[i] << " ";
        curCombination[ curState[i]/8 ] |= 1u << (curState[i]%8);
    }
    curCombinationValid=true;
    return curCombination;
}

void CombinatiorialGenerator::firstCombination() {
    // Initialize current state properly.
    for(unsigned j = 0; j<down; j++){
        curState[j] = j;   
    }
}

bool CombinatiorialGenerator::next() {
    // Very special case, down=0, at the end right now.
    if (down<=0){
        return false;
    }
    
    // If the generator was not started, calling
    // the first next() moves it to the first combination.
    if (!started){
        firstCombination();
        started=true;
        counter+=1;
        return true;
    }
    
    // Move internal state to the next combination.
    bool inext = internalNext();
    
    // If is already terminated, nothing to do next.
    if (!inext) return inext;
    
    // Move happened.
    // Invalidate current bit representation of the combination.
    curCombinationValid=false;
    counter+=1;    
    return true;
}


bool CombinatiorialGenerator::internalNext() {    
    // If we are on the end, 
    curState[down-1]+=1;
    
    // If the last counter overflowed, switch is needed.
    if (curState[down-1] > up-1){
        
        long x=down-1;
        while(curState[x] >= up-down+x && x>=0){
            x-=1;        
        }

        // Terminating algorithm? All combinations were generated.
        if (x<0){
            return false;
        }

        // If some shift happened, do the shifting.
        if (x!=down-1){
            curState[x]+=1;
            for(unsigned y=1; y<=(down-1-x); y++){
                curState[x+y]=curState[x]+y;
            }
        }
    }
    
    return true;
    
//    //
//    // Knuth combinatorial generator was here but it is not working!
//    //
//    
//    // T2 of the algorithm.
//    if (j>0){
//        x = j;
//        c[j] = x; 
//        j -= 1;
//        cout << "t2; j="<<j<<"; x="<< x << "; c[j+1]=" << c[j+1] <<  endl;
//        return true;
//    }
//    
//    // T3[Easy case?] 
//    if ((c[1]+1) < c[2]){
//        c[1] += 1;
//        cout << "easy" << endl;
//        return true;
//    } else {
//        cout << "nnn" << endl;
//        j = 2;
//    }
//    
//    // T4[Find j].
//    do {
//        c[j-1] = j-2;
//        x      = c[j]+1;
//        if (x == c[j+1]){
//            cout << "Equals x=" << x << endl;
//            j += 1;
//        }
//        cout << "cc" << endl;
//    } while(x == c[j+1]);
//    cout << "current j: " << j << endl;
//    
//    // T5[Done?]
//    if (j>down){
//        cout << "Terminating " << endl;
//        return false;
//    }
//    
//    // T6. [Increase cj].
//    c[j] = x; 
//    cout << " increasing c[" <<j<<"]=" << x << endl;
//    j -= 1;
//    return true;
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
    
    res += getQuadIdx(N-1-x1, x2, x3);
    return res;
}