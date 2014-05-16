/* 
 * File:   CombinatiorialGenerator.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 2:05 PM
 */

#ifndef COMBINATIORIALGENERATOR_H
#define	COMBINATIORIALGENERATOR_H

#include "base.h"
#include <stdio.h>
#include <limits.h>

class CombinatiorialGenerator {
private:
    /**
     * Number of possible elements to choose.
     */
    ULONG up;
    
    /**
     * Number of elements we want to choose.
     */
    ULONG down;
    
    /**
     * Total number of all combinations.
     */
    ULONG totalNum;
    
    /**
     * Number of bytes to represent the combination. 
     */
    ULONG byteWidth;  
    
    /**
     * if false then the generator is before beginning, next has to be called to
     * get to the first combination.
     */
    bool started;
    
    /**
     * Ordering number of the current combination.
     */
    ULONG counter;
    
    /**
     * Current state = indexes of the individual choices (e.g., {1,2,3} for initial combination X over 3).
     * Length of this array = down.
     */
    ULONG * curState;
    
    /**
     * Current combination generated - bit array of given size with given combination.
     * Length of this array = CEIL(up/8).
     */
    uchar * curCombination;
    bool curCombinationValid;
    
    /**
     * Internal state of the algorithm for generating all combinations.
     */
    ULONG * c;
    //long j;
    //long x;
    
    void firstCombination();
    bool internalNext();
    
public:
    CombinatiorialGenerator(ULONG up, ULONG down);
    virtual ~CombinatiorialGenerator();
    
    // Reset the internal state.
    void reset();
    
    // Getters.
    inline ULONG getUp()                { return this->up;              }
    inline ULONG getDown()              { return this->down;            }
    inline ULONG getTotalNum()          { return this->totalNum;        }
    inline ULONG getByteWidth()         { return this->byteWidth;       }
    inline ULONG getCounter()           { return this->counter;         }
    inline ULONG * getCurState()        { return this->curState;        }
    inline ULONG * getC()               { return this->c;  }
    uchar * getCurCombination();
    
    // Move to the next.
    bool next();
    
    /**
     * Compute binomial coefficient.
     * 
     * @param n
     * @param k
     * @return 
     */
    static ULONG binomial(ULONG n, ULONG k);
    
};

#endif	/* COMBINATIORIALGENERATOR_H */

