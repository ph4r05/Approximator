/* 
 * File:   KeccakFull.h
 * Author: ph4r05
 *
 * Created on June 19, 2014, 2:23 PM
 */

#ifndef KECCAKFULL_H
#define	KECCAKFULL_H


#include <inttypes.h>
#include <string.h>
#include "base.h"
#include "ICipher.h"

#ifndef uint8_t
typedef unsigned char uint8_t;
#endif

#ifndef ROTL64
#define ROTL64(x, y) (((x) << (y)) | ((x) >> (64 - (y))))
#endif

class KeccakFull : public ICipher {
private:
    int rounds, r, c;
    
public:
    KeccakFull();
    virtual ~KeccakFull();
        
    virtual unsigned getId() const { return 2; }
    virtual unsigned getInputBlockSize() const        { return 200; };
    virtual unsigned getOutputBlockSize() const       { return 200; };
    virtual unsigned getKeyBlockSize()  const         { return 0; };
    virtual int setNumRounds(int rounds)              { this->rounds = rounds; return 1; };
    virtual int getNumRounds() const                  { return this->rounds; }
        
    void set_parameters(int r, int c);
    
    // update the state with given number of rounds
    void keccakf(uint64_t st[25], int rounds) const;

    // compute a keccak hash (md) of given byte length from "in"
    int keccak(const uint8_t *in, int inlen, uint8_t *md, int mdlen) const;
    
    virtual int evaluate(const unsigned char * input, unsigned char * output) const;
    virtual int evaluate(const unsigned char *input, const unsigned char *key, unsigned char *output ) const;
    
    inline virtual int prepareKey(const unsigned char * key) 
    { return 1; }
    
    inline virtual int evaluateWithPreparedKey(const unsigned char * input, unsigned char * output) const
    { return evaluate(input, NULL, output); return 1;}
};

#endif	/* KECCAKFULL_H */

