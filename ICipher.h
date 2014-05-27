/* 
 * File:   ICipher.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 11:37 AM
 */

#ifndef ICIPHER_H
#define	ICIPHER_H

#include "base.h"

#define GETBIT(bitPos, inp) ((((inp)+(bitPos/8)) & (ULONG1 << bitPis%8)) > 0)

class ICipher {
public:
    ICipher();
    ICipher(const ICipher& orig);
    virtual ~ICipher();
    
    virtual unsigned getInputBlockSize() = 0;
    virtual unsigned getOutputBlockSize() = 0;
    virtual unsigned getKeyBlockSize() = 0;
    virtual int setNumRounds(int rounds) = 0;
    virtual int evaluate(const unsigned char * input, unsigned char * output) = 0;
    virtual int evaluate(const unsigned char * input, const unsigned char * key, unsigned char * output) = 0;
private:

};

#endif	/* ICIPHER_H */

