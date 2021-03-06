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
    
    virtual unsigned getId() const = 0;
    virtual unsigned getInputBlockSize() const = 0;
    virtual unsigned getOutputBlockSize() const = 0;
    virtual unsigned getKeyBlockSize() const = 0;
    virtual int setNumRounds(int rounds) = 0;
    virtual int getNumRounds() const = 0;
    
    virtual int evaluate(const unsigned char * input, unsigned char * output) const = 0;
    virtual int evaluate(const unsigned char * input, const unsigned char * key, unsigned char * output) const = 0;
    
    virtual int prepareKey(const unsigned char * key) = 0;
    virtual int evaluateWithPreparedKey(const unsigned char * input, unsigned char * output) const = 0;
private:

};

#endif	/* ICIPHER_H */

