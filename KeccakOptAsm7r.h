/* 
 * File:   KeccakOptAsm7r.h
 * Author: ph4r05
 *
 * Created on June 23, 2014, 2:40 PM
 */

#ifndef KECCAKOPTASM7R_H
#define	KECCAKOPTASM7R_H

#include <inttypes.h>
#include <string.h>
#include "base.h"
#include "ICipher.h"

class KeccakOptAsm7r : public ICipher {
private:
    unsigned rounds;
    
public:
    KeccakOptAsm7r();
    virtual ~KeccakOptAsm7r();
        
    virtual unsigned getId() const { return 3; }
    virtual unsigned getInputBlockSize() const        { return 16; };
    virtual unsigned getOutputBlockSize() const       { return 16; };
    virtual unsigned getKeyBlockSize()  const         { return 16; };
    virtual int getNumRounds() const                  { return 7; }
    virtual int setNumRounds(int rounds);
    
    virtual int evaluate(const unsigned char * input, unsigned char * output) const;
    virtual int evaluate(const unsigned char *input, const unsigned char *key, unsigned char *output ) const;
    
    inline virtual int prepareKey(const unsigned char * key) 
    { return 1; }
    
    inline virtual int evaluateWithPreparedKey(const unsigned char * input, unsigned char * output) const
    { return evaluate(input, NULL, output); return 1;}
};

#endif	/* KECCAKOPTASM7R_H */
