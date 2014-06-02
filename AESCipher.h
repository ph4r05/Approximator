/* 
 * File:   AESCipher.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 11:42 AM
 */

#ifndef AESCIPHER_H
#define	AESCIPHER_H

#include "base.h"
#include "ICipher.h"
#include "aes.h"
#include <openssl/aes.h>

class AESCipher : public ICipher {
private:
        int rounds;
        unsigned int key_schedule[60];
        
public:
    AESCipher();
    AESCipher(const AESCipher& orig);
    virtual ~AESCipher();
    
    virtual unsigned getInputBlockSize()        { return AES_BLOCK_SIZE; };
    virtual unsigned getOutputBlockSize()       { return AES_BLOCK_SIZE; };
    virtual unsigned getKeyBlockSize()          { return AES_BLOCK_SIZE; };
    virtual int setNumRounds(int rounds)        { this->rounds = rounds; return 1; };
    
    virtual int evaluate(const unsigned char * input, unsigned char * output);
    virtual int evaluate(const unsigned char * input, const unsigned char * key, unsigned char * output);
    
    inline virtual int prepareKey(const unsigned char * key) 
    { KeyExpansion(key,key_schedule, 128); return 1; }
    
    inline virtual int evaluateWithPreparedKey(const unsigned char * input, unsigned char * output) 
    { aes_encrypt(input, output, this->key_schedule, 128, this->rounds); return 1;}

};

#endif	/* AESCIPHER_H */

