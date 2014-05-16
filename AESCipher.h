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
#include <openssl/aes.h>

class AESCipher : public ICipher {
public:
    AESCipher();
    AESCipher(const AESCipher& orig);
    virtual ~AESCipher();
    
    virtual unsigned getInputBlockSize()        { return AES_BLOCK_SIZE; };
    virtual unsigned getOutputBlockSize()       { return AES_BLOCK_SIZE; };
    virtual unsigned getKeyBlockSize()          { return AES_BLOCK_SIZE; };
    
    virtual int evaluate(const unsigned char * input, unsigned char * output);
    virtual int evaluate(const unsigned char * input, const unsigned char * key, unsigned char * output);
private:

};

#endif	/* AESCIPHER_H */

