/* 
 * File:   AESCipherOpenSSL.h
 * Author: ph4r05
 *
 * Created on May 27, 2014, 10:53 AM
 */

#ifndef AESCIPHEROPENSSL_H
#define	AESCIPHEROPENSSL_H

#include "base.h"
#include "ICipher.h"
#include <openssl/aes.h>

class AESCipherOpenSSL {
public:
    AESCipherOpenSSL();
    AESCipherOpenSSL(const AESCipherOpenSSL& orig);
    virtual ~AESCipherOpenSSL();
    
    virtual unsigned getId() const { return 0; }
    virtual unsigned getInputBlockSize()        { return AES_BLOCK_SIZE; };
    virtual unsigned getOutputBlockSize()       { return AES_BLOCK_SIZE; };
    virtual unsigned getKeyBlockSize()          { return AES_BLOCK_SIZE; };
    virtual int setNumRounds(int rounds)        { return -1;             };
    virtual int getNumRounds() const            { return -1; }
    
    virtual int evaluate(const unsigned char * input, unsigned char * output);
    virtual int evaluate(const unsigned char * input, const unsigned char * key, unsigned char * output);
private:

};

#endif	/* AESCIPHEROPENSSL_H */

