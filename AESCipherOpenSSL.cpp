/* 
 * File:   AESCipherOpenSSL.cpp
 * Author: ph4r05
 * 
 * Created on May 27, 2014, 10:53 AM
 */

#include "AESCipherOpenSSL.h"
#include <string.h>

AESCipherOpenSSL::AESCipherOpenSSL() {
}

AESCipherOpenSSL::AESCipherOpenSSL(const AESCipherOpenSSL& orig) {
}

AESCipherOpenSSL::~AESCipherOpenSSL() {
}

int AESCipherOpenSSL::evaluate(const unsigned char* input, unsigned char* output) {
    unsigned char inp[AES_BLOCK_SIZE];
    unsigned char key[AES_BLOCK_SIZE];
    memcpy(inp, input,                  AES_BLOCK_SIZE);
    memcpy(key, input + AES_BLOCK_SIZE, AES_BLOCK_SIZE);
    
    return evaluate(inp, key, output);
}

int AESCipherOpenSSL::evaluate(const unsigned char* input, const unsigned char* key, unsigned char* output) {
    AES_KEY enc_key;
    AES_set_encrypt_key(key, getKeyBlockSize()*8, &enc_key);
    AES_encrypt(input, output, &enc_key);
    return 1;
}

