/* 
 * File:   AESCipher.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 11:42 AM
 */

#include <string.h>

#include "AESCipher.h"

AESCipher::AESCipher() {
    
}

AESCipher::AESCipher(const AESCipher& orig) {
    
}

AESCipher::~AESCipher() {
    
}

int AESCipher::evaluate(const unsigned char* input, unsigned char* output) {
    unsigned char inp[AES_BLOCK_SIZE];
    unsigned char key[AES_BLOCK_SIZE];
    memcpy(inp, input,                  AES_BLOCK_SIZE);
    memcpy(key, input + AES_BLOCK_SIZE, AES_BLOCK_SIZE);
    
    return evaluate(inp, key, output);
}

int AESCipher::evaluate(const unsigned char* input, const unsigned char* key, unsigned char* output) {
    AES_KEY enc_key;
    AES_set_encrypt_key(key, getKeyBlockSize()*8, &enc_key);
    AES_encrypt(input, output, &enc_key);
    return 1;
}

