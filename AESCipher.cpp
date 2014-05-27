/* 
 * File:   AESCipher.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 11:42 AM
 */

#include <string.h>

#include "aes.h"
#include "AESCipher.h"

AESCipher::AESCipher() : rounds(-1) {
    
}

AESCipher::AESCipher(const AESCipher& orig) : rounds(orig.rounds) {
    
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
   unsigned int key_schedule[60];
   
   // Generate enc key.
   KeyExpansion(key,key_schedule, 128);
   
   // Encrypt.
   aes_encrypt(input, output, key_schedule, 128, this->rounds);
   
   return 1;
}

