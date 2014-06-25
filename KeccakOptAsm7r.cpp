/* 
 * File:   KeccakOptAsm7r.cpp
 * Author: ph4r05
 * 
 * Created on June 23, 2014, 2:40 PM
 */

#include "KeccakOptAsm7r.h"
extern "C" {
#include "sha3/hash_functions/Keccak64_common/KeccakF-1600-interface.h"
#include "sha3/hash_functions/Keccak64_common/KeccakSponge.h"
}

KeccakOptAsm7r::KeccakOptAsm7r() : rounds(7) {
}

KeccakOptAsm7r::~KeccakOptAsm7r() {
}

int KeccakOptAsm7r::setNumRounds(int rounds) { 
    if (rounds!=1 && rounds!=2 && rounds!=3 && rounds!=7 && rounds!=8 && rounds!=24){
        return -1;
    }
    
    this->rounds=rounds; 
    return 1; 
}

int KeccakOptAsm7r::evaluate(const unsigned char* input, unsigned char* output) const {
    unsigned char ou[128] = {0};
    spongeState state;
    InitSponge(&state, 1024, 576);
    Absorb(&state, input, 32*8, 7);
    Squeeze(&state, ou, 1024, 7);
    memcpy(output, ou, 16);
    return 1;
}

int KeccakOptAsm7r::evaluate(const unsigned char* input, const unsigned char* key, unsigned char* output) const {
    unsigned char in[32] = {0};
    unsigned char ou[128] = {0};
    memcpy(in, input, 16);
    memcpy(in, key, 16);
    
    spongeState state;
    InitSponge(&state, 1024, 576);
    Absorb(&state, in, 32*8, 7);
    Squeeze(&state, ou, 1024, 7);
    memcpy(output, ou, 16);
    return 1;
}
