/* 
 * File:   KeccakOptAsm.cpp
 * Author: ph4r05
 * 
 * Created on June 23, 2014, 2:40 PM
 */

#include "KeccakOptAsm.h"
extern "C" {
#include "sha3/hash_functions/Keccak64_common/KeccakF-1600-interface.h"
#include "sha3/hash_functions/Keccak64_common/KeccakSponge.h"
}

KeccakOptAsm::KeccakOptAsm() : rounds(7) {
}

KeccakOptAsm::~KeccakOptAsm() {
}

int KeccakOptAsm::setNumRounds(int rounds) { 
    if (rounds!=24 && (rounds<1 || rounds>8)){
        return -1;
    }
    
    this->rounds=rounds; 
    return 1; 
}

int KeccakOptAsm::evaluate(const unsigned char* input, unsigned char* output) const {
    unsigned char ou[128] = {0};
    spongeState state;
    InitSponge(&state, 1024, 576);
    Absorb(&state, input, 32*8, rounds);
    Squeeze(&state, ou, 1024, rounds);
    memcpy(output, ou, 16);
    return 1;
}

int KeccakOptAsm::evaluate(const unsigned char* input, const unsigned char* key, unsigned char* output) const {
    unsigned char in[32] = {0};
    unsigned char ou[128] = {0};
    memcpy(in,    input, 16);
    memcpy(in+16, key, 16);
    
    spongeState state;
    InitSponge(&state, 1024, 576);
    Absorb(&state, in, 32*8, rounds);
    Squeeze(&state, ou, 1024, rounds);
    memcpy(output, ou, 16);
    return 1;
}
