/* 
 * File:   KeccakTest.cpp
 * Author: ph4r05
 *
 * Created on June 25, 2014, 1:50 PM
 */

#include <cstdlib>
#include "base.h"

extern "C" {
#include "sha3/hash_functions/Keccak64_common/KeccakF-1600-interface.h"
#include "sha3/hash_functions/Keccak64_common/KeccakSponge.h"
}

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    unsigned int i=0, j=0, rounds, seed;
    
    if (argc == 3) {
        seed = atoi(argv[1]);
        rounds = atoi(argv[2]);
    } else {
        fprintf(stderr, "Error, 2 parameters needed: $seed $rounds\n");
        return 1;
    }
    
    fprintf(stderr, "Seed=%u, rounds=%u\n", seed, rounds);
    srandom(seed);
    for(i=0; i<128; i++){
        unsigned char ou[128] = {0};
        unsigned char input[200] = {0};
        unsigned char output[16] = {0};
        
        // Random init of the input
        for(j=0; j<32; j++){
            input[j] = rand() % 0x100;
        }
        
        spongeState state;
        InitSponge(&state, 1024, 576);
        Absorb(&state, input, 32*8, rounds);
        Squeeze(&state, ou, 1024, rounds);
        memcpy(output, ou, 16);
        
        // Dump
        dumpHex(cout, input, 32);
        dumpHex(cout, output, 16);
        cout << endl;
    }
    
    return 0;
}

