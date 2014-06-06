#include "base.h"
#include <iostream>
#include <iomanip>
using namespace std;

void dumpUchar(std::ostream & c, const uchar* inp, unsigned int size, bool endl) {
    for(unsigned int i = 0; i < size; i++){
        c << (unsigned int) inp[i] << " ";
    }
    
    if (endl){
        c << std::endl;
    }
}

void readUcharToUlong(const uchar * input, uint size, ULONG * iBuff) {
    // At first reset memory of the big buffer.
    memset(iBuff, 0, size);
    // And read particular parts of the input to the big buffer.
    for(uint x = 0; x < size; x++){
        iBuff[x/SIZEOF_ULONG] = READ_TERM_1(iBuff[x/SIZEOF_ULONG], input[x], x%SIZEOF_ULONG);
    }
}

void readUlongToUchar(uchar* output, uint size, const ULONG* iBuff) {
    for(uint x=0; x<size; x++){
        output[x] = (iBuff[x/SIZEOF_ULONG] >> (8* (x % SIZEOF_ULONG))) & ((unsigned char)0xffu);
    }
}

void randomBuffer(uchar * buffer, uint size){
    for(unsigned int k=0; k<size; k++){ 
        buffer[k] = (rand() % (0x100ul)); 
    }
}
