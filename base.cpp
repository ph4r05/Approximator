#include "base.h"
#include <iostream>
#include <iomanip>
using namespace std;

void dumpUcharHex(std::ostream & c, uchar* inp, unsigned int size) {
    c << showbase // show the 0x prefix
         << internal // fill between the prefix and the number
         << setfill('0'); // fill with 0s
    
    for(unsigned int i = 0; i < size; i++){
        c << hex << setw(4) << (unsigned int) inp[i] << " ";
    }
    
    c << endl;
}

void dumpUlongHex(std::ostream & c, ULONG * inp, unsigned int size) {
    c << showbase // show the 0x prefix
         << internal // fill between the prefix and the number
         << setfill('0'); // fill with 0s
    
    for(unsigned int i = 0; i < size; i++){
        c << hex << setw(8) << (ULONG) inp[i] << " ";
    }
    
    c << endl;
}

void dumpUchar(std::ostream & c, uchar* inp, unsigned int size) {
    for(unsigned int i = 0; i < size; i++){
        c << (unsigned int) inp[i] << " ";
    }
    
    c << endl;
}

