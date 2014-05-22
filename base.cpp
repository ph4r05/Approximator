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

