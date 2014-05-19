/* 
 * File:   main.cpp
 * Author: ph4r05
 *
 * Created on May 16, 2014, 10:20 AM
 */

#include <cstdlib>

#include "Approximation.h"
#include "AESCipher.h"
#include "CombinatiorialGenerator.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    AESCipher c;
    Approximation ap;
    ap.setCipher(&c);
    
    ap.work();
    ap.genMessages();
    
    /*cout << "Start" << endl;
    for(uint i=0; i<1000; i++){
        CombinatiorialGenerator cgen(256, 2);
        for(uint k=0; cgen.next(); k++){
            ;
        }
    }
    cout << "end" << endl;*/
    
    return 0;
}

