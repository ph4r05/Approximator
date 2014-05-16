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
    
    return 0;
}

