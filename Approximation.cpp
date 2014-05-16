/* 
 * File:   Approximation.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 10:21 AM
 */

#include "Approximation.h"
#include "CombinatiorialGenerator.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include <cstring>
#include <cstdlib> 
#include <ctime> 

using namespace std;

Approximation::Approximation() : cip(NULL) {
}

Approximation::~Approximation() {
}

void Approximation::work() {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    finput         = new uchar[byteWidth];
    uchar * key    = finput + cip->getInputBlockSize();
    
    // Generate ciphertext for constant    
    memset(finput, 0, byteWidth);
    cip->evaluate(finput, key, output);
    
    // Dump
    ofstream hist("cip_0.txt");
    dumpUchar(hist, output, cip->getOutputBlockSize());
    hist.close();
    
    // Generate order1 .. order3 cipher data
    for(unsigned order=1; order<=3; order++){
        CombinatiorialGenerator cg(byteWidth*8, order);
        
        cout << "Starting with order: " << order << endl;
        ofstream cip1(std::string("cip_") + std::to_string(order) + ".txt");
        for(ULONG ctr=0; cg.next(); ctr++){
            uchar * input = cg.getCurCombination();
            
            // Evaluate cipher on current combinations.
            cip->evaluate(input, input + cip->getInputBlockSize(), output);
            
            // Dump
            dumpUchar(cip1, output, cip->getOutputBlockSize());
        }
        
        cip1.close();
    }
    
    cout << "Generating finished" << endl;
}

void Approximation::genMessages() {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    finput         = new uchar[byteWidth];
    uchar * key    = new uchar[cip->getKeyBlockSize()];
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    ofstream cip3("cip_msgs.txt");
    
    // Generate key at random.
    for(unsigned int i=0; i<cip->getKeyBlockSize(); i++){ 
        key[i] = (rand() % (0xffu+1)); 
    }
    
    dumpUchar(cip3, key, cip->getKeyBlockSize());
    
    // Generate tons of random messages.
    for(unsigned long i=0; i<16777216ul; i++){
        // Generate message at random.
        for(unsigned int i=0; i<cip->getInputBlockSize(); i++){ 
            finput[i] = (rand() % (0xffu+1)); 
        }
        
        cip->evaluate(finput, key, output);
        
        // Dump
        dumpUchar(cip3, output, cip->getOutputBlockSize());
    }
    
    delete[] key;
    delete[] output;
}

