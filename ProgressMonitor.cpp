/* 
 * File:   ProgressMonitor.cpp
 * Author: ph4r05
 * 
 * Created on May 19, 2014, 1:52 PM
 */

#include "ProgressMonitor.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

ProgressMonitor::~ProgressMonitor() {
}

void ProgressMonitor::setCur(double current) {
    assert(step>0);
    
    if (current >= (last+step)){
        // How many dost to echo ?
        const unsigned toEcho = floor((current+0.00001-last) / step);
        //cout << "curr="<<current<<"; last="<<last<<"; toecho="<<toEcho<<"; step="<<step<<endl;
        
        for(unsigned i=0; i<toEcho; i++){
            cout << "." << flush;
        } 
        cout << flush;
        
        // Set last to the low multiple of step.
        last = floor((current+0.00001) / step) * step;
    }
}

