/* 
 * File:   Approximation.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 10:21 AM
 */

#ifndef APPROXIMATION_H
#define	APPROXIMATION_H

#include <iostream>
#include <stdio.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstring>

#include "base.h"
#include "ICipher.h"

// Maximum order of terms.
#define MAX_ORDER 7

class Approximation {
private:
    ICipher  * cip;
    
    // Key & input parameters prepared for cipher input.
    uchar * finput;
    
public:
    Approximation();
    virtual ~Approximation();
    
    void work();
    
    ICipher * getCipher(){ return cip; }
    void      setCipher(ICipher * cip) { this->cip = cip; }
    
    void genMessages();
};

#endif	/* APPROXIMATION_H */

