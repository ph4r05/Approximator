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

Approximation::Approximation() : cip(NULL), finput(NULL) {
}

Approximation::~Approximation() {
}

void Approximation::work() {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    
    // Boundary on the term order to store term coefficients.
    const uint order2store=3;
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    finput         = new uchar[byteWidth];
    uchar * key    = finput + cip->getInputBlockSize();
    
    // Allocating space for the coefficients.
    for(unsigned int i = 0; i<4; i++){
        coefficients[i] = new std::vector<bool>[8 * cip->getOutputBlockSize()];
    }
    
    // Generate ciphertext for constant    
    memset(finput, 0, byteWidth);
    cip->evaluate(finput, key, output);
    
    // Dump
    ofstream hist("cip_0.txt");
    dumpUchar(hist, output, cip->getOutputBlockSize());
    hist.close();
    
    // Dump polynomials to the files, one file per polynomial.
    ofstream ** coefs = new ofstream*[8 * cip->getOutputBlockSize()];
    for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
        coefs[i] = new ofstream(std::string("poly_") + std::to_string(i) + ".txt");
    }
    
    // Read output to constants.
    for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
        coefficients[0][i].resize(1, 0);
        coefficients[0][i].at(0) = (output[i/8] & (1u<<(i%8))) > 0;
        (*coefs[i]) << ((uint)coefficients[0][i].at(0)) << endl;
    }
    
    // Generate order1 .. order3 cipher data
    for(unsigned order=1; order<=order2store; order++){
        CombinatiorialGenerator cg(byteWidth*8, order);
        
        cout << "Starting with order: " 
                << order 
                << "; number of bytes to store coefficients: " 
                << (8 * cip->getOutputBlockSize() * order * cg.getTotalNum() / 8)
                << endl;
        
        // Generate coefficients for this order & store it to the vector.
        for(unsigned int i = 0; order<=3 && i < 8 * cip->getOutputBlockSize(); i++){
            coefficients[order][i].resize(order*(cg.getTotalNum()+1), 0);
        }
        
        ofstream cip1(std::string("ciphertexts_order_") + std::to_string(order) + ".txt");
        for(ULONG ctr=0; cg.next(); ctr++){
            uchar * input = cg.getCurCombination();
            
            // Evaluate cipher on current combinations.
            cip->evaluate(input, input + cip->getInputBlockSize(), output);
            
            // Dump
            dumpUchar(cip1, output, cip->getOutputBlockSize());
            
            // Generate coefficients for this order & store it to the vector.
            for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
                // Evaluate coefficients for each polynomial in the cipher
                // for the given term specified by the state of the combinatorial generator. 
                //
                // Current term value: all previous terms including enabled bits XOR ciphertext.
                // For example, term to determine: x1x6x9:
                //   constant                XOR
                //   x1   XOR x6   XOR x9    XOR
                //   x1x6 XOR x1x9 XOR x6x9
                //
                // Obtain ciphertext bit corresponding to the current polynomial. 
                bool curValue = (output[i/8] & (1u<<(i%8))) > 0;
                // XOR constant
                curValue ^= coefficients[0][i].at(0);
                // XOR monomials included in current combination, if applicable.
                for(uint ti=0; order>1 && ti<order; ti++){
                    curValue ^= coefficients[1][i].at(cg.getCurState()[ti]);
                }
                // XOR all quadratic terms, if applicable.
                // Using combinations generator to generate all quadratic terms 
                // that are possible to construct using variables present in current term;
                if (order>2){
                    CombinatiorialGenerator cgQuadratic(order, 2);
                    for(; cgQuadratic.next(); ){
                        ULONG idx = CombinatiorialGenerator::getQuadIdx(
                                8*byteWidth, 
                                cg.getCurState()[cgQuadratic.getCurState()[0]],  // first var. idx. in quadr. term. 
                                cg.getCurState()[cgQuadratic.getCurState()[1]]); // second var. idx. in quadr. term.
                        
                        curValue ^= coefficients[2][i].at(idx);
                    }
                }
                // XOR all cubic terms, if applicable.
                // Using combinations generator to generate all cubic terms 
                // that are possible to construct using variables present in current term;
                if (order>3 && order2store>=3){
                    CombinatiorialGenerator cgCubic(order, 3);
                    for(; cgCubic.next(); ){
                        ULONG idx = CombinatiorialGenerator::getCubeIdx(
                                8*byteWidth, 
                                cg.getCurState()[cgCubic.getCurState()[0]],  // first var. idx. in quadr. term. 
                                cg.getCurState()[cgCubic.getCurState()[1]],  // second var. idx. in quadr. term.
                                cg.getCurState()[cgCubic.getCurState()[2]]); // third var. idx. in quadr. term.
                        
                        curValue ^= coefficients[3][i].at(idx);
                    }
                }
                
                // Store only for limited level...
                for(uint ti=0; order<=order2store && ti<order; ti++){
                    coefficients[order][i].at(order*cg.getCounter() + ti) = curValue;
                }
                
                // If term is too high, cannot continue since we don not have lower terms stored.
                if (order>=order2store+1){
                    break;
                }
                
                // Dump terms to the file, only if term is present.
                if (order<=order2store+1 && curValue){
                    for(uint ti=0; ti<order; ti++){
                        (*coefs[i]) << "x_" << setw(4) << setfill('0') << right << (uint)(*(cg.getCurState()+ti));
                    }
                    
                    (*coefs[i]) << " + ";
                }
            }
        }
        
        cip1.close();
        
        // Make a newline at the end of the current order.
        for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
            (*coefs[i]) << endl;
        }
    }
    
    // Free memory allocated for coefficients.
    for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
        coefs[i]->close();
        delete coefs[i];
        coefs[i] = NULL;
    }
    
    delete[] coefs;
    
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

