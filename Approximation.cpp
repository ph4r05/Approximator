/* 
 * File:   Approximation.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 10:21 AM
 */

#include "Approximation.h"
#include "CombinatiorialGenerator.h"
#include "ProgressMonitor.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include <cstring>
#include <cstdlib> 
#include <ctime> 
#include <cassert>
#include <unistd.h>

using namespace std;

Approximation::Approximation() : cip(NULL), finput(NULL), dumpCoefsToFile(false), cubeBinomialSums(NULL) {
}

Approximation::~Approximation() {
    if (finput!=NULL){
        delete[] finput;
        finput=NULL;
    }
    
    if (cubeBinomialSums!=NULL){
        delete[] cubeBinomialSums;
        cubeBinomialSums=NULL;
    }
}

void Approximation::setCipher(ICipher* cip) {
    this->cip = cip;
    this->byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
}


void Approximation::init() {
    assert(cip!=NULL);
    
    // Pre-compute cube binomial coefficients
    cubeBinomialSums = new ULONG[8*byteWidth];
    cubeBinomialSums[0] = 0ul;
    
    ULONG res = 0;
    for(uint i = 1; i<8*byteWidth; i++){
        res += CombinatiorialGenerator::binomial(8*byteWidth-i, 2);
        cubeBinomialSums[i] = res;
    }
}

ULONG Approximation::getCubeIdx(ULONG x1, ULONG x2, ULONG x3) {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    return cubeBinomialSums[x1] + CombinatiorialGenerator::getQuadIdx(byteWidth-1-x1, x2, x3);
}

void Approximation::work() {    
    // Boundary on the term order to store term coefficients.
    orderLimit=2;
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    finput         = new uchar[byteWidth];
    uchar * key    = finput + cip->getInputBlockSize();
    
    // Further pre-computation & initialization.
    init();
    
    // Allocating space for the coefficients, for each output polynomial we
    // allocate separate coefficients vector.
    for(unsigned int i = 0; i<4; i++){
        coefficients[i] = new std::vector<bool>[8 * cip->getOutputBlockSize()];
    }
    
    // Generate ciphertext for calculating constant term (key=0, message=0).
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
    
    // Read output of encryption and obtain constant terms.
    for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
        coefficients[0][i].resize(1, 0);
        coefficients[0][i][0] = (output[i/8] & (1u<<(i%8))) > 0;
        (*coefs[i]) << ((uint)coefficients[0][i][0]) << endl;
    }
    
    // Generate order1 .. order3 cipher data
    for(unsigned order=1; order<=orderLimit; order++){
        CombinatiorialGenerator cg(byteWidth*8, order);
        CombinatiorialGenerator cgQuadratic(order > 2 ? order : 2, 2);
        CombinatiorialGenerator cgCubic(order > 3 ? order : 3, 3);
        
        cout << "Starting with order: " 
                << order 
                << "; combinations: " << cg.getTotalNum()
                << "; number of bytes to store coefficients: " 
                << (8 * cip->getOutputBlockSize() * cg.getTotalNum() / 8)
                << endl;
        cout << " ";
        
        // Generate coefficients for this order & store it to the vector.
        for(unsigned int i = 0; order<=3 && i < 8 * cip->getOutputBlockSize(); i++){
            coefficients[order][i].resize(cg.getTotalNum()+1, 0);
        }
        
        ProgressMonitor pm(0.01);
        ofstream cip1(std::string("ciphertexts_order_") + std::to_string(order) + ".txt");
        for(; cg.next(); ){
            const uchar * input = cg.getCurCombination();
            
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
                
                // 1. XOR constant
                curValue ^= coefficients[0][i][0];
                
                // 2. XOR monomials included in current combination, if applicable.
                for(uint ti=0; order>1 && ti<order; ti++){
                    curValue ^= coefficients[1][i][cg.getCurState()[ti]];
                }
                
                // 3. XOR all quadratic terms, if applicable.
                // Using combinations generator to generate all quadratic terms 
                // that are possible to construct using variables present in current term;
                if (order>2){
                    for(cgQuadratic.reset(); cgQuadratic.next(); ){
                        ULONG idx = CombinatiorialGenerator::getQuadIdx(
                                8*byteWidth, 
                                cg.getCurState()[cgQuadratic.getCurState()[0]],  // first var. idx. in quadr. term. 
                                cg.getCurState()[cgQuadratic.getCurState()[1]]); // second var. idx. in quadr. term.
                        
                        curValue ^= coefficients[2][i][idx];
                    }
                }
                
                // 4. XOR all cubic terms, if applicable.
                // Using combinations generator to generate all cubic terms 
                // that are possible to construct using variables present in current term;
                if (order>3 && orderLimit>=3){
                    for(cgCubic.reset(); cgCubic.next(); ){
                        ULONG idx = getCubeIdx(
                                cg.getCurState()[cgCubic.getCurState()[0]],  // first var. idx. in quadr. term. 
                                cg.getCurState()[cgCubic.getCurState()[1]],  // second var. idx. in quadr. term.
                                cg.getCurState()[cgCubic.getCurState()[2]]); // third var. idx. in quadr. term.
                        
                        curValue ^= coefficients[3][i][idx];
                    }
                }
                
                // Store value of the current coefficient on his place in the coef. vector.
                coefficients[order][i][cg.getCounter()] = curValue;
                
                // If term is too high, cannot continue since we don not have lower terms stored.
                if (order>=orderLimit+1){
                    break;
                }
                
                // Dump terms to the file, only if term is present.
                if (dumpCoefsToFile && order<=orderLimit+1 && curValue){
                    for(uint ti=0; ti<order; ti++){
                        (*coefs[i]) << "x_" << setw(4) << setfill('0') << right << (uint)(*(cg.getCurState()+ti));
                    }
                    
                    (*coefs[i]) << " + ";
                }
            }
            
            // Progress monitoring.
            double cProg = (double)cg.getCounter() / (double)cg.getTotalNum();
            pm.setCur(cProg);
        }
        pm.setCur(1.0);
        
        cout << endl;
        cip1.close();
        
        // Make a newline at the end of the current order.
        for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
            (*coefs[i]) << endl;
        }
    }
    
    // Here we test the accuracy of the high order approximation, random key,
    // random message. 
    cout << "Testing approximation quality: " << endl << " ";
    sleep(1);
    testPolynomialApproximation();
    
    // Free memory allocated for coefficients.
    for(unsigned int i = 0; i < 8 * cip->getOutputBlockSize(); i++){
        coefs[i]->close();
        delete coefs[i];
        coefs[i] = NULL;
    }
    
    delete[] coefs;
    delete[] output;
    
    cout << "Generating finished" << endl;
}

int Approximation::testPolynomialApproximation() {
     // Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * outputPol = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    ULONG * hits   = new ULONG[8*byteWidth];
    const ULONG genLimit = 100ul;
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    for(uint p=0; p<8*byteWidth; p++){
        hits[p] = 0;
    }
    
    // Generate 2^20 random messages and keys, evaluate it
    // both on cipher and polynomials.
    ProgressMonitor pm(0.01);
    for(unsigned long i=0; i<genLimit; i++){
        // Generate cipher input at random.
        for(unsigned int k=0; k<byteWidth; k++){ 
            input[k] = (rand() % (0xffu+1)); 
        }
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(input, outputPol);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*byteWidth; p++){
            hits[p] += (outputCip[p/8] & (1u << (p%8))) == (outputPol[p/8] & (1u << (p%8)));
        }
        
        // Progress monitoring.
        double cProg = (double)i / (double)genLimit;
        pm.setCur(cProg);
    }
    pm.setCur(1.0);
    
    cout << endl << "Approximation quality test finished." << endl;
    for(uint p=0; p<8*byteWidth; p++){
        cout << "  f_" << setw(4) << setfill('0') << right << p << " = ";
        cout << ((double)hits[p] / (double)genLimit) << endl;
    }
    
    // Free the memory.
    delete[] hits;
    delete[] outputCip;
    delete[] outputPol;
    delete[] input;
    
    return 1;
}


int Approximation::evaluateCoefficients(const unsigned char* input, unsigned char* output) {
    // We can assume that approximate half of the coefficients are enabled/present
    // in the resulting polynomial, thus evaluation is based on the iteration of 
    // the combinatorial generator and reading coefficient by coefficient.
    const uint numPolynomials = 8*cip->getOutputBlockSize();
    const uint bitWidth = 8*byteWidth;
    
    // Reset output buffer, only ones will be set here, has to be set to zero.
    memset(output, 0, cip->getOutputBlockSize());
    
    for(uint pidx=0; pidx < numPolynomials; pidx++){
        // Evaluate polynomial pidx on the input provided.
        
        // 1. Use constant term for initialization.
        bool res = coefficients[0][pidx][0];
        
        // 2. linear terms
        for(uint j=0; j<bitWidth; j++){
            res ^= coefficients[1][pidx][j] & ((input[j/8] & (1u << (j%8))) > 0);
        }
        
        // 3. quadratic and cubic terms, quartic and higher if applicable.
        for(uint order=2; order<=orderLimit; order++){
            CombinatiorialGenerator cgen(bitWidth, order);
            for(; cgen.next(); ){
                const ULONG ctr = cgen.getCounter();
                const ULONG * state = cgen.getCurState();
                
                if (coefficients[order][pidx][ctr]==0) {
                    // This coefficient is zero, term would not have any contribution.
                    continue;
                }
                
                // Coefficient is one, evaluate this term on the input data.
                bool termVal=1;
                
                // State is array of size $order, in each element it contains 
                // index of the variable present in the term.
                // Example: the first state for order=3 is {0,1,2}
                // What corresponds to x_0x_1x_2, thus bits 0,1,2 has to be multiplied.
                for(uint i=0; i<order && termVal; i++){
                    termVal &= ((input[state[i]/8] & (1u << (state[i]%8))) > 0);
                }
                
                res ^= termVal;
            }
        }
        
        // Store result to the output buffer.
        if (res){
            output[pidx/8] |= (1ul << (pidx%8));
        }
    }
    
    return 0;
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

