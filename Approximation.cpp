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

Approximation::Approximation() : cip(NULL), finput(NULL), 
        ulongInp(NULL), ulongOut(NULL), dumpCoefsToFile(false), 
        outputWidthUlong(0), inputWidthUlong(0), cubeBinomialSums(NULL) {
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
    
    if (ulongOut!=NULL){
        delete[] ulongOut;
        ulongOut = NULL;
    }
    
    if (ulongInp!=NULL){
        delete[] ulongInp;
        ulongInp = NULL;
    }
}

void Approximation::setCipher(ICipher* cip) {
    this->cip = cip;
    this->byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    this->outputWidthUlong = OWN_CEIL((double)cip->getOutputBlockSize() / (double)SIZEOF_ULONG);
    this->inputWidthUlong  = OWN_CEIL((double)byteWidth / (double)SIZEOF_ULONG);
}

void Approximation::init() {
    assert(cip!=NULL);
    assert(outputWidthUlong>0 && inputWidthUlong>0);
    
    // Pre-compute cube binomial coefficients
    cubeBinomialSums = new ULONG[8*byteWidth];
    cubeBinomialSums[0] = 0ul;
    
    // Allocating & computing binomial sums for cubic term coefficients.
    ULONG res = 0;
    for(uint i = 1; i<8*byteWidth; i++){
        res += CombinatiorialGenerator::binomial(8*byteWidth-i, 2);
        cubeBinomialSums[i] = res;
    }
    
    // Allocating space for the coefficients, for each output polynomial we
    // allocate separate coefficients vector.
    for(unsigned int order = 0; order<=orderLimit; order++){
        // Compute the size of coefficient array.
        // Binomial(8*byteWidth, order) * outputWidthUlong.
        //
        ULONG vecSize = CombinatiorialGenerator::binomial(8*byteWidth, order) * outputWidthUlong;
        coefficients[order].assign(vecSize, (ULONG) 0);
    }
    
    ulongOut = new ULONG[outputWidthUlong];
    ulongInp = new ULONG[inputWidthUlong];
}

ULONG Approximation::getCubeIdx(ULONG x1, ULONG x2, ULONG x3) {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    return cubeBinomialSums[x1] + CombinatiorialGenerator::getQuadIdx(byteWidth-1-x1, x2, x3);
}

void Approximation::work() {    
    // Boundary on the term order to store term coefficients.
    orderLimit=3;
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    finput         = new uchar[byteWidth];
    uchar * key    = finput + cip->getInputBlockSize();
    
    // Further pre-computation & initialization.
    init();
    
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
    for(unsigned int i = 0; i < cip->getOutputBlockSize(); i++){
        coefficients[0][i/SIZEOF_ULONG] = READ_TERM_1(coefficients[0][i/SIZEOF_ULONG], output[i], i%SIZEOF_ULONG);
        
        // File dumping is temporarily disabled due to representation switch.
        //(*coefs[i]) << ((uint)coefficients[0][i][0]) << endl;
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
        
        ProgressMonitor pm(0.01);
        ofstream cip1(std::string("ciphertexts_order_") + std::to_string(order) + ".txt");
        for(; cg.next(); ){
            const uchar * input = cg.getCurCombination();
            
            // Evaluate cipher on current combinations.
            cip->evaluate(input, input + cip->getInputBlockSize(), output);
            
            // Transform output to the ULONG array
            // For better memory handling and XORing in an one big register.
            memset(ulongOut, 0, sizeof(ULONG) * outputWidthUlong);
            for(uint x=0; x<cip->getOutputBlockSize(); x++){
                ulongOut[x/SIZEOF_ULONG] = READ_TERM_1(ulongOut[x/SIZEOF_ULONG], output[x], x%SIZEOF_ULONG);
            }
            
            // Dump
            //dumpUchar(cip1, output, cip->getOutputBlockSize());
            
            // Generate coefficients of this order, perform it simultaneously
            // on one ULONG type.
            for(uint ulongCtr=0; ulongCtr<outputWidthUlong; ulongCtr++){
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
                ULONG curValue = ulongOut[ulongCtr];
                
                // 1. XOR constant
                curValue ^= coefficients[0][ulongCtr];
                
                // 2. XOR monomials included in current combination, if applicable.
                for(uint ti=0; order>1 && orderLimit>=1 && ti<order; ti++){
                    curValue ^= coefficients[1][outputWidthUlong*cg.getCurState()[ti] + ulongCtr];
                }
                
                // 3. XOR all quadratic terms, if applicable.
                // Using combinations generator to generate all quadratic terms 
                // that are possible to construct using variables present in current term;
                if (order>2 && orderLimit>=2){
                    for(cgQuadratic.reset(); cgQuadratic.next(); ){
                        ULONG idx = CombinatiorialGenerator::getQuadIdx(
                                8*byteWidth, 
                                cg.getCurState()[cgQuadratic.getCurState()[0]],  // first var. idx. in quadr. term. 
                                cg.getCurState()[cgQuadratic.getCurState()[1]]); // second var. idx. in quadr. term.
                        
                        curValue ^= coefficients[2][outputWidthUlong*idx + ulongCtr];
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
                        
                        curValue ^= coefficients[3][outputWidthUlong*idx + ulongCtr];
                    }
                }
                
                // Store value of the current coefficient on his place in the coef. vector.
                coefficients[order][outputWidthUlong*cg.getCounter() + ulongCtr] = curValue;
                
                // If term is too high, cannot continue since we don not have lower terms stored.
                if (order>=orderLimit+1){
                    break;
                }
                
                // Dump terms to the file, only if term is present.
                // File dumping is temporarily disabled due to representation switch.
                /*if (dumpCoefsToFile && order<=orderLimit+1 && curValue){
                    for(uint ti=0; ti<order; ti++){
                        (*coefs[i]) << "x_" << setw(4) << setfill('0') << right << (uint)(*(cg.getCurState()+ti));
                    }
                    
                    (*coefs[i]) << " + ";
                }*/
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
    
    // Test approximation correctness
    cout << "Testing approximation correctness: " << endl << " ";
    selftestApproximation();
    
    // Here we test the accuracy of the high order approximation, random key,
    // random message. 
    cout << "Testing approximation quality: " << endl << " ";
    testPolynomialApproximation();
    
    delete[] output;    
    cout << "Generating finished" << endl;
}

int Approximation::selftestApproximation() {
// Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * outputPol = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    ULONG * hits   = new ULONG[8*cip->getOutputBlockSize()];
    const ULONG genLimit = 100ul;
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    memset(hits, 0, sizeof(ULONG)*8*cip->getOutputBlockSize());
    
    // Generate messages and keys, evaluate it both on cipher and polynomials.
    ProgressMonitor pm(0.01);
    for(unsigned long i=0; i<genLimit; i++){
        // Generate cipher input at random.
        memset(input, 0, sizeof(uchar) * byteWidth);
        
        // Generate input block randomly with hamming weight smaller than maximal
        // precomputed order.
        uint randOrder = 1 + (rand() % orderLimit);
        for(uint k=0; k<randOrder; k++){
            const uint randIdx = rand() % (8*byteWidth);
            input[(randIdx/8) % byteWidth] |= 1ul << (randIdx%8);
        }
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(input, outputPol);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            hits[p] += (outputCip[p/8] & (1u << (p%8))) == (outputPol[p/8] & (1u << (p%8)));
        }
        
        // Progress monitoring.
        double cProg = (double)i / (double)genLimit;
        pm.setCur(cProg);
    }
    pm.setCur(1.0);
    
    // Determine test success
    bool success=true;
    for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
        if (hits[p]!=genLimit){
            success=false;
            break;
        }
    }
    
    cout << endl << "Self test finished: ";
    if (success){
        cout << " [  OK  ]" << endl;
    } else {
        cout << " [ FAIL ]" << endl << "Frequencies: " << endl;
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            cout << dec << "  f_" << setw(4) << setfill('0') << right << p << " = ";
            cout << ((double)hits[p] / (double)genLimit) << endl;
        }
    }
    
    // Free the memory.
    delete[] hits;
    delete[] outputCip;
    delete[] outputPol;
    delete[] input;
    return success;
}

int Approximation::testPolynomialApproximation() {
    // Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * outputPol = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    ULONG * hits   = new ULONG[8*cip->getOutputBlockSize()];
    const ULONG genLimit = 1000ul;
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    memset(hits, 0, sizeof(ULONG)*8*cip->getOutputBlockSize());
    
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
        
        //cout << "final: " << endl;
        //dumpUcharHex(cout, outputCip, 16);
        //dumpUcharHex(cout, outputPol, 16);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            hits[p] += (outputCip[p/8] & (1u << (p%8))) == (outputPol[p/8] & (1u << (p%8)));
        }
        
        // Progress monitoring.
        double cProg = (double)i / (double)genLimit;
        pm.setCur(cProg);
    }
    pm.setCur(1.0);
    
    cout << endl << "Approximation quality test finished." << endl;
    for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
        cout << dec << "  f_" << setw(4) << setfill('0') << right << p << " = ";
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
    const uint bitWidth = 8*byteWidth;
    
    // Reset output buffer, only ones will be set here, has to be set to zero.
    memset(ulongInp, 0, sizeof(ULONG) * inputWidthUlong);
    
    // Copy input to ULONG for better manipulation.
    for(uint x=0; x<byteWidth; x++){
        ulongInp[x/SIZEOF_ULONG] = READ_TERM_1(ulongInp[x/SIZEOF_ULONG], input[x], x%SIZEOF_ULONG);
    }
    
    // Evaluation on ULONGs.
    for(uint ulongCtr=0; ulongCtr<outputWidthUlong; ulongCtr++){
        // Evaluate SIZEOF_ULONG polynomials simultaneously.
        // 1. Use constant term for initialization.
        ulongOut[ulongCtr] = coefficients[0][ulongCtr];
    }
    
    // 2. linear, quadratic and cubic terms, quartic and higher if applicable.
    for(uint order=1; order<=orderLimit; order++){
        CombinatiorialGenerator cgen(bitWidth, order);
        for(; cgen.next(); ){
            const ULONG ctr = cgen.getCounter();
            
            // Now term being evaluated is fixed, defined by the state
            // of the combinatorial generator. 
            //
            // Get bit-mask with those bits enabled corresponding to variables in
            // the particular term determined by cgen.
            const ULONG * comb  = cgen.getCurUlongCombination();
            
            //
            // Evaluate particular term on the input.
            //
            // Some elements of the array can be zero, but this does not mean
            // the term itself is zero, it just can be defined on the end of the
            // array. 
            //
            // For example x_127 in 64-bit architecture would be comb[0]=0, comb[1]=highest bit.
            //
            bool termEval=true;
            for(uint uctr2=0; uctr2<inputWidthUlong; uctr2++){
                termEval &= (comb[uctr2]==0) ? 1 : (comb[uctr2] & ulongInp[uctr2]) == comb[uctr2];
            }
            
            // If term is null, nothing to do here, go evaluate next one.
            if (!termEval){
                continue;
            }
            
            // Term is evaluated to 1, thus XOR it to the result - where it is present.
            for(uint uctr2=0; uctr2<outputWidthUlong; uctr2++){
                ulongOut[uctr2] ^= coefficients[order][outputWidthUlong*ctr + uctr2];
            }
        }
    }
    
    // Transform ULONG to output.
    for(uint x=0; x<cip->getOutputBlockSize(); x++){
        output[x] = (ulongOut[x/SIZEOF_ULONG] >> (8* (x % SIZEOF_ULONG))) & ((unsigned char)0xffu);
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

