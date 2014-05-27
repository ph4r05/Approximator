/* 
 * File:   Approximation.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 10:21 AM
 * TODO: multithreaded!!!
 * TODO: MPI???
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
#include <NTL/mat_GF2.h>
#include <openssl/md5.h>
#include "NTLUtils.h"
#include <boost/unordered_map.hpp>


// The following macro should be 1 to call FGb modulo a prime number.
#define LIBMODE 1  
#define CALL_FGB_DO_NOT_DEFINE
#include "call_fgb.h"
#define FGb_MAXI_BASE 100000

NTL_CLIENT
using namespace std;
using namespace NTL;

Approximation::Approximation(uint orderLimit) : cip(NULL), dumpCoefsToFile(false), 
        outputWidthUlong(0), inputWidthUlong(0), binomialSums(NULL),
        threadCount(1), varNames(NULL), keybitsToZero(0),
        poly2take(NULL), numPolyActive(0), fgbFile(NULL) {
    
    this->orderLimit = orderLimit;
}

Approximation::~Approximation() {    
    if (binomialSums!=NULL){
        for(uint order = 1; order <= orderLimit; order++){
            delete[] binomialSums[order-1];
            binomialSums[order-1]=NULL;
        }
        
        delete[] binomialSums;
        binomialSums = NULL;
    }
    
    if (varNames!=NULL){
        for(uint var = 0; var <= 8*byteWidth; var++){
            delete[] varNames[var];
            varNames[var]=NULL;
        }
        
        delete[] varNames;
        varNames = NULL;
    }
    
    if (poly2take!=NULL){
        delete[] poly2take;
        poly2take = NULL;
    }
    
    if (fgbFile!=NULL){
        fclose(fgbFile);
        fgbFile = NULL;
    }
}

void Approximation::setCipher(ICipher* cip) {
    this->cip = cip;
    this->byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    this->logBitInputWidth = ceil(log2(this->byteWidth*8));
    this->outputWidthUlong = OWN_CEIL((double)cip->getOutputBlockSize() / (double)SIZEOF_ULONG);
    this->inputWidthUlong  = OWN_CEIL((double)byteWidth / (double)SIZEOF_ULONG);
}

void Approximation::init() {
    assert(cip!=NULL);
    assert(outputWidthUlong>0 && inputWidthUlong>0);
    assert(orderLimit>=0 && orderLimit<MAX_ORDER);
    assert(SIZEOF_ULONG == sizeof(ULONG));
    
    // Pre-compute binomial sums for combination index computation.
    if (orderLimit >= 1){
        binomialSums = new ULONG * [orderLimit];
        for(uint order = 1; order <= orderLimit; order++){
            const uint idx = order-1;
            binomialSums[idx]    = new ULONG[8*byteWidth];
            binomialSums[idx][0] = 0ul;

            ULONG res = 0;
            for(uint i = 1; i<8*byteWidth; i++){
                res += CombinatiorialGenerator::binomial(8*byteWidth-i, order-1);
                binomialSums[idx][i] = res;
            }
        }
    }
    
    // Variable names for FGb.
    varNames = new char * [8*byteWidth];
    for(uint var=0; var<8*byteWidth; var++){
        varNames[var] = new char[10];
        snprintf(varNames[var], 10, "x_%03d", var);
    }
    
    // Init polynomial-to-take bitmap.
    poly2take = new ULONG[outputWidthUlong];
    numPolyActive = 8*cip->getOutputBlockSize();
    memset(poly2take, 0xff, SIZEOF_ULONG * outputWidthUlong);
    
    // Open logfile for library
    fgbFile = fopen("fgb.log", "a+");
}

ULONG Approximation::getCubeIdx(ULONG x1, ULONG x2, ULONG x3) const {
    return orderLimit < 3 ? 0 : binomialSums[2][x1] + CombinatiorialGenerator::getQuadIdx(8*byteWidth-1-x1, x2-x1-1, x3-x1-1);
}

ULONG Approximation::getCombinationIdx(uint order, const ULONG* xs, uint xsOffset, ULONG Noffset, ULONG combOffset) const {
    assert(order<=orderLimit);
    
    // Small order.
    if (order==0) return 0;
    if (order==1) return xs[0+xsOffset] - combOffset;
    if (order==2) return CombinatiorialGenerator::getQuadIdx(8*byteWidth-Noffset, xs[0+xsOffset]-combOffset, xs[1+xsOffset]-combOffset);
    
    // Order 3.
    // Take recursive parameters into consideration.
    const ULONG x1 = xs[0+xsOffset] - combOffset;
    
    // This is simple optimization to remove 1 recursion step.
    if (order==3) {
        const ULONG sum = binomialSums[2][x1+Noffset] - binomialSums[2][Noffset];
        return sum + CombinatiorialGenerator::getQuadIdx(
                8*byteWidth-1-x1-Noffset, 
                xs[1+xsOffset]-combOffset-x1-1, 
                xs[2+xsOffset]-combOffset-x1-1);
    }
    
    // Order 3 and higher.
    // Sum of the previous combinations is shifted due to Noffset.
    // N is reduced thus the binomial sum is shifted (originating point is smaller).
    return (binomialSums[order-1][x1+Noffset] - binomialSums[order-1][Noffset]) + 
            getCombinationIdx(
            order-1, 
            xs, 
            xsOffset+1, 
            Noffset+1+x1, 
            combOffset+1+x1
           );
}

int Approximation::getCombinationFromIdx(uint order, ULONG* xs, ULONG idx) const {    
    int nOffset=0;
    for(int x=order-1; x>=1; x--){
        int i=8*byteWidth-1-nOffset;
        for(; i>=0; i--){
            const ULONG biSum = (binomialSums[x][i+nOffset] - binomialSums[x][nOffset]);
            // If current index is bigger than the current sum, the previous sum
            // helps us to determine index of the order.
            if (idx >= biSum){
                break;
            }
        }

        idx-=(binomialSums[x][i+nOffset] - binomialSums[x][nOffset]);
        
        xs[order-x-1] = i+nOffset;
        nOffset = xs[order-x-1]+1;
    }
    
    // Final element is already determined.
    xs[order-1] = idx + nOffset;
    return 1;
}

ULONG Approximation::getCombinationULong(uint order, const ULONG* xs) const {
    assert(order <= 0xf);
    ULONG res = order & 0xf;
    
    uint offset = 4;  // order.
    for(uint x=0; x<order; x++){
        res |= xs[x] << offset;
        offset += logBitInputWidth;
    }
    
    return res;
}

int Approximation::getCombinationFromULong(ULONG* xs, ULONG combUlong) const {
    uint order = combUlong & 0x7;
    assert(order <= MAX_ORDER);
    
    combUlong = combUlong >> 4; // remove order.
    for(uint x=0; x<order; x++){
        xs[x] = combUlong & ((1u << logBitInputWidth)-1);
        combUlong = combUlong >> logBitInputWidth;
    }
    
    return 1;
}

bool Approximation::isPoly2Take(uint polyIdx) const {
    return (poly2take[polyIdx / (8*SIZEOF_ULONG)] & (ULONG1 << (polyIdx % (8*SIZEOF_ULONG)))) > 0;
}

void Approximation::setPoly2Take(const std::vector<std::string> & map) {
    // Set by default to zero.
    memset(poly2take, 0x0, SIZEOF_ULONG * outputWidthUlong);
    // And iterate over enabled functions. 
    numPolyActive=0;
    for (std::vector<std::string>::const_iterator it = map.begin() ; it != map.end(); ++it){
        const std::string cur = *it;
        const int idx = std::stoi(cur);
        // Possible duplicates.
        if (poly2take[idx / (8*SIZEOF_ULONG)] & ULONG1 << (idx % (8*SIZEOF_ULONG))) continue;
        // Set to map & update counter.
        poly2take[idx / (8*SIZEOF_ULONG)] |= ULONG1 << (idx % (8*SIZEOF_ULONG));
        numPolyActive+=1;
    }
}

uint Approximation::getNumVariables() const {
    return cip->getKeyBlockSize()*8-keybitsToZero;
}

void Approximation::computeCoefficients() {
    ULONG * ulongOut = new ULONG[outputWidthUlong];
    ULONG * ulongInp = new ULONG[inputWidthUlong];
    uchar * finput   = new uchar[byteWidth];
    CombinatiorialGenerator ** cgenerators = new CombinatiorialGenerator * [orderLimit];
    ULONG * tmpCombination = new ULONG[orderLimit+1];
    
    // Allocating space for the coefficients.
    for(unsigned int order = 0; order<=orderLimit; order++){
        // Compute the size of coefficient array.
        // Binomial(8*byteWidth, order) * outputWidthUlong.
        //
        ULONG vecSize = CombinatiorialGenerator::binomial(8*byteWidth, order) * outputWidthUlong;
        cout << "  Allocating coefficient storage for order " << order << "; Bytes=" << vecSize * sizeof(ULONG) << endl;
        coefficients[order].assign(vecSize, (ULONG) 0);
    }
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    uchar * key    = finput + cip->getInputBlockSize();
    
    // Generate ciphertext for calculating constant term (key=0, message=0).
    memset(finput, 0, byteWidth);
    cip->evaluate(finput, key, output);
    
    // Read output of encryption and obtain constant terms.
    for(unsigned int i = 0; i < cip->getOutputBlockSize(); i++){
        coefficients[0][i/SIZEOF_ULONG] = READ_TERM_1(coefficients[0][i/SIZEOF_ULONG], output[i], i%SIZEOF_ULONG);
    }
    
    // Find polynomial terms coefficient for order 1..orderLimit.
    for(uint order=1; order<=orderLimit; order++){
        CombinatiorialGenerator cg(byteWidth*8, order);
        
        // Combinatorial generators for computing XOR indices.
        for(uint i=0; (i+2) < order; i++){
            cgenerators[i] = new CombinatiorialGenerator(order, i+2);
        }
        
        //CombinatiorialGenerator cgQuadratic(order > 2 ? order : 2, 2);
        //CombinatiorialGenerator cgCubic(order > 3 ? order : 3, 3);
        
        cout << "Starting with order: " 
                << order 
                << "; combinations: " << cg.getTotalNum()
                << "; number of bytes to store coefficients: " 
                << (8 * cip->getOutputBlockSize() * cg.getTotalNum() / 8)
                << endl;
        cout << " ";
        
        // Here is the point for parallelization.
        // Each thread/computing node can compute x-th combination from the generator
        // effectively partitioning combination space. 
        // Synchronization barrier is needed after finishing particular order
        // because in order to compute order N we need to have coefficients of
        // terms of order N-1 and less.
        ProgressMonitor pm(0.01);
        for(; cg.next(); ){
            const uchar * input = cg.getCurCombination();
            
            // Evaluate cipher on current combination.
            cip->evaluate(input, input + cip->getInputBlockSize(), output);
            
            // Transform output to the ULONG array
            // For better memory handling and XORing in an one big register.
            this->readUcharToUlong(output, cip->getOutputBlockSize(), ulongOut);
            
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
                
                // 3. XOR all quadratic, cubic, quartic, etc.. terms if applicable.
                // In order to determine current term coefficient all lower terms
                // that can be obtained from this one has to be taken into account.
                // and XORed into the result.
                for(uint xorOrder=2; xorOrder<order && order<=orderLimit; xorOrder++){
                    
                    // Obtain a shortcut reference to the combinatorial generator
                    // for this xorOrder. This generator is used to generate 
                    // combinations of variables from to original term to construct
                    // a lower term.
                    CombinatiorialGenerator * const xorOrderGen = cgenerators[xorOrder-2];
                    
                    // Iterate over all possible variable combinations to the new 
                    // resulting term of a lower order.
                    for(xorOrderGen->reset(); xorOrderGen->next(); ){
                        
                        // Construct xorOrder term representation for term
                        // index computation.
                        for(uint tmpCombCtr = 0; tmpCombCtr<xorOrder; tmpCombCtr++){
                            tmpCombination[tmpCombCtr] = cg.getCurState()[xorOrderGen->getCurState()[tmpCombCtr]];
                        }
                        
                        // Term index computation.
                        ULONG idx = getCombinationIdx(xorOrder, tmpCombination);
                        
                        // XOR current result register with coefficient register.
                        // If low order term is present in the approximation function,
                        // it has to be taken into account.
                        curValue ^= coefficients[xorOrder][outputWidthUlong*idx + ulongCtr];
                    }
                }
                
                // Store value of the current coefficient on his place in the coef. vector.
                coefficients[order][outputWidthUlong*cg.getCounter() + ulongCtr] = curValue;
                
                // If term is too high, cannot continue since we don not have lower terms stored.
                if (order>=orderLimit+1){
                    break;
                }
            }
            
            // Progress monitoring.
            double cProg = (double)cg.getCounter() / (double)cg.getTotalNum();
            pm.setCur(cProg);
        }
        pm.setCur(1.0);
        
        cout << endl;
        
        // Combinatorial generators destruction.
        for(uint i=0; i < order; i++){
            delete cgenerators[i];
        }
    }
    
    delete[] tmpCombination;
    delete[] ulongInp;
    delete[] ulongOut;
    delete[] output;  
    delete[] finput;
    delete[] cgenerators;
}

int Approximation::selftestApproximation(unsigned long numSamples) const {
// Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * outputPol = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    ULONG * hits   = new ULONG[8*cip->getOutputBlockSize()];
    ULONG * iBuff = new ULONG[this->inputWidthUlong];   // Input ULONG buffer
    ULONG * variablesValueMask = new ULONG[this->inputWidthUlong];
    ULONG * ulongOut = new ULONG[outputWidthUlong];
    ULONG * ulongInp = new ULONG[inputWidthUlong];
    const ULONG genLimit = numSamples;
    uint matchErrors=0;
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    memset(hits, 0, sizeof(ULONG)*8*cip->getOutputBlockSize());
    memset(variablesValueMask, 0xff, sizeof(ULONG)*this->inputWidthUlong);
    
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
            input[(randIdx/8)] |= ULONG1 << (randIdx%8);
        }
        
        this->readUcharToUlong(input, byteWidth, iBuff);
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(input, outputPol, ulongInp, ulongOut);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            hits[p] += (outputCip[p/8] & (ULONG1 << (p%8))) == (outputPol[p/8] & (ULONG1 << (p%8)));
        }
        
        // Test partial evaluation for correctness, has to be in match with full evaluation.
        std::vector<ULONG> coeffEval[MAX_ORDER];
        coeffEval[0].assign(this->inputWidthUlong, 0ul);
        this->partialEvaluation(this->byteWidth*8, variablesValueMask, iBuff, coeffEval);
        
        // Check result of partial evaluation w.r.t. full evaluation. Has to match!
        bool partEvalError=false;
        for(uint x = 0; x < cip->getOutputBlockSize()*8; x++){
            const bool bitFullEval = ((outputPol[x/8] & (ULONG1 << (x%8))) > 0);
            const bool bitPartEval = (((coeffEval[0][x/(8*SIZEOF_ULONG)]) & (ULONG1 << (x % (8*SIZEOF_ULONG)))) > 0);
            if (bitFullEval != bitPartEval){
                partEvalError=true;
            }
        }
        
        if (partEvalError){
            matchErrors+=1;
            //cout << " Error in evaluation!" << endl;
            //dumpBin(cout, coeffEval[0], this->outputWidthUlong);
            //dumpBin(cout, outputPol, cip->getOutputBlockSize());
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
    
    cout << "Matching test:      ";
    if (matchErrors==0){
        cout << " [  OK  ]" << endl;
    } else {
        cout << " [  FAIL  ]   numFails=" << matchErrors << endl;   
    }
    
    // Free the memory.
    delete[] hits;
    delete[] outputCip;
    delete[] outputPol;
    delete[] input;
    delete[] iBuff;
    delete[] variablesValueMask;
    delete[] ulongInp;
    delete[] ulongOut;
    return success;
}

int Approximation::testPolynomialApproximation(unsigned long numSamples) const {
    // Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * outputPol = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    ULONG * hits   = new ULONG[8*cip->getOutputBlockSize()];
    ULONG * ulongOut = new ULONG[outputWidthUlong];
    ULONG * ulongInp = new ULONG[inputWidthUlong];
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    memset(hits, 0, sizeof(ULONG)*8*cip->getOutputBlockSize());
    
    // Generate 2^20 random messages and keys, evaluate it
    // both on cipher and polynomials.
    ProgressMonitor pm(0.01);
    for(unsigned long i=0; i<numSamples; i++){
        // Generate cipher input at random.
        for(unsigned int k=0; k<byteWidth; k++){ 
            input[k] = (rand() % (0xffu+1)); 
        }
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(input, outputPol, ulongInp, ulongOut);
        
        //cout << "final: " << endl;
        //dumpUcharHex(cout, outputCip, 16);
        //dumpUcharHex(cout, outputPol, 16);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            hits[p] += (outputCip[p/8] & (ULONG1 << (p%8))) == (outputPol[p/8] & (ULONG1 << (p%8)));
        }
        
        // Progress monitoring.
        double cProg = (double)i / (double)numSamples;
        pm.setCur(cProg);
    }
    pm.setCur(1.0);
    
    cout << endl << "Approximation quality test finished." << endl;
    for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
        cout << dec << "  f_" << setw(4) << setfill('0') << right << p << " = ";
        cout << ((double)hits[p] / (double)numSamples) << endl;
    }
    
    // Free the memory.
    delete[] hits;
    delete[] outputCip;
    delete[] outputPol;
    delete[] input;
    delete[] ulongInp;
    delete[] ulongOut;
    
    return 1;
}

int Approximation::evaluateCoefficients(const unsigned char* input, unsigned char* output, ULONG * iBuff, ULONG * oBuff) const {
    // We can assume that approximate half of the coefficients are enabled/present
    // in the resulting polynomial, thus evaluation is based on the iteration of 
    // the combinatorial generator and reading coefficient by coefficient.
    const uint bitWidth = 8*byteWidth;
    
    // Reset output buffer, only ones will be set here, has to be set to zero
    // and copy to bigger buffer for better manipulation & speed.
    this->readUcharToUlong(input, byteWidth, iBuff);
    
    // Evaluation on ULONGs.
    for(uint ulongCtr=0; ulongCtr<outputWidthUlong; ulongCtr++){
        // Evaluate SIZEOF_ULONG polynomials simultaneously.
        // 1. Use constant term for initialization.
        oBuff[ulongCtr] = coefficients[0][ulongCtr];
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
                termEval &= (comb[uctr2]==0) ? 1 : (comb[uctr2] & iBuff[uctr2]) == comb[uctr2];
            }
            
            // If term is null, nothing to do here, go evaluate next one.
            if (!termEval){
                continue;
            }
            
            // Term is evaluated to 1, thus XOR it to the result - where it is present.
            for(uint uctr2=0; uctr2<outputWidthUlong; uctr2++){
                oBuff[uctr2] ^= coefficients[order][outputWidthUlong*ctr + uctr2];
            }
        }
    }
    
    // Transform ULONG to output.
    this->readUlongToUchar(output, cip->getOutputBlockSize(), oBuff);
    return 0;
}

void Approximation::genMessages() {
    const unsigned byteWidth = cip->getInputBlockSize() + cip->getKeyBlockSize();
    
    // Allocate input & key buffers
    uchar * output = new uchar[cip->getOutputBlockSize()];
    uchar * finput = new uchar[byteWidth];
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
    
    delete[] finput;
    delete[] key;
    delete[] output;
}

ULONG Approximation::numberOfTerms(ULONG variables, ULONG maxOrder) const {
    ULONG res = 0;
    for(ULONG ord=0; ord<=maxOrder; ord++){
        res += CombinatiorialGenerator::binomial(variables, ord);
    }
    
    return res;
}

int Approximation::selftestIndexing() const {
    const uint bitWidth = 8*this->byteWidth;
    
    // Test quadratic indexing equations.
    cout << "Testing quadratic indexing." << endl;
    CombinatiorialGenerator cg2(bitWidth, 2);
    for(; cg2.next(); ){
        const ULONG ctr = cg2.getCounter();
        const ULONG * state = cg2.getCurState();
        const ULONG ctrComputed = CombinatiorialGenerator::getQuadIdx(bitWidth, state[0], state[1]);
        if (ctrComputed != ctr){
            dumpUlongHex(cerr, state, 2);
            cerr << "Invalid index for order 2 ctr="<<ctr<<"; computed: " << ctrComputed << endl;
        }
    }
    
    // Test cubic indexing equations.
    cout << "Testing cubic indexing." << endl << " ";
    CombinatiorialGenerator cg3(bitWidth, 3);
    ProgressMonitor pm3(0.01);
    for(; cg3.next(); ){
        const ULONG ctr = cg3.getCounter();
        const ULONG * state = cg3.getCurState();
        const ULONG ctrComputed1 = CombinatiorialGenerator::getCubeIdx(bitWidth,  state[0], state[1], state[2]);
        if (ctrComputed1 != ctr){
            dumpUlongHex(cerr, state, 3);
            cerr << "Invalid index for order 3 ctr="<<ctr<<"; computed1: " << ctrComputed1 << endl;
        }
        
        // Progress monitoring.
        double cProg = (double)ctr / (double)cg3.getTotalNum();
        pm3.setCur(cProg);
    }
    pm3.setCur(1.0);
    cout << endl;
    
    // Test general indexing.
    cout << "Testing general indexing." << endl;
    for(uint order=0; order<=orderLimit; order++){
        cout << " Testing order " << order << endl << " ";

        CombinatiorialGenerator cg(bitWidth, order);
        
        ProgressMonitor pm(0.01);
        const ULONG total = cg.getTotalNum();
        for(; cg.next(); ){
            const ULONG ctr     = cg.getCounter();            
            const ULONG * state = cg.getCurState();
            const ULONG ctrComputed = getCombinationIdx(order, state);
            if (ctrComputed != ctr){
                dumpUlongHex(cerr, state, order);
                cerr << "Invalid index for order " << order << " ctr=" << ctr << "; computed1: " << ctrComputed << endl;
            }
            
            // Progress monitoring.
            double cProg = (double)ctr / (double)total;
            pm.setCur(cProg);
        }
        pm.setCur(1.0);
        cout << endl;
    }
    
    // Randomized test of the combination indexing inversion.
    const uint testTarget = 10000;
    cout << "Testing indexing inversion" << endl << " ";
    ProgressMonitor pm(0.01);
    for(uint x=0; x<testTarget; x++){
        // Randomize order selection
        uint order2test = (rand() % orderLimit) + 1;
        ULONG comb[MAX_ORDER];
        ULONG combx[MAX_ORDER];
        
        comb[0] = rand() % (bitWidth-order2test);
        for(uint i=1; i<order2test; i++){
            const int mod = bitWidth-order2test-comb[i-1]-1;
            if (mod==0){
                comb[i] = comb[i-1]+1;
            } else {
                comb[i] = comb[i-1]+1+(rand() % mod);
            }
        }
        
        ULONG idx = getCombinationIdx(order2test, comb);
        getCombinationFromIdx(order2test, combx, idx);
        ULONG idx2 = getCombinationIdx(order2test, combx);
        
        if (idx!=idx2){
            cout << " Problem in combination order=" << order2test << "; idx=" << idx << endl;
            dumpHex(cout, comb, order2test);
            dumpHex(cout, combx, order2test);
        }
        
        // Test Ulong inversion
        idx = getCombinationULong(order2test, comb);
        getCombinationFromULong(combx, idx);
        idx2 = getCombinationULong(order2test, combx);
        
        if (idx!=idx2){
            cout << " Problem in combination order=" << order2test << "; uidx=" << idx << endl;
            dumpHex(cout, comb, order2test);
            dumpHex(cout, combx, order2test);
        }
        
        pm.setCur((double)x / (double)testTarget);
    }
    pm.setCur(1.0);
    
    cout << endl << "Test completed" << endl;
    return 0;
}

void Approximation::solveKeyGrobner(uint samples, bool dumpInputBase, bool selfTest, int basisReduction) const {
    // Allocate input & key buffers
    uchar * outputCip = new uchar[cip->getOutputBlockSize()];
    uchar * input  = new uchar[byteWidth];
    uchar * key    = new uchar[cip->getKeyBlockSize()];
    uchar * keySol = new uchar[cip->getKeyBlockSize()];
    const uint numPolynomials = numPolyActive;
    const uint numVariables = this->getNumVariables();
    
    // Seed (primitive).
    srand((unsigned)time(0)); 
    
    // Bit-mask of variables for which we have a valid value.
    ULONG * variablesValueMask = new ULONG[this->inputWidthUlong];
    memset(variablesValueMask, 0x0, sizeof(ULONG) * inputWidthUlong);
    // Set plaintext bits to 1 --> they will be evaluated.
    for(uint i=0; i<8*cip->getInputBlockSize(); i++){
        variablesValueMask[i/(8*SIZEOF_ULONG)] |= ULONG1 << (i % (8*SIZEOF_ULONG));
    }
    
    // Zero key bits are known to us, thus evaluate them with zero during partial evaluation.
    for(uint i=0; i<keybitsToZero; i++){
        const uint idx = 8*cip->getInputBlockSize() + i;
        variablesValueMask[idx/(8*SIZEOF_ULONG)] |= ULONG1 << (idx % (8*SIZEOF_ULONG));
    }
    
    // Input ULONG buffer
    ULONG * iBuff = new ULONG[this->inputWidthUlong];
    ULONG * iTmpBuff = new ULONG[this->inputWidthUlong];
    ULONG * oTmpBuff = new ULONG[this->outputWidthUlong];
    
    // FGb polynomials. We have here specific number of polynomials.
    Dpol * inputBasis1  = new Dpol[samples*numPolynomials];
    Dpol * inputBasis2 = new Dpol[samples*numPolynomials];
    Dpol * inputBasis = inputBasis1;
    Dpol * outputBasis = new Dpol[FGb_MAXI_BASE];
    memset(inputBasis1, 0, sizeof(Dpol_INT) * samples * numPolynomials);
    memset(inputBasis2, 0, sizeof(Dpol_INT) * samples * numPolynomials);
    
    // Generate key at random, once for the cipher.
    // From now it will behave as a black-box and we assume key is unknown to us.
    memset(key, 0, sizeof(uchar) * cip->getKeyBlockSize());
    for(uint i=keybitsToZero; i<8*cip->getKeyBlockSize(); i++){ 
        key[i/8] |= (rand() % 2) ? (1u << i%8) : 0; 
    }
    
    // Copy key to the input field for selftest.
    memset(input, 0, byteWidth);
    for(uint i=0; i<cip->getKeyBlockSize(); i++){
        input[i+cip->getInputBlockSize()] = key[i];
    }
    
    cout << "Generated secret key: " << endl;
    dumpHex(cout, key, cip->getKeyBlockSize());
    dumpBin(cout, key, cip->getKeyBlockSize());
    
    // Dump how polynomial-take-map looks like.
    cout << " NumVariables=" << numVariables << endl;
    cout << " NumPoly="<<numPolynomials<<"; Polymap: ";
    dumpHex(cout, poly2take, outputWidthUlong);
    cout << " Expected basis size=" << dec << (samples*numPolynomials) << endl;
    cout << " Variables to evaluate mask="; 
    dumpHex(cout, variablesValueMask, this->inputWidthUlong);
    if (selfTest){
        cout << " Self-testing of the solveGb()" << endl;
    }
    
    // Whether to show some information during sample computation or not...
    bool interSampleOutput = samples < 4;
    
    // Init FGb library.
    int step0=-1;
    int bk0=0;
    
    // Polynomial hashes
    boost::unordered_map<ULONG, uint> polynomialHashes;
    
    // Input basis size.
    uint inputBasisSize=0;
  
    // Generate tons of random messages.
    ProgressMonitor pmSample(0.01);
    for(ulong sample=0; sample<samples; sample++){
        if (interSampleOutput){
            cout << endl << " [+] Starting with sample="<<sample<<endl;
        }
        
        // Generate message at random.
        // Fix plaintext variables to the generated ones. 
        // Now we obtain system of equations with key variables.
        // 128 equations with (128-keybitsToZero) unknown bits.
        for(uint i=0; samples > 1 && i<cip->getInputBlockSize(); i++){ 
            input[i] = (rand() % 0x100); 
        }
        
        // Input variables, only masked are taken into consideration during partial evaluation.
        this->readUcharToUlong(input, cip->getInputBlockSize(), iBuff);
        this->readUcharToUlong(key,   cip->getKeyBlockSize(),   iBuff + OWN_CEIL((double)cip->getInputBlockSize() / (double)SIZEOF_ULONG));
        
        // Evaluation.
        if (selfTest){
            // If we are using self test, evaluate this on the approximation function.
            this->evaluateCoefficients(input, outputCip, iTmpBuff, oTmpBuff);
        } else {
            // Evaluate cipher.
            cip->evaluate(input, key, outputCip);       
        }
        
        // New function is stored in another coefficient array.
        // Allocate space for the coefficients.
        std::vector<ULONG> coeffEval[MAX_ORDER];
        for(unsigned int order = 0; order<=orderLimit; order++){
            ULONG vecSize = CombinatiorialGenerator::binomial(numVariables, order) * outputWidthUlong;
            if (interSampleOutput){
                cout << " Allocating pEval function, order=" << dec << order << "; vecSize=" << vecSize << endl;
            }
            
            coeffEval[order].assign(vecSize, (ULONG) 0);
        }
        
        // Partial evaluation = reduces terms with evaluated variables.
        if (interSampleOutput){
            cout << " Going to partially evaluate approximating function." << endl;
        }
        partialEvaluation(numVariables, variablesValueMask, iBuff, coeffEval);
        
        // Add ciphertext values to the polynomials to obtain system of equations.
        // For particular input message.
        for(uint x=0; x<cip->getOutputBlockSize(); x++){
            coeffEval[0][x/SIZEOF_ULONG] ^= ((ULONG)outputCip[x])<<(8 * (x % SIZEOF_ULONG));
        }
        
        // Proceed polynomial per polynomial, build input basis for FGb.
        if (interSampleOutput){
            cout << " Generating input basis." << endl << " ";
        }
        
        ProgressMonitor pmBasis(0.01);
        ULONG numTermsSum=0;
        uint polyCtr=0;
        for(uint poly=0; poly<numPolynomials; poly++){
            // If this polynomial is not selected, do not add it in the input base.
            if (isPoly2Take(poly)==false){
                continue;
            }
            
            // Convert out internal polynomial representation to FGb representation.
            ULONG numTerms = 0;
            ULONG hash = 0;
            uint curPolyIdx = sample*numPolynomials + polyCtr;
            inputBasis1[curPolyIdx] = polynomial2FGb(numVariables, coeffEval, orderLimit, poly, &numTerms, &hash);
            numTermsSum += numTerms;
            
            // If zero terms, it is null, do not add it to the base since 
            // it does not increase basis dimension.
            if (numTerms==0){
                continue;
            }
            
            // Polynomial duplicity check, do not add duplicates to the input base.
            // Since it does not increase basis dimension.
            if (polynomialHashes.count(hash)>0){
                //cout << "Polynomial with idx=" << (curPolyIdx) << " is already present in the basis, idx=" << polynomialHashes[hash] << "; hash=" << hex << hash << endl;
                continue;
            } else {
                polynomialHashes.insert(std::pair<ULONG,uint>(hash, curPolyIdx));
            }
            
            // Polynomial is unique. Add to the basis.
            if (interSampleOutput){
                //cout << setw(4) << right << curPolyIdx << " is f_" << poly << endl;
                pmBasis.setCur((double)poly / double(numPolynomials));
            }
            
            polyCtr+=1;
        }
        
        // Update complete basis size.
        inputBasisSize+=polyCtr;
        
        // If verbose mode, finish progress bar & write average number of terms in polynomials.
        if (interSampleOutput){
            pmBasis.setCur(1.0);
            cout << "; Number of terms on average: " << ((double)numTermsSum / (double)numPolynomials) << endl;
        } else {
            // If there is no intersample output, show a general progress bar.
            pmSample.setCur((double)sample / (double)samples);
        }
    }
    
    // Hack: we don't want over determined system for now, take last
    // numVariables equations to the input base.
    if (basisReduction > 0 && inputBasisSize > numVariables){
        for(uint i=0; i<numVariables; i++){
            inputBasis2[i] = inputBasis1[inputBasisSize-1-i];
        }
        
        cout << " Number of equations reduced to " << numVariables << endl;
        if (dumpInputBase){
            cout << " Original basis, size=" << inputBasisSize << endl;
            dumpBasis(numVariables, inputBasis1, inputBasisSize);
        } 
        
        inputBasisSize = numVariables;
        inputBasis = inputBasis2;
    }
    
    // Only if general progressbar is shown.
    if (!interSampleOutput){
        pmSample.setCur(1.0);
        cout << endl;
    }
    
    // Print out input basis
    if (dumpInputBase){
        cout << " Input basis, size=" << inputBasisSize << endl;
        dumpBasis(numVariables, inputBasis, inputBasisSize);
    }
    
    // Compute Gb.
    {
    int nb;
    double t0;
    const int n_input=inputBasisSize; /* we have X polynomials on input */
    struct sFGB_Comp_Desc Env;
    FGB_Comp_Desc env=&Env;
    env->_compute=FGB_COMPUTE_GBASIS; /* The following function can be used to compute Gb, NormalForms, RR, ... 
                                         Here we want to compute a Groebner Basis */
    env->_nb=0; /* parameter is used when computing NormalForms (see an example in bug_prog2.c */
    env->_force_elim=0; /* if force_elim=1 then return only the result of the elimination 
                           (need to define a monomial ordering DRL(k1,k2) with k2>0 ) */
    env->_off=0;       /* should be 0 for modulo p computation	*/

    env->_index=900000; /* This is is the maximal size of the matrices generated by F4 
                          you can increase this value according to your memory */
    env->_zone=0;    /* should be 0 */
    env->_memory=0;  /* should be 0 */
    /* Other parameters :
       t0 is the CPU time (reference to a double)
       bk0 : should be 0 
       step0 : this is the number primes for the first step
               if step0<0 then this parameter is automatically set by the library
     */
    cout << " [+] Going to compute GB, n_input="<<dec<<n_input<<endl;
    nb=FGB(groebner)(inputBasis,n_input,outputBasis,1,0,&t0,bk0,step0,0,env);

    // For now just print out the Grobner basis.
    dumpBasis(numVariables, outputBasis, nb);
    cout << "Input basis dimension=" << n_input << endl;
    cout << "Basis dimension=" << nb << endl;
    cout << "Number of variables=" << numVariables << endl;
    cout << "Stats: cpu=" << t0 << "; memory=" << env->_memory << "; zone=" << env->_zone << endl;
    
    // TODO: use linearization trick.
    // TODO: use http://icm.mcs.kent.edu/reports/1995/gb.pdf to solve the system.
    
    cout << " [+] Going to solve GB" << endl;
    int haveSol = solveGb(numVariables, outputBasis, nb, keySol);
    if (haveSol >= 0){
        double hits = 0.0;
        double hitRatio;
        for(uint i=0; i<numVariables; i++){
            // If number of variables is reduced, take this into account
            // for precise key bit hit ratio (real key is shifted).
            const uint keyIdx = 8*cip->getKeyBlockSize() - numVariables + i;
            hits += ((key[keyIdx/8] & (1u << (keyIdx%8))) > 0) == ((keySol[i/8] & (1u << (i%8))) > 0);
        }
        
        hitRatio = hits / (double)numVariables;
        cout << "Key bit-hit ratio: " << hitRatio << endl;
    }
    }
    
    delete[] key;
    delete[] keySol;
    delete[] input;
    delete[] outputCip;
    delete[] variablesValueMask;
    delete[] iBuff;
    delete[] iTmpBuff;
    delete[] oTmpBuff;
    delete[] inputBasis1;
    delete[] inputBasis2;
    delete[] outputBasis;
}

int Approximation::solveGb(uint numVariables, Dpol* basis, uint numPoly, uchar * solvedKey) const {
    // Detect if system has no solution.
    // This happens if Gb = <1>.
    if (numPoly==1){
        const I32 nb_mons = FGB(nb_terms)(basis[0]);
        if (nb_mons==1){
            cout << "   System has only 1 variable, 1 term. Probably Gb=<1> what means there is no solution. NumVariables=" << numVariables << endl;
            return -2;
        }
        
        cout << "   System has only 1 equation. Cannot solve such under-determined system. NumVariables=" << numVariables << endl;
        return -3;
    }
    
    // If system is too under-determined, do not solve it. 
    if (numVariables > (numPoly+4)){
        cout << "   Cannot solve (under-determined) system with " << numVariables << " and " << numPoly << " equations right now right now" << endl;
        return -1;
    }
    
    // TODO: Use linearization trick, now do not use it.
    // In order to solve the system, we need to separate constant terms from the
    // system of equation to a side vector b
    vec_GF2 b(INIT_SIZE, numPoly);
    
    // TODO: in linearization trick, it is needed to substitute complex terms 
    // with a new linear variable. We have to count the nonspecific terms
    // and assign an unique linear representation to each of them.
    // Thus size of the matrix is known only AFTER this process of linearization.
    mat_GF2 systm(INIT_SIZE, numPoly, numVariables);
    
    // Add polynomial to the systm matrix, constant term to the b vector.
    for(uint polyIdx = 0; polyIdx < numPoly; polyIdx++){
        const I32 nb_mons = FGB(nb_terms)(basis[polyIdx]);  // Number of Monomials.
        I32* Mons = new I32[numVariables * nb_mons];        // Exponents for variables in terms.
        I32* Cfs = new I32[nb_mons];                        // Coefficients for terms.
        FGB(export_poly)(numVariables, nb_mons, Mons, Cfs, basis[polyIdx]);
        I32 j;
        for (j = 0; j < nb_mons; j++) {
            UI32 k, is_one = 1;
            I32* ei = Mons + j*numVariables;

            // In GF(2) all non-zero coefficients are 1.
            // But check for null anyway.
            if (Cfs[j]==0) {
                continue;
            }
            
            for (k = 0; k < numVariables; k++){
                // Exponent of the variable k in the term j.
                if (ei[k]==0) continue;
                
                // Coping with the exponents that should not be here - XOR
                if (IsZero(systm.get(polyIdx, k))){
                    systm.put(polyIdx, k, 1ul);
                } else {
                    systm.put(polyIdx, k, 0ul);
                }
                
                is_one = 0;
            }
            
            // Constant term, should be only 1 in the polynomial.
            if (is_one) {
                b.put(polyIdx, 1ul);
            }
        }

        delete[] Mons;
        delete[] Cfs;
    }
    
    //dumpVector(b);
    //dumpMatrix(systm);
    
    // Try to solve the system
    GF2 determinant;
    vec_GF2 solution;
    
    solve(determinant, solution, systm, b);
    if (IsZero(determinant)){
        cout << "   Determinant is zero, cannot solve this sytem." << endl;
        return -4;
    }
    
    // We have solution, convert it to the uchar array and dump in hexa and in the binary.
    const uint solByteSize = (uint)ceil(numVariables/8.0);
    memset(solvedKey, 0x0, solByteSize);
    
    for(uint i=0; i<numVariables; i++){
        if (IsOne(solution.get(i))){
            solvedKey[i/8] |= 1u << (i%8);
        }
    }
    
    cout << "   We have the solution:" << endl;
    dumpHex(cout, solvedKey, solByteSize);
    dumpBin(cout, solvedKey, solByteSize);
    return 1;
}

void Approximation::dumpFGbPoly(uint numVariables, Dpol poly) const {
    // Import the internal representation of each polynomial computed by FGb.
    const I32 nb_mons = FGB(nb_terms)(poly);        // Number of Monomials.
    I32* Mons = new I32[numVariables * nb_mons];    // (UI32*) (malloc(sizeof (UI32) * numVariables * nb_mons));
    I32* Cfs = new I32[nb_mons];                    // (I32*) (malloc(sizeof (I32) * nb_mons));
    FGB(export_poly)(numVariables, nb_mons, Mons, Cfs, poly);
    I32 j;
    for (j = 0; j < nb_mons; j++) {

        UI32 k, is_one = 1;
        I32* ei = Mons + j*numVariables;

        if (j > 0) {
            cout << "+";
        }

        cout << Cfs[j];
        for (k = 0; k < numVariables; k++)
            if (ei[k]) {
                if (ei[k] == 1) {
                    cout << varNames[k];
                } else {
                    cout << varNames[k] << "^" << ei[k];
                }
                is_one = 0;
            }
        if (is_one) {
            cout << "*1";
        }
    }

    delete[] Mons;
    delete[] Cfs;
}

void Approximation::dumpBasis(uint numVariables, Dpol* basis, uint numPoly) const {
    cout << "[ len=" << numPoly << endl;
    for (uint i = 0; i < numPoly; i++) {
        // Use this fuction to print the result.
        // FGB(see_Dpol)(outputBasis[i]);
        
        cout << "# " << dec << setw(4) << setfill('0') << right << i << ": ";
        dumpFGbPoly(numVariables, basis[i]);
        if (i < (numPoly - 1)) {
            cout << endl;
        }
    }
    cout << "]" << endl;
}

void Approximation::initFGb(uint numVariables) const {
    FGB(enter)(); /* First thing to do : GMP original memory allocators are saved */
    
    // Do not change the following parameters (change will cause runtime error):
    //   2 is the number of bytes of each coefficients so the
    //     maximal prime is < 2^16
    //   2 is the number of bytes of each exponent:
    //     it means that each exponent should be < 2^15
    FGB(init_urgent)(2,MAPLE_FGB_BIGNNI,"DRLDRL",FGb_MAXI_BASE,0);
    
    FGB(init)(1,1,0,fgbFile); /* do not change */
    {
      UI32 pr[]={(UI32)(2)}; /* We compute in GF(2)[x1,x2,...,x_{8*byteWidth}] */
      FGB(reset_coeffs)(1,pr);
    }
    
    {
      FGB(reset_expos)(numVariables,0,varNames);  /* Define the monomial ordering: DRL(k1,k2) where 
                                    k1 is the size of the 1st block of variables 
                                    k2 is the size of the 2nd block of variables 
                                    and k1+k2=nb_vars is the total number of variables
                                   */
    }
}

void Approximation::deinitFGb() const {
    resetFGb();
    FGB(exit)(); /* restore original GMP allocators */
}

void Approximation::resetFGb() const {
    FGB(reset_memory)(); /* to reset Memory */
}


int Approximation::partialEvaluation(uint numVariables, ULONG * variablesValueMask, ULONG * iBuff, std::vector<ULONG> * coeffEval) const{
    assert(numVariables <= 8*byteWidth);
    ULONG * newTerm = new ULONG[orderLimit+1];
    for(uint order = 0; order <= orderLimit; order++){
        CombinatiorialGenerator oldCg(8*byteWidth, order);
        
        // Iterate over coefficients for the higher function and update coefficients for
        // the lower function terms.
        for(; oldCg.next(); ){
            const ULONG ctr = oldCg.getCounter();
            // Now term being evaluated is fixed, defined by the state
            // of the combinatorial generator. 
            //
            // Get bit-mask with those bits enabled corresponding to variables in
            // the particular term determined by oldCg.
            const ULONG * comb  = oldCg.getCurUlongCombination();
            
            // Evaluate particular term on the input.
            // Here we evaluate only specified variables on the input. Result
            // of this evaluation becomes a new term coefficient for the low function term.
            bool termEval=true;
            for(uint uctr2=0; uctr2<inputWidthUlong; uctr2++){
                // Term to evaluate will be obtained after masking with provided mask.
                const ULONG term = comb[uctr2] & variablesValueMask[uctr2];
                // Evaluate masked term on the masked input (only for masked values 
                // the input has sense).
                termEval &= (term & iBuff[uctr2] & variablesValueMask[uctr2]) == term;
            }
            
            // If term is null thus cannot be in the low function under this input.
            if (!termEval){
                continue;
            }
            
            // This term is present in the low function.
            // Determine the low term index manually.
            uint torder = 0;
            uint curRemoved = 0; // number of removed variables in interval 0..current.
            memset(newTerm, 0, orderLimit * SIZEOF_ULONG);
            for(uint uctr2=0; uctr2<inputWidthUlong; uctr2++){
                // Examine each bit separately.
                // If we already have term of the highest possible order, stop it.
                for(uint x=0; x<SIZEOF_ULONG*8 && torder<=orderLimit; x++){
                    // Optimization - skipping large consecutive blocks.
                    if ((x & 0x7) == 0){
                        // If there is a full skip (all variables evaluated, skip.)
                        if (((variablesValueMask[uctr2]>>x) & 0xff) == 0xff){
                            curRemoved+=8;
                            x+=7;
                            continue;
                        }
                    }
                    
                    // Future optimization: compute hamming weight to an external table.
                    // If current bit is set to 1 in the mask -> was evaluated already (not a variable anymore).
                    if (variablesValueMask[uctr2] & (ULONG1 << x)){
                        curRemoved+=1;
                        continue;
                    }
                    
                    // Mask is zero here -> it is a variable. If it is present in 
                    // the termRest, we have the winner.
                    if ((comb[uctr2] & (ULONG1 << x)) > 0){
                        // This bit is set, add to the combination.
                        newTerm[torder++] = (uctr2*SIZEOF_ULONG*8)+x-curRemoved;
                    }
                }
            }
            assert(torder <= orderLimit);
            assert(curRemoved <= 8*byteWidth);
            
            // Compute index of this term and toggle it in the low function.
            // Number of variables is reduced, thus set the offset on the number combinations.
            ULONG newTermIdx = getCombinationIdx(torder, newTerm, 0, 8*byteWidth-numVariables);
            
            // Debugging code:
            //cout << " termOrder="<<torder<<"; from term order="<<order<<";  original term=";
            //dumpHex(cout, oldCg.getCurState(), order);
            //cout << "  new=";
            //dumpHex(cout, newTerm, torder);
            
            // After partial evaluation term has coefficient 1 and can be present
            // in polynomials under this input.
            // Thus if the higher term (originating one) was present in the particular
            // polynomial, add this low term (reduced/evaluated) to the corresponding polynomial.
            for(uint uctr2=0; uctr2<outputWidthUlong; uctr2++){
                coeffEval[torder][outputWidthUlong*newTermIdx + uctr2] ^= coefficients[order][outputWidthUlong*ctr + uctr2];
            }
        }
    }
    
    delete[] newTerm;
    return 0;
}

Dpol_INT Approximation::polynomial2FGb(uint numVariables, std::vector<ULONG>* coefs, uint maxOrder, uint polyIdx, ULONG * numTerms, ULONG * hash) const {
    ULONG termsEnabled = 0;
    Dpol_INT prev;
    I32 * termRepresentation = new I32[numVariables];
    MD5_CTX md5Ctx;
    
    // Has to determine exact number of enabled terms in the polynomial.
    const uint coefSegment = polyIdx / (SIZEOF_ULONG * 8);
    const uint coefOffset = polyIdx % (SIZEOF_ULONG * 8);
    for(uint order = 0; order <= orderLimit; order++){
        CombinatiorialGenerator cg(numVariables, order);
        for(; cg.next(); ){
            const ULONG ctr = cg.getCounter(); 
            termsEnabled += (coefs[order][outputWidthUlong*ctr+coefSegment] & (ULONG1<<(coefOffset))) > 0;
        }
    }
    
    if (hash!=NULL){
        MD5_Init(&md5Ctx);
        MD5_Update(&md5Ctx, &termsEnabled, SIZEOF_ULONG);
    }
    
    // Create an empty polynomial with specified number of terms. 
    prev=FGB(creat_poly)(termsEnabled);
    
    // Iterate again - now construct the polynomial.
    uint termCounter=0;
    for(uint order = 0; order <= orderLimit && termCounter < termsEnabled; order++){
        CombinatiorialGenerator cg(numVariables, order);
        for(; cg.next() && termCounter < termsEnabled; ){
            const ULONG ctr = cg.getCounter();
            if ((coefs[order][outputWidthUlong*ctr+coefSegment] & (ULONG1<<(coefOffset))) == 0) continue;
            
            // Update hash value if we want to compute it. 
            if (hash!=NULL){
                MD5_Update(&md5Ctx, &ctr, SIZEOF_ULONG);
            }
            
            // Coefficient is present, set term variables to representation.
            const ULONG * state = cg.getCurState();
            memset(termRepresentation, 0, sizeof(I32) * numVariables);
            for(uint x=0; x<order; x++){
                termRepresentation[state[x]]=1;
            }
            
            // Set term to the polynomial.
            FGB(set_expos2)(prev, termCounter, termRepresentation, numVariables);
            
            // Set term coefficient (simple, always 1 since we are in GF(2)).
            FGB(set_coeff_I32)(prev, termCounter, 1);
            
            termCounter+=1;
        }
    }
    
    // Final hash computation
    if (hash!=NULL){
        ULONG hash1, hash2;
        uchar md[MD5_DIGEST_LENGTH];
        
        MD5_Final(md, &md5Ctx);
        readUcharToUlong(md,              SIZEOF_ULONG, &hash1);
        readUcharToUlong(md+SIZEOF_ULONG, SIZEOF_ULONG, &hash2);
        
        *hash = hash1 ^ hash2;
    }
    
    // Sorting for GB, has to be done and is very slowly...
    FGB(full_sort_poly2)(prev);
    
    if (numTerms!=NULL){
        *numTerms = termsEnabled;
    }
    
    delete[] termRepresentation;
    return prev;
}

void Approximation::readUcharToUlong(const uchar * input, uint size, ULONG * iBuff) const {
    // At first reset memory of the big buffer.
    memset(iBuff, 0, size);
    // And read particular parts of the input to the big buffer.
    for(uint x = 0; x < size; x++){
        iBuff[x/SIZEOF_ULONG] = READ_TERM_1(iBuff[x/SIZEOF_ULONG], input[x], x%SIZEOF_ULONG);
    }
}

void Approximation::readUlongToUchar(uchar* output, uint size, const ULONG* iBuff) const {
    for(uint x=0; x<size; x++){
        output[x] = (iBuff[x/SIZEOF_ULONG] >> (8* (x % SIZEOF_ULONG))) & ((unsigned char)0xffu);
    }
}

//
// FGb helper
//
extern "C" {
    //FILE* log_output;
    void info_Maple(const char* s)
    {
      /* 
         if (verbose)
         {
         fprintf(stderr,"%s",s);
         fflush(stderr);
         }
      */
    }

    void FGb_int_error_Maple(const char* s)
    {
      fprintf(stderr,"Error: %s",s);
      fflush(stderr);
      exit(3);
    }

    void FGb_error_Maple(const char* s)
    {
      FGb_int_error_Maple(s);
    }

    void FGb_checkInterrupt()
    {
        
    }

    void FGb_int_checkInterrupt()
    {
    }

    void FGb_push_gmp_alloc_fnct(void *(*alloc_func) (size_t),
                                 void *(*realloc_func) (void *, size_t, size_t),
                                 void (*free_func) (void *, size_t))
    {
    }

    void FGb_pop_gmp_alloc_fnct()
    {
    }
}
