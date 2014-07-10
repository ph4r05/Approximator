/* 
 * File:   Approximation.cpp
 * Author: ph4r05
 * 
 * Created on May 16, 2014, 10:21 AM
 * TODO: multithreaded!!!
 * TODO: MPI???
 */
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
#include <boost/unordered_map.hpp>
#include <thread>

#include "Approximation.h"
#include "CombinatiorialGenerator.h"
#include "ProgressMonitor.h"
#include "NTLUtils.h"
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

// Signal blocking
#include <sys/types.h>
#include <signal.h>
#include <stdlib.h>

NTL_CLIENT
using namespace std;
using namespace NTL;

Approximation::Approximation(uint orderLimit) : cip(NULL), dumpCoefsToFile(false), 
        outputWidthUlong(0), inputWidthUlong(0),
        threadCount(1), keybitsToZero(0),
        poly2take(NULL), numPolyActive(0),
        verboseLvl(1) {
    
    this->orderLimit = orderLimit;
}

Approximation::~Approximation() {       
    if (poly2take!=NULL){
        delete[] poly2take;
        poly2take = NULL;
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
    
    // Init combinatorial indexer.
    combIndexer.init(8*byteWidth, orderLimit);
    
    // Init FGb helper
    fgb.init(byteWidth, orderLimit, cip->getOutputBlockSize()*8);
    
    // Init polynomial-to-take bitmap.
    poly2take = new ULONG[outputWidthUlong];
    numPolyActive = 8*cip->getOutputBlockSize();
    memset(poly2take, 0xff, SIZEOF_ULONG * outputWidthUlong);
    
    // Initialize signal blocking.
    sigemptyset(&pendingSignals);
    sigemptyset(&blockingMask);
    sigaddset(&blockingMask, SIGHUP);
    sigaddset(&blockingMask, SIGTERM);
    sigaddset(&blockingMask, SIGQUIT);
    sigaddset(&blockingMask, SIGINT);   // CTRL+C
    sigaddset(&blockingMask, SIGABRT);
    sigaddset(&blockingMask, SIGUSR1);
    sigaddset(&blockingMask, SIGUSR2);
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

void Approximation::computeCoefficients(std::vector<ULONG> * coefficients) {
    ULONG * ulongOut = new ULONG[outputWidthUlong];
    ULONG * ulongInp = new ULONG[inputWidthUlong];
    uchar * finput   = new uchar[byteWidth];
    uchar * output   = new uchar[cip->getOutputBlockSize()];
    
    CombinatiorialGenerator ** cgenerators = new CombinatiorialGenerator * [orderLimit+1];
    ULONG * tmpCombination = new ULONG[orderLimit+1];
    
    // Allocating space for the coefficients.
    for(unsigned int order = 0; order<=orderLimit; order++){
        ULONG vecSize = CombinatiorialGenerator::binomial(8*byteWidth, order) * outputWidthUlong;
        cout << "  Allocating coefficient storage for order " << order << "; Bytes=" << vecSize * sizeof(ULONG) << endl;
        coefficients[order].assign(vecSize, (ULONG) 0);
    }
    
    // TODO: for order 3 and higher use multiple threads to parallelize the computation!
    // TODO: regarding the progress bar, only thread num 0 shows it. Does not mind
    // since each thread has space to search of the same size, progress should be
    // approximately the same for each thread.
    
    // Find polynomial terms coefficient for order 1..orderLimit.
    for(uint order=0; order<=orderLimit; order++){
        CombinatiorialGenerator cg(byteWidth*8, order);
        
        // Combinatorial generators for computing XOR indices.
        for(uint i=0; i < order; i++){
            cgenerators[i] = new CombinatiorialGenerator(order, i);
        }
        
        cout << "Starting with order: " 
                << order 
                << "; combinations: " << cg.getTotalNum()
                << "; number of bytes to store coefficients: " 
                << (8 * cip->getOutputBlockSize() * cg.getTotalNum() / 8)
                << endl << " ";
        
        // Here is the point for parallelization.
        // Each thread/computing node can compute x-th combination from the generator
        // effectively partitioning combination space. 
        // Synchronization barrier is needed after finishing particular order
        // because in order to compute order N we need to have coefficients of
        // terms of order N-1 and less.
        ULONG combCtr=0;
        ProgressMonitor pm(0.01);
        for(; cg.next(); combCtr++){
            const uchar * input = cg.getCurCombination();
            
            // Evaluate cipher on current combination.
            cip->evaluate(input, input + cip->getInputBlockSize(), output);
            
            // Transform output to the ULONG array
            // For better memory handling and XORing in an one big register.
            readUcharToUlong(output, cip->getOutputBlockSize(), ulongOut);
            
            // Evaluate coefficients for each polynomial in the cipher
            // for the given term specified by the state of the combinatorial generator. 
            //
            // Current term value: all previous terms including enabled bits XOR ciphertext.
            // For example, term to determine: x1x6x9:
            //   constant                XOR
            //   x1   XOR x6   XOR x9    XOR
            //   x1x6 XOR x1x9 XOR x6x9
            //
            // XOR all constant, linear, quadratic, cubic, etc.. terms if applicable.
            // In order to determine current term coefficient all lower terms
            // that can be obtained from this one has to be taken into account.
            // and XORed into the result.
            for(uint xorOrder=0; xorOrder<order && order<=orderLimit; xorOrder++){
                // Obtain a shortcut reference to the combinatorial generator
                // for this xorOrder. This generator is used to generate 
                // combinations of variables from to original term to construct
                // a lower term.
                CombinatiorialGenerator * const xorOrderGen = cgenerators[xorOrder];

                // Iterate over all possible variable combinations to the new 
                // resulting term of a lower order.
                for(xorOrderGen->reset(); xorOrderGen->next(); ){
                    // Construct xorOrder term representation for term
                    // index computation.
                    for(uint tmpCombCtr = 0; tmpCombCtr<xorOrder; tmpCombCtr++){
                        tmpCombination[tmpCombCtr] = cg.getCurState()[xorOrderGen->getCurState()[tmpCombCtr]];
                    }

                    // Term index computation.
                    ULONG idx = combIndexer.getCombinationIdx(xorOrder, tmpCombination);

                    // XOR current result register with coefficient register.
                    // If low order term is present in the approximation function,
                    // it has to be taken into account.
                    for(uint ulongCtr=0; ulongCtr<outputWidthUlong; ulongCtr++){
                        ulongOut[ulongCtr] ^= coefficients[xorOrder][outputWidthUlong*idx + ulongCtr];
                    }
                }
            }

            // Store value of the current coefficient on his place in the coef. vector.
            for(uint ulongCtr=0; ulongCtr<outputWidthUlong; ulongCtr++){
                coefficients[order][outputWidthUlong*cg.getCounter() + ulongCtr] = ulongOut[ulongCtr];
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
        
        readUcharToUlong(input, byteWidth, iBuff);
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(this->coefficients, input, outputPol, ulongInp, ulongOut);
        
        // Compute statistics - number of hits for individual polynomial.
        for(uint p=0; p<8*cip->getOutputBlockSize(); p++){
            hits[p] += (outputCip[p/8] & (ULONG1 << (p%8))) == (outputPol[p/8] & (ULONG1 << (p%8)));
        }
        
        // Test partial evaluation for correctness, has to be in match with full evaluation.
        std::vector<ULONG> coeffEval[MAX_ORDER];
        coeffEval[0].assign(this->inputWidthUlong, 0ul);
        this->partialEvaluation(this->coefficients, this->byteWidth*8, variablesValueMask, iBuff, coeffEval);
        
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
        randomBuffer(input, byteWidth);
        
        // Evaluate cipher.
        cip->evaluate(input, input + cip->getInputBlockSize(), outputCip);
        
        // Evaluate polynomial.
        this->evaluateCoefficients(this->coefficients, input, outputPol, ulongInp, ulongOut);
        
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

int Approximation::evaluateCoefficients(const std::vector<ULONG> * coefficients,
        const unsigned char* input, unsigned char* output, ULONG * iBuff, ULONG * oBuff) const {
    // We can assume that approximate half of the coefficients are enabled/present
    // in the resulting polynomial, thus evaluation is based on the iteration of 
    // the combinatorial generator and reading coefficient by coefficient.
    const uint bitWidth = 8*byteWidth;
    
    // Reset output buffer, only ones will be set here, has to be set to zero
    // and copy to bigger buffer for better manipulation & speed.
    readUcharToUlong(input, byteWidth, iBuff);
    
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
    readUlongToUchar(output, cip->getOutputBlockSize(), oBuff);
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
    randomBuffer(key, cip->getKeyBlockSize());
    dumpUchar(cip3, key, cip->getKeyBlockSize());
    
    // Generate tons of random messages.
    for(unsigned long i=0; i<16777216ul; i++){
        // Generate message at random.
        randomBuffer(finput, cip->getInputBlockSize());
        
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
        ULONG outerCtr = 0;
        for(; cg.next(); outerCtr++){
            const ULONG ctr     = cg.getCounter();            
            const ULONG * state = cg.getCurState();
            const ULONG ctrComputed = combIndexer.getCombinationIdx(order, state);
            if (ctrComputed != ctr){
                dumpUlongHex(cerr, state, order);
                cerr << "Invalid index for order " << order << " ctr=" << ctr << "; computed1: " << ctrComputed << endl;
            }
            
            // Progress monitoring.
            double cProg = (double)ctr / (double)total;
            pm.setCur(cProg);
        }
        if (outerCtr != total){
            cerr << "Invalid number of iterations for configuration (" << bitWidth << ", " << order << ")!" << endl;
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
        
        ULONG idx = combIndexer.getCombinationIdx(order2test, comb);
        combIndexer.getCombinationFromIdx(order2test, combx, idx);
        ULONG idx2 = combIndexer.getCombinationIdx(order2test, combx);
        
        if (idx!=idx2){
            cout << " Problem in combination order=" << order2test << "; idx=" << idx << endl;
            dumpHex(cout, comb, order2test);
            dumpHex(cout, combx, order2test);
        }
        
        // Test Ulong inversion
        idx = combIndexer.getCombinationULong(order2test, comb);
        combIndexer.getCombinationFromULong(combx, idx);
        idx2 = combIndexer.getCombinationULong(order2test, combx);
        
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
    cout << " NumVariables=" << numVariables << "; zeroKeyBits=" << keybitsToZero << endl;
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
        readUcharToUlong(input, cip->getInputBlockSize(), iBuff);
        readUcharToUlong(key,   cip->getKeyBlockSize(),   iBuff + OWN_CEIL((double)cip->getInputBlockSize() / (double)SIZEOF_ULONG));
        
        // Evaluation.
        if (selfTest){
            // If we are using self test, evaluate this on the approximation function.
            this->evaluateCoefficients(this->coefficients, input, outputCip, iTmpBuff, oTmpBuff);
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
        partialEvaluation(this->coefficients, numVariables, variablesValueMask, iBuff, coeffEval);
        
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
            inputBasis1[curPolyIdx] = fgb.polynomial2FGb(numVariables, coeffEval, orderLimit, poly, &numTerms, &hash);
            numTermsSum += numTerms;
            
            // If zero terms, it is null, do not add it to the base since 
            // it does not increase basis dimension.
            if (numTerms==0 || fgb.isPoly1(inputBasis1[curPolyIdx], numVariables)){
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
    
    // Hack: we don't want over-determined system for now, take last
    // numVariables equations to the input base.
    if (basisReduction > 0 && inputBasisSize > numVariables){
        for(uint i=0; i<numVariables; i++){
            inputBasis2[i] = inputBasis1[inputBasisSize-1-i];
        }
        
        cout << " Number of equations reduced to " << numVariables << endl;
        if (dumpInputBase){
            cout << " Original basis, size=" << inputBasisSize << endl;
            fgb.dumpBasis(numVariables, inputBasis1, inputBasisSize);
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
        fgb.dumpBasis(numVariables, inputBasis, inputBasisSize);
    }
    
    // Compute Gb.
    cout << " [+] Going to compute GB, n_input="<<dec<<inputBasisSize<<endl;
    double t0 = 0.0;
    int nb = fgb.computeFGb(inputBasisSize, inputBasis, outputBasis, &t0);
    
    // For now just print out the Grobner basis.
    fgb.dumpBasis(numVariables, outputBasis, nb);
    cout << "Input basis dimension=" << inputBasisSize << endl;
    cout << "Basis dimension=" << nb << endl;
    cout << "Number of variables=" << numVariables << endl;
    cout << "Stats: cpu=" << t0 << endl;
    
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
        if (fgb.isPoly1(basis[0], numVariables)){
            cout << "   System generates ideal=<1> what means there is no solution. NumVariables=" << numVariables << endl;
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
        const I32 nb_mons = fgb.getNumberOfTerms(basis[polyIdx]);  // Number of Monomials.
        I32* Mons = new I32[numVariables * nb_mons];        // Exponents for variables in terms.
        I32* Cfs = new I32[nb_mons];                        // Coefficients for terms.
        I32 j;
        
        fgb.exportPolynomial(numVariables, nb_mons, Mons, Cfs, basis[polyIdx]);
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
    if (numVariables > systm.NumRows()){
        cout << "   System is under-determined, rows=" << (systm.NumRows()) << endl;
        return -8;
    }
    
    // Try to solve the system
    mat_GF2 gaussed;
    long rank = gaussPh4r05(gaussed, systm, b, numVariables);
    if (rank < numVariables){
        cout << "   Determinant is zero, cannot solve this system. Rank=" << rank << endl;
        return -4;
    }
    
    // We have solution, convert it to the uchar array and dump in hexa and in the binary.
    const uint solByteSize = (uint)ceil(numVariables/8.0);
    memset(solvedKey, 0x0, solByteSize);
    
    for(uint i=0; i<numVariables; i++){
        if (IsOne(gaussed.get(i, numVariables))){
            solvedKey[i/8] |= 1u << (i%8);
        }
    }
    
    cout << "   We have the solution:" << endl;
    dumpHex(cout, solvedKey, solByteSize);
    dumpBin(cout, solvedKey, solByteSize);
    return 1;
}

int Approximation::partialEvaluation(const std::vector<ULONG> * coefficients,
        uint numVariables, ULONG * variablesValueMask, ULONG * iBuff, std::vector<ULONG> * coeffEval) const{
    
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
            ULONG newTermIdx = combIndexer.getCombinationIdx(torder, newTerm, 0, 8*byteWidth-numVariables);
            
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

uint Approximation::genBitPosMask(uint* termBitPositions, uint bitWidth, const ULONG* termMask, uint termWeight) const {
    for(uint pos=0, bpos=0; pos<bitWidth ; pos++){
        const ULONG bmask = (ULONG1 << (pos % (8*SIZEOF_ULONG)));
        if ((termMask[(pos / (8*SIZEOF_ULONG))] & bmask) != bmask) continue;
        
        assert(bpos < termWeight);
        termBitPositions[bpos++] = pos;
    }
    
    return 1;
}

int Approximation::subCubeTerm(uint termWeight, const ULONG* termMask, const uchar* finput, 
        ULONG* subcube, uint step, uint offset, uint subCubes, bool precompKey) const {
    const uint bitWidth = 8*byteWidth;
    const uint outByteWidth = cip->getOutputBlockSize();
    
    // Function input/output for evaluation.
    uchar output[cip->getOutputBlockSize()];
    uchar input[byteWidth];
    
    // Local ULONG output buffer (on stack).
    ULONG sCube[outputWidthUlong * (subCubes+1)];
    ULONG tmpOut[outputWidthUlong];
    
    // Buffer stores mapping to the term bit positions present in term.
    // Used for mapping from termWeight combinations to numVariables combinations.
    uint termBitPositions[termWeight];
    genBitPosMask(termBitPositions, bitWidth, termMask, termWeight);
    
    // Pre-computation of the key if it is possible, speed optimization.
    // If the last bit position is not in the key block,
    // this optimization can be performed.
    // Key is pre-computed in the calling method in order to avoid race conditions.
    bool precomputedKey = precompKey && (termWeight==0 || (termBitPositions[termWeight-1] < (8*cip->getInputBlockSize())));
    
    // Reset output cubes.
    memset(sCube, 0, SIZEOF_ULONG * outputWidthUlong * (subCubes+1));
    
    // Start computing a cube from the available variables.
    // Global combination counter for parallelization.
    ULONG globalCtr = 0u; 
    // Start from order 0 (constant) and continue up to termWeight.
    for(uint orderCtr=0; orderCtr <= termWeight; orderCtr++){
        // Current order to XOR is orderCtr.
        CombinatiorialGenerator cg(termWeight, orderCtr);
        for(; cg.next(); globalCtr++){
            // Check if is configured to skip this combination.
            if (step>1 && ((globalCtr % step) != offset)) continue;
            // Array of the combination. Each array element contains set 
            // element index chosen for combination.
            const ULONG * comb = cg.getCurState();
            
            // Build the input buffer for the function to evaluate.
            // As the base, use finput variable.
            // Term bits have to be set to zero in finput!
            memcpy(input, finput, byteWidth);
            // Reflect current combination to the finput.
            // Turn bits specified by current combination to 1.
            // Mapping to the bit positions has to be used.
            for(uint tmpOrder=0; tmpOrder<orderCtr; tmpOrder++){
                const uint bitIdx = termBitPositions[comb[tmpOrder]];
                input[bitIdx / 8] |= 1u << (bitIdx % 8);
            }
            
            // Evaluate target function
            if (precomputedKey){
                cip->evaluateWithPreparedKey(input, output);
            } else {
                cip->evaluate(input, input + cip->getInputBlockSize(), output);
            }
            
            // XOR result of this evaluation to the cube block.
            // Each bit corresponds to a different polynomial. For example: f0, f1, ..., f_128 for AES.
            memset(tmpOut, 0, outByteWidth);
            for(uint i=0; i<outByteWidth; i++){
                tmpOut[i/SIZEOF_ULONG] ^= ((ULONG)output[i]) << (8*(i%SIZEOF_ULONG));
            }
            for(uint i=0; i<outputWidthUlong; i++){
                sCube[i] ^= tmpOut[i];
            }
            
            // SubCubes: XOR to the cube that does not contain particular element.
            for(uint subCubeIdx=0; subCubeIdx<subCubes; subCubeIdx++){
                // Subcube 0 is the subcube that does not contain 1 element from the master cube.
                const uchar bitMask = 1u << (termBitPositions[subCubeIdx] % 8);
                const uint  bitPos  = termBitPositions[subCubeIdx] / 8;
                
                if ((input[bitPos] & bitMask) == 0u){      
                    for(uint i=0; i<outputWidthUlong; i++){
                        const uint cbIdx = (subCubeIdx+1)*outputWidthUlong + i;
                        sCube[cbIdx] ^= tmpOut[i];
                    }
                }
            }
        }
    }
    
    memcpy(subcube, sCube, outputWidthUlong * (subCubes+1) * SIZEOF_ULONG);
    return 1;
}

int Approximation::subCubeTermThreaded(uint termWeight, const ULONG* termMask, 
        const uchar* finput, ULONG* subcube, uint subCubes) const
{    
    // Size of the buffer for one thread. Takes subcubes computation into consideration.
    const uint oneThreadBuffSize = outputWidthUlong * (subCubes+1);
    
    // Each thread will have separate ULONG output buffer.
    ULONG * oThreadBuff = new ULONG[oneThreadBuffSize * threadCount];
    
    // Compute cubes in a parallel fashion.
    // Split work among several threads.
    // If only 1 thread is to be used, do not use threads at all.
    // This is good for debugging computation code, since segfault in
    // thread is more difficult to debug.
    memset(subcube, 0, SIZEOF_ULONG * oneThreadBuffSize);
    memset(oThreadBuff, 0, SIZEOF_ULONG * oneThreadBuffSize * threadCount);

    // Pre-compute key, always the same.
    cip->prepareKey(finput + cip->getInputBlockSize());

    std::vector<std::thread> threads;
    for(uint tidx = 0; threadCount > 1 && tidx < threadCount; tidx++){
        // Create & start a new thread with its work partition.
        // Definition of the thread function is using lambda expressions.
        // See http://en.cppreference.com/w/cpp/language/lambda 
        // for lambda expressions.
        threads.push_back(std::thread(
             [=](){
                 this->subCubeTerm(
                         termWeight, termMask, finput,
                         oThreadBuff+(oneThreadBuffSize*tidx),      // Thread output buffer.
                         this->threadCount,                         // Step = thread count.
                         tidx,                                      // Offset = thread index.
                         subCubes,                                  // Compute also n-1 subcubes?
                         false);                                    // Precomputed key?
            }));
    }

    // Join on threads, wait for the result.
    if (threadCount>1){
        for(auto& cthread : threads){
            if (cthread.joinable()){
                cthread.join();
            } else {
                cerr << " .!thread not joinable" << endl;
            }
        }
    } else {
        this->subCubeTerm(termWeight, termMask, finput, oThreadBuff, 1, 0, subCubes, false);
    }

    // Assemble thread results to one single cube result.
    // XOR all resulting sub-cubes into one cube.
    for(uint tidx = 0; tidx < threadCount; tidx++){
        //cout << " ...SubResult = ";
        //dumpHex(cout, oThreadBuff + tidx*outputWidthUlong, outputWidthUlong);

        for(uint octr=0; octr < outputWidthUlong; octr++){
            subcube[octr] ^= oThreadBuff[tidx*oneThreadBuffSize + octr];
        }
    }
    
    delete[] oThreadBuff;
    return 1;
}

int Approximation::keyCube(uint wPlain, uint wKey, uint startOrder, uint stopOrder,
        const ULONG * termMask, std::vector<ULONG> * keyCubes, ULONG * isSuperpoly) const 
{
    const uint numPlainBits = cip->getInputBlockSize() * 8;
    const uint numKeyBits = cip->getKeyBlockSize() * 8;
    
    uchar input[byteWidth];
    ULONG oBuff[outputWidthUlong];    
    ULONG tmpCombination[orderLimit+1];
    ULONG globalCtr = 0;
    for(uint orderCtr=startOrder; orderCtr <= stopOrder; orderCtr++){
        CombinatiorialGenerator cg(numKeyBits, orderCtr);            
        for(; cg.next(); globalCtr++){
            // Array of the key variables.
            const ULONG * keyVars = cg.getCurState();
            // Generate finput for the cube process.
            memset(input, 0, byteWidth);
            // Reflect current key combination to the input.
            // Turn bits specified by current combination to 1.
            // Mapping to the bit positions has to be used.
            for(uint tmpOrder=0; tmpOrder<orderCtr; tmpOrder++){
                const uint bitIdx = numPlainBits + keyVars[tmpOrder];
                input[bitIdx / 8] |= 1u << (bitIdx % 8);
            }

            /*cout << "Starting order=" << orderCtr << "; keyCombinationIdx=" << cg.getCounter() << "; all=" << cg.getTotalNum() <<  endl;
            cout << "Key = ";
            dumpHex(cout, cg.getCurState(), orderCtr);
            cout << "Input = ";
            dumpHex(cout, input, byteWidth);//*/

            // Compute cubes in a parallel fashion. No subcubes allowed.
            if (threadCount>1){
                this->subCubeTermThreaded(wPlain, termMask, input, oBuff, 0);
            } else {
                cip->prepareKey(input + cip->getInputBlockSize());
                this->subCubeTerm(wPlain, termMask, input, oBuff, 1, 0, 0, false);
            }
            
            // XOR with key sub cubes already stored for this plaintext.
            // For example if we compute linear key terms, they has to be XORed
            // with the constant term.
            for(uint xorOrder=0; xorOrder<orderCtr; xorOrder++){
                CombinatiorialGenerator xorOrderGen = CombinatiorialGenerator(orderCtr, xorOrder);
                for(; xorOrderGen.next(); ){
                    for(uint tmpCombCtr = 0; tmpCombCtr<xorOrder; tmpCombCtr++){
                        tmpCombination[tmpCombCtr] = cg.getCurState()[xorOrderGen.getCurState()[tmpCombCtr]];
                    }

                    ULONG idx = combIndexer.getCombinationIdx(xorOrder, tmpCombination);
                    for(uint octr=0; octr < outputWidthUlong; octr++){
                        oBuff[octr] ^= keyCubes[xorOrder][outputWidthUlong*idx + octr];
                    }
                }
            }
            
            // Store the keyCube result for higher cubes computation.
            for(uint octr=0; octr < outputWidthUlong; octr++){
                keyCubes[orderCtr][outputWidthUlong*cg.getCounter() + octr] = oBuff[octr];
            }

            // This relation is interesting only if there is at least one
            // linear relation in it...
            for(uint octr=0; orderCtr>0 && octr<outputWidthUlong; octr++){
                isSuperpoly[octr] |= oBuff[octr];
            }
        }
    }
    
    return 1;
}

std::string Approximation::getCubeCacheName(uint wPlain, uint wKey) const {
    std::string fcacheNameStr = ("cube_a" + to_string(cip->getId()) \
                            + "_r" + to_string(cip->getNumRounds()) \
                            + "_cp" + to_string(wPlain) \
                            + "_ck" + to_string(wKey) \
                            + ".xml");
    return fcacheNameStr;
}

int Approximation::readCubeArchive(const char* fname, CubeRelations_vector& vct) const {
    cout << "Reading cache file=" << fname << endl;
    std::ifstream ifs(fname); 
    if (ifs.good()){
        boost::archive::xml_iarchive ia(ifs);
        ia >> BOOST_SERIALIZATION_NVP(vct);
        return 1;
    }
    
    return 0;
}

int Approximation::writeCubeArchive(const char* fname, CubeRelations_vector& vct) const {
    sigset_t originalMask;
    sigemptyset(&originalMask);
    sigprocmask(SIG_BLOCK, &blockingMask, &originalMask);
    cout << "<save" << flush;

    {
        std::ofstream ofs(fname); //assert(ofs.good());
        boost::archive::xml_oarchive oa(ofs);
        oa << BOOST_SERIALIZATION_NVP(vct);
    }

    cout << "/>" << endl;
    sigprocmask(SIG_SETMASK, &originalMask, NULL);
    
    return 1;
}

int Approximation::writeRelationToArchive(const char* fname, CubeRelations_vector& vct, 
        uint wKey, ULONG* termMask, std::vector<ULONG>* keyCubes, ULONG* isSuperpoly) const 
{
    const uint sizePlaintextUlong = OWN_CEIL((double)cip->getInputBlockSize() / (double)SIZEOF_ULONG);
    
    CubeRelations_t toInsert, *cur;
    std::vector<CubeRelations_t> & keyRelations = vct.get();
    
    keyRelations.push_back(toInsert);
    cur = &(keyRelations.at(keyRelations.size()-1));
    cur->numSuperpolys = hamming_weight_array(isSuperpoly, outputWidthUlong);
    cur->wkey = wKey;

    // Copy plaintext.
    for(uint octr=0; octr < sizePlaintextUlong; octr++){
        cur->termMask.push_back(termMask[octr]);
    }

    // Copy superpoly bitmask.
    for(uint octr=0; octr < outputWidthUlong; octr++){
        cur->isSuperpoly.push_back(isSuperpoly[octr]);
    }

    // Copy superpoly vectors.
    for(uint orderCtr=0; orderCtr <= wKey; orderCtr++){
        cur->superpolys[orderCtr] = keyCubes[orderCtr];
    }

    // Rewrite the archive with the updated keyRelations.
    // Block signals during save in order to avoid storage corruption.
    return writeCubeArchive(fname, vct);
}

ULONG Approximation::dumpCoefficients(std::ostream& c, const std::vector<ULONG>* coefficients, 
        uint maxOrder, uint numVariables, uint polyIdx, uint fmt) const 
{
    const ULONG valMask = ULONG1 << (polyIdx % (SIZEOF_ULONG*8));
    const uint  valPos  =            polyIdx / (SIZEOF_ULONG*8);
    ULONG numterms = 0;
    uchar charBuffer = 0;
    uchar charCtr = 0;
    
    for(uint orderCtr=0; orderCtr <= maxOrder; orderCtr++){
        // Current order to XOR is orderCtr.
        CombinatiorialGenerator cg(numVariables, orderCtr);            
        for(; cg.next(); ){
            // Determine if the current coefficient is enabled for given polynomial.        
            bool isOn = (((coefficients[orderCtr][outputWidthUlong*cg.getCounter() + valPos]) & valMask) != ((ULONG)0));
            numterms+=isOn;
            
            // In binary form, just dump binary coefficient.
            if (fmt==2){
                c << (isOn? "1":"0");
                continue;
            }
            
            // Real binary form.
            if (fmt==3){
                if (isOn){
                    charBuffer |= 1u << (charCtr);
                }                
                
                charCtr+=1;
                if (charCtr>=8){
                    c << charBuffer;
                    charCtr=0u;
                    charBuffer=0u;
                
                }
                
                continue;
            }
            
            // Textual representations, skip.
            if (!isOn){
                continue;
            }
            
            // Term is present.
            if (orderCtr==0){
                c << "1 ";
                continue;
            }
            
            // High order terms.
            const ULONG * keyVars = cg.getCurState();
            for(uint i=0; i<orderCtr; i++){
                c << "x_"<<keyVars[i];
                if ((orderCtr-1) > i) c << "*";
                else c << " ";
            }
        }
    }
    
    // CharBuffer is not flushed completely on purpose since it may contain null 
    // bytes that are not desired to be here for statistical randomness testing.
    if (charCtr==8){
        c << charBuffer;
    }
    
    return numterms;
}

ULONG Approximation::dumpPlaintextCube(std::ostream & c, ULONG * pcube, uint pcubeSize) const {
    uint maxCubeBit = pcubeSize * SIZEOF_ULONG * 8;
    for(uint curBit=0; curBit < maxCubeBit; curBit++){
        const uint bitSeg = curBit / (8*SIZEOF_ULONG);
        const uint bitPos = curBit % (8*SIZEOF_ULONG);
        const ULONG bitMask = ULONG1 << bitPos;
        if ((pcube[bitSeg] & bitMask) != bitMask) continue;
        c << "p_" << dec << curBit << " ";
    }
    
    return 0;
}

ULONG Approximation::dumpOutputFunctions(std::ostream& c, const std::vector<ULONG>* coefficients, 
        uint maxOrder, uint numVariables, uint numPoly, bool nonNullOnly, uint fmt) const 
{
    ULONG totalTerms = 0;
    sigset_t originalMask;
    sigemptyset(&originalMask);
    sigprocmask(SIG_BLOCK, &blockingMask, &originalMask);
    
    for(uint s=0; s<numPoly; s++){
        stringstream ss;
        ULONG cNumTerms = dumpCoefficients(ss, coefficients, maxOrder, numVariables, s, fmt);
        totalTerms += cNumTerms;
        
        if (nonNullOnly && cNumTerms==0) continue;
        if (fmt==1){
            c << " f_="<<setw(3)<< s <<" = ";
            c << ss.str();
            c << endl;
        } else if (fmt==2){
            c << ss.str() << endl;
        } else if (fmt==3){
            c << ss.str();
        }
    }
    
    sigprocmask(SIG_SETMASK, &originalMask, NULL);
    return totalTerms;
}

int Approximation::keyCubePart(uint wPlain, uint wKey, uint orderCtr,
        uint subCubesLimit, uint curSubCube, uint step, uint offset, 
        const ULONG* termMask, std::vector<ULONG>* keyCubes, ULONG* isSuperpoly) const 
{
    const uint numKeyBits   = cip->getKeyBlockSize() * 8;
    const uint numPlainBits = cip->getInputBlockSize() * 8;
    const uint resultSize = outputWidthUlong * (subCubesLimit+1); // Size of the master cube + N x (N-1) sub-cubes result block.
    
    uchar * input   = new uchar[byteWidth];
    ULONG * oBuff   = new ULONG[resultSize];
    ULONG * tmpCombination = new ULONG[orderCtr+1];
    CombinatiorialGenerator ** cgenerators = new CombinatiorialGenerator * [orderCtr];
    ULONG globalCtr = 0;
    double totalComb = 0.0;
    bool doProgressMonitoring = offset==0 && verboseLvl>0 && wPlain>12;
    ProgressMonitor pm(0.01);
    
    // Init combination generators needed for XORing with sub terms.
    for(uint i=0; i < orderCtr; i++){
        cgenerators[i] = new CombinatiorialGenerator(orderCtr, i);
    }
    
    // Current order to XOR is orderCtr.
    CombinatiorialGenerator cg(numKeyBits, orderCtr);      
    totalComb = cg.getTotalNum();
    
    for(; cg.next(); globalCtr++){
        // Parallelization step.
        if (step>1 && ((globalCtr % step) != offset)) continue;
        // Array of the key variables.
        const ULONG * keyVars = cg.getCurState();
        // Generate finput for the cube process.
        memset(input, 0, byteWidth);
        // Reflect current key combination to the input.
        // Turn bits specified by current combination to 1.
        // Mapping to the bit positions has to be used.
        for(uint tmpOrder=0; tmpOrder<orderCtr; tmpOrder++){
            const uint bitIdx = numPlainBits + keyVars[tmpOrder];
            input[bitIdx / 8] |= 1u << (bitIdx % 8);
        }

        /*cout << "Starting order[sub="<<curSubCube<<"]=" << orderCtr << "; keyCombinationIdx=" << cg.getCounter() << "; all=" << cg.getTotalNum() <<  endl;
        cout << "Key = ";
        dumpHex(cout, cg.getCurState(), orderCtr);
        cout << "Input = ";
        dumpHex(cout, input, byteWidth);//*/

        // Compute cubes in a parallel fashion.
        // Code computes also N x (N-1) subcubes directly.
        cip->prepareKey(input + cip->getInputBlockSize());
        this->subCubeTerm(wPlain, termMask, input, oBuff, 1, 0, subCubesLimit, false);
        // Threaded version:
        //this->subCubeTermThreaded(wPlain, termMask, input, oBuff, subCubesLimit);
        
        // XOR with key sub cubes already stored for this plaintext.
        // For example if we compute linear key terms, they has to be XORed
        // with the constant term.
        for(uint xorOrder=0; xorOrder<orderCtr; xorOrder++){
            CombinatiorialGenerator & xorOrderGen = *cgenerators[xorOrder];
            for(xorOrderGen.reset(); xorOrderGen.next(); ){
                for(uint tmpCombCtr = 0; tmpCombCtr<xorOrder; tmpCombCtr++){
                    tmpCombination[tmpCombCtr] = cg.getCurState()[xorOrderGen.getCurState()[tmpCombCtr]];
                }

                ULONG idx = combIndexer.getCombinationIdx(xorOrder, tmpCombination);
                for(uint octr=0; octr < outputWidthUlong; octr++){
                    oBuff[octr] ^= keyCubes[xorOrder][outputWidthUlong*idx + octr];
                }

                // Subsubes
                for(uint subCubeIdx=0, sOffset=wKey+curSubCube+1; 
                        subCubeIdx < subCubesLimit; 
                        subCubeIdx++, sOffset+=wKey+2+curSubCube)
                {
                    for(uint octr=0; octr < outputWidthUlong; octr++){
                        oBuff[outputWidthUlong*(subCubeIdx+1) + octr] ^= keyCubes[sOffset + xorOrder][outputWidthUlong*idx + octr];
                    }
                }
            }
        }

        // Store the keyCube result for higher cubes computation.
        for(uint octr=0; octr < outputWidthUlong; octr++){
            keyCubes[orderCtr][outputWidthUlong*cg.getCounter() + octr] = oBuff[octr];
        }

        // Subsubes
        for(uint subCubeIdx=0, sOffset=wKey+curSubCube+1; 
                subCubeIdx < subCubesLimit; 
                subCubeIdx++, sOffset+=wKey+2+curSubCube)
        {
            for(uint octr=0; octr < outputWidthUlong; octr++){
                keyCubes[sOffset + orderCtr][outputWidthUlong*cg.getCounter() + octr] = oBuff[outputWidthUlong*(subCubeIdx+1) + octr];
            }
        }

        //
        // This relation is interesting only if there is at least one
        // linear relation in it...
        for(uint octr=0; orderCtr>0 && octr<resultSize; octr++){
            isSuperpoly[octr] |= oBuff[octr];
        }
        
        // Progress monitoring.
        if (doProgressMonitoring){
            pm.setCur((double)cg.getCounter() / totalComb);
        }
    }
    
    // Progress monitoring.
    if (doProgressMonitoring){
        pm.setCur(1.0);
        cout << endl;
    }
    
    // Destroy combination generators.
    for(uint i=0; i<orderCtr; i++){
        if (cgenerators[i]!=NULL){
            delete cgenerators[i];
            cgenerators[i] = NULL;
        }
    }
    
    delete[] input;
    delete[] oBuff;
    delete[] tmpCombination;
    delete[] cgenerators;
    
    return 1;
}

int Approximation::cubeAttack(uint wPlain, uint wKey, uint numRelations, uint subCubesLimit, 
        bool saveRelations, bool dumpAllRelations) const 
{
    const uint bitWidth = 8*byteWidth;
    const uint numKeyBits = cip->getKeyBlockSize() * 8;
    
    // Function input/output for evaluation.
    bool doProgressMonitoring = verboseLvl>0 && wPlain>12;
    uchar * output = new uchar[cip->getOutputBlockSize()];
    uchar * key    = new uchar[cip->getKeyBlockSize()];
    uchar * input  = new uchar[byteWidth];
    uchar * solvedKey = new uchar[cip->getKeyBlockSize()];
    // Number of subcubes to compute. 
    uint subCubes = subCubesLimit;
    assert(wPlain >= subCubes);
    // Size of the master cube + N x (N-1) sub-cubes result block.
    const uint resultSize = outputWidthUlong * (subCubes+1);
    // Local ULONG output buffer (on stack).
    ULONG * oBuff = new ULONG[resultSize];
    // Is superpoly present in the given cube?
    ULONG * isSuperpoly = new ULONG[resultSize];
    // Plaintext buffer.
    ULONG * termMask = new ULONG[inputWidthUlong];
    ULONG * subTermMask = new ULONG[inputWidthUlong];
    // SubCubes for key bits - result of approximation for one particular plaintext.
    std::vector<ULONG> * keyCubes = new std::vector<ULONG>[wKey+1+subCubes*(wKey+2)];
    // Interesting key bits relations are stored in this array, serializes to a file.
    CubeRelations_vector keyRelationsVector;
    // Random variable bit vector = for constructing random terms.
    // In cube attack we are using plaintext variables.
    std::vector<uint> vars;
    for(uint i=0; i<(cip->getInputBlockSize()*8); i++){
        vars.push_back(i);
    }
    ULONG * tmpCombination = new ULONG[orderLimit+1];
    
    cout << " Going to start cube, threads=" << threadCount << endl;
    uint nonzeroCoutner=0;
    uint totalRelations=0;
    
    // Cube attack caches found relations to the file since it is not that fast 
    // to find them and in order to enable interrupted computation.
    std::string fcacheNameStr = getCubeCacheName(wPlain, wKey);
    if (saveRelations){
        readCubeArchive(fcacheNameStr.c_str(), keyRelationsVector);
    }
    
    // Binary storage of the relations found.
    std::string relFile = getCubeCacheName(wPlain, wKey) + ".bin";
    ofstream relOf(relFile.c_str(), ios_base::app);
    
    std::vector<CubeRelations_t> & keyRelations = keyRelationsVector.get();
    nonzeroCoutner = keyRelations.size();
    totalRelations = keyRelationsVector.getTotal();
    cout << " Loaded " << nonzeroCoutner << " relations from the file; total=" << totalRelations << endl;
    
    // Generate multiple relations for the key variables from the black-box function.
    for(uint relationIdx = 0; nonzeroCoutner < numRelations; relationIdx++){
        if (wPlain>7 || true){
            cout << "." << flush;
        }
        
        // Collect relations for the key variables.
        // For this current relation, we generate randomly wPlain-bit 
        // plaintext. 
        random_shuffle(std::begin(vars), std::end(vars));
        for(uint subCube=0; subCube<=0; subCube++){                 // TODO: for now no subcubes, code artefact.
            memset(termMask, 0, SIZEOF_ULONG * inputWidthUlong);
            for(uint bitPos=0; bitPos < (wPlain-subCube); bitPos++){
                const uint bitIdx = vars[bitPos];
                termMask[bitIdx / (8*SIZEOF_ULONG)] |= ULONG1 << (bitIdx % (8*SIZEOF_ULONG));
            }
            
            // Convert to num. bit representation, needed for sub-cubes computation.
            uint termBitPositions[wPlain];
            genBitPosMask(termBitPositions, bitWidth, termMask, wPlain);

            // During cube computation on the key, you have to compute also cube
            // on the keys, so memorize lower key cubes for computation of the higher
            // key cubes, note, all the time plaintext is the same.
            // Initialize keycube vector space.
            for(uint order = 0; order<=(wKey+subCube); order++){
                ULONG vecSize = CombinatiorialGenerator::binomial(numKeyBits, order) * outputWidthUlong;
                keyCubes[order].assign(vecSize, (ULONG) 0);
            }
            
            // Allocation for subcubes.
            for(uint subCubeIdx=0, sOffset=wKey+subCube+1; 
                    subCubeIdx < subCubes; 
                    subCubeIdx++, sOffset+=wKey+2+subCube)
            {
                for(uint order = 0; order<=(wKey+1+subCube); order++){
                    ULONG vecSize = CombinatiorialGenerator::binomial(numKeyBits, order) * outputWidthUlong;
                    keyCubes[sOffset + order].assign(vecSize, (ULONG) 0);
                }
            }

            // Collects superpolys according to the paper, reset for each plaintext.
            memset(isSuperpoly, 0, SIZEOF_ULONG * resultSize);
            
            // For each key, one plaintext cube has to be computed.
            for(uint orderCtr=0; orderCtr <= (wKey+subCube); orderCtr++){ // TODO: subcube computation is not correct!
                if (doProgressMonitoring){
                    cout << " r=" << relationIdx 
                            << "; sub=" << subCube 
                            << "; orderCtr=" << orderCtr 
                            << "; time=" << currentDateTime() << endl;
                }
                
                // Current order to XOR is orderCtr.
                // If order is 0, do not parallelize.
                uint threadNum = orderCtr==0 ? 1 : this->threadCount;
                if (threadNum==1){
                    this->keyCubePart(
                        wPlain, wKey, orderCtr,
                        subCubes, subCube, threadNum, 0, 
                        termMask, keyCubes, isSuperpoly);
                } else {
                    std::vector<std::thread> threads;
                    for(uint tidx = 0; tidx < threadNum; tidx++){
                        // Create & start a new thread with its work partition.
                        // Definition of the thread function is using lambda expressions.
                        // See http://en.cppreference.com/w/cpp/language/lambda 
                        // for lambda expressions.
                        threads.push_back(std::thread(
                            [=](){
                                 this->keyCubePart(
                                        wPlain, wKey, orderCtr,
                                        subCubes, subCube, threadNum, tidx, 
                                        termMask, keyCubes, isSuperpoly);
                            }));
                    }

                    // Join on threads, wait for the result.
                    for(auto& cthread : threads){
                        if (cthread.joinable()){
                            cthread.join();
                        } else {
                            cerr << " .!thread not joinable" << endl;
                        }
                    }
                }
                
                // Optimization: do not compute quadratic & quartic tests
                // if there is no linear term. Fast optimization, linear terms
                // are important, and its computation is almost for free. If we have
                // them, we have to check quadratic/quartic terms, but this computation
                // is very expensive since we cannot recycle previous cubes anymore.
                uint numSuperpolys=hamming_weight_array(isSuperpoly, outputWidthUlong);
                if (orderCtr > 0 && subCube>0 && numSuperpolys==0){
                    cout << "/" << flush;
                    break;
                }
            }

            // Superpoly relations are stored to the keyRelations structure.
            // Format is prescribed by the serialization of the class (uses BOOST serialization):
            uint numSuperpolys=hamming_weight_array(isSuperpoly, outputWidthUlong);
            uint totalSuperpolys=hamming_weight_array(isSuperpoly, resultSize);
            if (totalSuperpolys>0){
                cout << "|" << numSuperpolys << ":" << totalSuperpolys << "[";
                for(uint i=0; i<subCubes; i++){
                    uint cc = hamming_weight_array(isSuperpoly+outputWidthUlong*(i+1), outputWidthUlong);
                    cout << cc << ",";
                } cout << "]" << flush;
            }
            
            // Save relation from the master cube.
            if (numSuperpolys > 0 || dumpAllRelations){
                cout << endl
                        << "    ----- Have superpoly --------; sub=" << subCube
                        << "; num=" << numSuperpolys
                        << "; soFar=" << nonzeroCoutner 
                        << "; total=" << totalRelations << endl;

                nonzeroCoutner+=1;
                totalRelations+=numSuperpolys;
                keyRelationsVector.setTotal(totalRelations);
                
                if (verboseLvl>1){
                    //dumpOutputFunctions(cout, keyCubes, wKey, numKeyBits, cip->getOutputBlockSize()*8, true, 2);
                    dumpOutputFunctions(relOf, keyCubes, wKey, numKeyBits, cip->getOutputBlockSize()*8, false, 3);
                }
                
                // Write to the archive.
                if (saveRelations){
                    writeRelationToArchive(fcacheNameStr.c_str(), keyRelationsVector, 
                        wKey+subCube, termMask, keyCubes, isSuperpoly);
                }
            }
            
            // Need to compute quadratic parts of the N-1 sized subcubes.
            // This is done only if subcube already contains at least one linear term.
            if(totalSuperpolys>0){
                for(uint sCubeIdx=0, sOffset=wKey+subCube+1; 
                        sCubeIdx < subCubes; 
                        sCubeIdx++, sOffset+=wKey+2+subCube)
                {
                    const uint isSuperpolyOffset = (sCubeIdx+1)*outputWidthUlong;
                    uint weight = hamming_weight_array(isSuperpoly+isSuperpolyOffset, outputWidthUlong);
                    if (weight<=1) continue;
                    
                    cout << "[" << weight << flush << endl;
                    
                    // Cancel given bit position in this sub mask
                    memcpy(subTermMask, termMask, SIZEOF_ULONG * inputWidthUlong);
                    subTermMask[termBitPositions[sCubeIdx]/(8*SIZEOF_ULONG)] &= ~(ULONG1<<(termBitPositions[sCubeIdx] % (8*SIZEOF_ULONG))); 
                    
                    // Compute quadratic key cubes.
                    const uint wSubPlain = wPlain-1-subCube;
                    assert(hamming_weight_array(subTermMask, outputWidthUlong) == wSubPlain);
                    
                    keyCube(wSubPlain, wKey+1+subCube, 2+subCube, 2+subCube, subTermMask, keyCubes+sOffset, isSuperpoly+isSuperpolyOffset);
                    
                    // Number of quads?
                    // Syso kontrapriklad:
                    ULONG termCount = CombinatiorialGenerator::binomial(numKeyBits, 2);
                    ULONG sysoCtr = 0;
                    for(uint syso=0; syso<termCount; syso++){
                        sysoCtr += hamming_weight_array(keyCubes[sOffset+2], syso*outputWidthUlong, outputWidthUlong);
                    }
                    
                    // Dumps subcube's polynomials in a readable form.
                    dumpOutputFunctions(cout, keyCubes + sOffset, 2, numKeyBits, cip->getOutputBlockSize()*8, true);
                    
                    cout << "2We have " << sysoCtr  
                            << " sysovin, total=" << ((double)sysoCtr / (termCount*outputWidthUlong*SIZEOF_ULONG*8.0)) 
                            << "; plaintext cube: ";
                    dumpPlaintextCube(cout, subTermMask, inputWidthUlong);
                    cout << endl;
                    
                    // Save relations.
                    nonzeroCoutner+=1;
                    totalRelations+=weight;
                    keyRelationsVector.setTotal(totalRelations);
                
                    // Write to the archive.
                    if (saveRelations){
                        writeRelationToArchive(fcacheNameStr.c_str(), keyRelationsVector, 
                            wKey+1+subCube, subTermMask, keyCubes+sOffset, isSuperpoly+isSuperpolyOffset);
                    }
                    
                    cout << "]" << flush;
                } 
                cout << endl;
            }
            
            // No additional sub-cubes. Not needed.
            if (totalSuperpolys>0){
                break;
            }
        }
    }
    
    cout << "Number of nonzero equations: " << dec << nonzeroCoutner << endl;
    
    // Process found relations and try to solve the system online.
    // Online-phase strategy: select random key (that we want to discover).
    randomBuffer(key, cip->getKeyBlockSize());
    // Prepare input - paste key into it.
    memset(input, 0, byteWidth);
    memcpy(input + cip->getInputBlockSize(), key, cip->getOutputBlockSize());
    cout << "Input ready="; dumpHex(cout, input, byteWidth);
    cout << "GenKey="; dumpHex(cout, key, cip->getOutputBlockSize());
     
    // Online phase of the attack. Key is fixed, goal is to determine it.
    // We have to determine b_t for each relation:
    //
    // \Sum_{v \in C_t} f(v,x) = b_t
    // a_1x_1 + a_2x_2 + \dots + a_nx_n + c = b_t       (this is for one relation)
    // 
    // From this we get system of n variables and more than n equations we want 
    // to solve with GB or Gaussian elimination.
    // Solve the system that is over-determined. There are probably a lot of 
    // linearly dependent equations in the base.
    long rank = cubeOnlineAttack(keyRelationsVector, input, solvedKey);
    if (rank >= numKeyBits){
        cout << "   We have the solution:" << endl;
        dumpHex(cout, solvedKey, cip->getKeyBlockSize());
        dumpHex(cout, key, cip->getKeyBlockSize());
        
        // Compute bit-hit ratio.
        uint matches = numBitMatches(solvedKey, key, 0, numKeyBits);
        cout << "Matches=" << matches << " that is " << (((double)matches/(double)numKeyBits)*100) << " % " << endl; 
    }
    
    delete[] input;
    delete[] key;
    delete[] output;
    delete[] oBuff;
    delete[] isSuperpoly;
    delete[] solvedKey;
    delete[] termMask;
    delete[] keyCubes;
    delete[] tmpCombination;
    delete[] subTermMask;
    return 1;
}

long Approximation::cubeOnlineAttack(CubeRelations_vector & keyRelationsVector, const uchar * input, uchar * solvedKey) const {
    const uint numKeyBits = cip->getKeyBlockSize() * 8;
    const uint sizePlaintextUlong = OWN_CEIL((double)cip->getInputBlockSize() / (double)SIZEOF_ULONG);
    
    ULONG * termMask      = new ULONG[inputWidthUlong];
    ULONG * oBuff         = new ULONG[outputWidthUlong];
    uint totalRelations   = keyRelationsVector.getTotal();
    std::vector<CubeRelations_t> & keyRelations = keyRelationsVector.get();
    
    // Solution vector
    vec_GF2 solution;
    // Right side of the system is stored here. It is c^b_t.
    vec_GF2 b(INIT_SIZE, totalRelations);
    // System of the equations to solve.
    mat_GF2 systm(INIT_SIZE, totalRelations, numKeyBits);
    // Working matrix for Gaussian elimination.
    mat_GF2 systmGauss(INIT_SIZE, totalRelations, numKeyBits+1);
    // Add polynomial to the systm matrix, constant term to the b vector.
    uint curRow=0;
    for (std::vector<CubeRelations_t>::iterator it = keyRelations.begin() ; it != keyRelations.end(); ++it){
        // Copy plaintext
        memset(termMask, 0, SIZEOF_ULONG * inputWidthUlong);
        for(uint tidx=0; tidx<sizePlaintextUlong; tidx++){
            termMask[tidx] = it->termMask[tidx];
        }
        
        uint wCurPlain = hamming_weight_array(termMask, inputWidthUlong);
        if (verboseLvl>2){
            cout << "Solving plaintext cube["<<wCurPlain<<"]="; 
            dumpHex(cout, termMask, inputWidthUlong);
        }
        
        // Has to cube f(v,x) for C_t in order to obtain b_t.
        this->subCubeTermThreaded(wCurPlain, termMask, input, oBuff, 0);
        
        // For each polynomial present
        for(uint polyIdx=0; polyIdx < numKeyBits; polyIdx++){
            // b_t is stored in vectorized form (for each f_i) in oBuff.
            std::string desc = "";
            uint polySegm  = polyIdx / (8*SIZEOF_ULONG);
            ULONG polyMask = ULONG1 << (polyIdx % (8*SIZEOF_ULONG));
            
            // Only for valid superpoly - there is at least one linear relation.
            if ((it->isSuperpoly[polySegm] & polyMask) != polyMask){
                continue;
            }
            
            const bool b_t = (oBuff[polySegm] & polyMask) == polyMask;
            const bool c   = ((it->superpolys[0][polySegm]) & polyMask) == polyMask;
            b.put(curRow, c ^ b_t);

            // Now init the system of equations.
            for(uint varIdx=0; varIdx < numKeyBits; varIdx++){
                const bool ai = ((it->superpolys[1][varIdx*outputWidthUlong+polySegm]) & polyMask) == polyMask;
                systm.put(curRow, varIdx, ai);
                desc += ai ? "+" : "_";
            }
            
            if (verboseLvl>1){
                cout << " p[" << setw(4) << right << polyIdx << "]=" << desc << endl;
            }
            
            curRow+=1;
        }
    }
    
    // Solve the system that is over-determined. There are probably a lot of 
    // linearly dependent equations in the base.
    long rank = gaussPh4r05(systmGauss, systm, b, numKeyBits);
    if (verboseLvl>1){
        cout << "Gaussian elimination performed on the system matrix, rank=" << dec << rank << "; totalRelations=" << totalRelations << endl;
    }
    
    if (rank < numKeyBits){
        if (verboseLvl>1){
            cout << " Cannot solve the system, rank is too low" << endl;
        }
    } else {
        // Transfer first numKeyBits rows to the systm matrix and b vector.
        solution.SetLength(numKeyBits);
        for(uint idx=0; idx < numKeyBits; idx++){
            solution.put(idx, systmGauss.get(idx, numKeyBits));
        }
        
        // We have solution, convert it to the uchar array and dump in hexa and in the binary.
        const uint solByteSize = (uint)ceil(numKeyBits/8.0);
        memset(solvedKey, 0x0, solByteSize);
        for(uint i=0; i<numKeyBits; i++){
            if (IsOne(solution.get(i))){
                solvedKey[i/8] |= 1u << (i%8);
            }
        }
    }
    
    delete[] termMask;
    delete[] oBuff;
    return rank;
}

const std::string Approximation::currentDateTime() const {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);

    return buf;
}


void Approximation::initFGb(uint numVariables) const {
    fgb.initFGb(numVariables);
}

void Approximation::deinitFGb() const {
    fgb.deinitFGb();
}

void Approximation::resetFGb() const {
    fgb.resetFGb();
}
