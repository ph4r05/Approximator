/* 
 * File:   FGbHelper.cpp
 * Author: ph4r05
 * 
 * Created on May 28, 2014, 11:51 AM
 */

#include "FGbHelper.h"
#include <iostream>
#include <iomanip>
#include <openssl/md5.h>

#include "CombinatiorialGenerator.h"

// The following macro should be 1 to call FGb modulo a prime number.
#define LIBMODE 1  
#define CALL_FGB_DO_NOT_DEFINE
#include "call_fgb.h"

using namespace std;

FGbHelper::FGbHelper() : varNames(NULL), fgbFile(NULL), orderLimit(0), byteWidth(0),
    outputWidthUlong(0) {
}
    
FGbHelper::~FGbHelper() {
    deinit();
}

void FGbHelper::init(ULONG byteWidth, uint orderLimit, uint outputBits) {
    this->byteWidth = byteWidth;
    this->orderLimit = orderLimit;
    this->outputWidthUlong = OWN_CEIL((double)outputBits / (double)(8.0 * SIZEOF_ULONG));
    
    // Variable names for FGb.
    varNames = new char * [8*byteWidth];
    for(uint var=0; var<8*byteWidth; var++){
        varNames[var] = new char[10];
        snprintf(varNames[var], 10, "x_%03d", var);
    }
    
    // Open logfile for library
    fgbFile = fopen("fgb.log", "a+");
}

void FGbHelper::deinit() {
    if (varNames!=NULL){
        for(uint var = 0; var <= 8*byteWidth; var++){
            delete[] varNames[var];
            varNames[var]=NULL;
        }
        
        delete[] varNames;
        varNames = NULL;
    }
    
    if (fgbFile!=NULL){
        fclose(fgbFile);
        fgbFile = NULL;
    }
}

I32 FGbHelper::getNumberOfTerms(Dpol poly) const {
    return FGB(nb_terms)(poly);
}

bool FGbHelper::isPolyNull(Dpol poly) const {
    return getNumberOfTerms(poly)==0;
}

bool FGbHelper::isPoly1(Dpol poly, I32 numVariables) const {
    const I32 nbTerms = getNumberOfTerms(poly);
    if (nbTerms!=1) return false;
    
    I32 cfs = 0;
    I32 mons[numVariables];
    
    FGB(export_poly)(numVariables, nbTerms, mons, &cfs, poly);
    if (cfs!=1) return false;
    
    {
        bool isOne=true;
        for(int i=0; i<numVariables; i++){
            if (mons[i]){
                isOne=false;
                break;
            }
        }
        
        return isOne;
    }
}

I32 FGbHelper::exportPolynomial(I32 numVariables, I32 numTerms, I32* exps, I32* coefs, Dpol polynomial) const {
    return FGB(export_poly)(numVariables, numTerms, exps, coefs, polynomial);
}

Dpol_INT FGbHelper::polynomial2FGb(uint numVariables, std::vector<ULONG>* coefs, uint maxOrder, uint polyIdx, ULONG * numTerms, ULONG * hash) const {
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
        if (hash!=NULL){
            MD5_Update(&md5Ctx, &order, sizeof(uint));
            MD5_Update(&md5Ctx, "<HASH-NEW-ORDER>", sizeof(char));
        }
        
        for(; cg.next() && termCounter < termsEnabled; ){
            const ULONG ctr = cg.getCounter();
            if ((coefs[order][outputWidthUlong*ctr+coefSegment] & (ULONG1<<(coefOffset))) == 0) continue;
            
            // Update hash value if we want to compute it. 
            if (hash!=NULL){
                MD5_Update(&md5Ctx, &ctr, SIZEOF_ULONG);
                MD5_Update(&md5Ctx, "<HASH-NEW-TERM>", sizeof(char));
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

void FGbHelper::dumpFGbPoly(uint numVariables, Dpol poly) const {
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

void FGbHelper::dumpBasis(uint numVariables, Dpol* basis, uint numPoly) const {
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

void FGbHelper::initFGb(uint numVariables) const {
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

void FGbHelper::deinitFGb() const {
    resetFGb();
    FGB(exit)(); /* restore original GMP allocators */
}

void FGbHelper::resetFGb() const {
    FGB(reset_memory)(); /* to reset Memory */
}

int FGbHelper::computeFGb(int n_input, Dpol* inputBasis, Dpol* outputBasis, double* t0) const {
    int step0=-1;
    int bk0=0;
    int nb;
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
    nb=FGB(groebner)(inputBasis,n_input,outputBasis,1,0,t0,bk0,step0,0,env);
    return nb;
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

