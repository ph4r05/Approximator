#include "base.h"
#include <iostream>
#include <iomanip>
NTL_CLIENT

using namespace std;
using namespace NTL;

void dumpUchar(std::ostream & c, const uchar* inp, unsigned int size, bool endl) {
    for(unsigned int i = 0; i < size; i++){
        c << (unsigned int) inp[i] << " ";
    }
    
    if (endl){
        c << std::endl;
    }
}

void readUcharToUlong(const uchar * input, uint size, ULONG * iBuff) {
    // At first reset memory of the big buffer.
    memset(iBuff, 0, size);
    // And read particular parts of the input to the big buffer.
    for(uint x = 0; x < size; x++){
        iBuff[x/SIZEOF_ULONG] = READ_TERM_1(iBuff[x/SIZEOF_ULONG], input[x], x%SIZEOF_ULONG);
    }
}

void readUlongToUchar(uchar* output, uint size, const ULONG* iBuff) {
    for(uint x=0; x<size; x++){
        output[x] = (iBuff[x/SIZEOF_ULONG] >> (8* (x % SIZEOF_ULONG))) & ((unsigned char)0xffu);
    }
}

void randomBuffer(uchar * buffer, uint size){
    for(unsigned int k=0; k<size; k++){ 
        buffer[k] = (rand() % (0x100ul)); 
    }
}

uint numBitMatches(const uchar * a, const uchar * b, uint bitPosStart, uint bitPosEnd, uint offsetA, uint offsetB){
    uint hits = 0;
    for(uint i=bitPosStart; i<bitPosEnd; i++){
        const uint keyIdxA = i+offsetA;
        const uint keyIdxB = i+offsetB;
        hits += ((a[keyIdxA/8] & (1u << (keyIdxA%8))) > 0) == ((b[keyIdxB/8] & (1u << (keyIdxB%8))) > 0);
    }

    return hits;
}

long gaussPh4r05(mat_GF2& M, const mat_GF2& A, const vec_GF2 & b, long w)
{
   long k, l;
   long i, j;
   long pos;

   long n = A.NumRows();
   long m = A.NumCols()+1;

   if (w < 0 || w > n)
      Error("gauss: bad args");

   M.SetDims(n, m);
   vec_GF2 aa;
   aa.SetLength(m);

   for (i = 0; i < n; i++) {
      aa = A[i];
      aa.SetLength(m);
      aa.put(m-1, b.get(i));
      M[i] = aa;
   }
   
   long wm = (m + NTL_BITS_PER_LONG - 1)/NTL_BITS_PER_LONG;   
   
   l = 0;
   for (k = 0; k < w && l < n; k++) {
      long wk = k/NTL_BITS_PER_LONG;
      long bk = k - wk*NTL_BITS_PER_LONG;
      _ntl_ulong k_mask = 1UL << bk;


      pos = -1;
      for (i = l; i < n; i++) {
         if (M[i].rep.elts()[wk] & k_mask) {
            pos = i;
            break;
         }
      }

      if (pos != -1) {
         if (l != pos)
            swap(M[pos], M[l]);

         _ntl_ulong *y = M[l].rep.elts();

         // Null rows below and above.
         for (i = 0; i < n; i++) {
            // M[i] = M[i] + M[l]*M[i,k]
            if (i==l) continue;
            
            if (M[i].rep.elts()[wk] & k_mask) {
               _ntl_ulong *x = M[i].rep.elts();

               for (j = wk; j < wm; j++)
                  x[j] ^= y[j];
            }
         }

         l++;
      }
   }
   
   return l;
}


