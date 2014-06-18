/* 
 * File:   base.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 2:06 PM
 */

#ifndef BASE_H
#define	BASE_H
#include <iostream>
#include <iomanip>
#include <cstring>
#include <bitset>
#include <vector>
#include <random>       // std::default_random_engine
#include <algorithm>    // std::move_backward
#include <array>
#include <iterator>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>

typedef unsigned char uchar;
typedef unsigned int  uint;
typedef unsigned long ulong;

// Determine whether we are building for a 64-bit platform.
// _LP64: http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
#if defined(_M_X64) || defined(__amd64__) || defined(_LP64) || defined(_ILP64)
#define COMPILER_X64
#endif

// Define ulong type.
#if defined(COMPILER_X64)
// Compilation for 64 bit platform.
typedef unsigned long long ULONG;
#define SIZEOF_ULONG 8
#define FULL_ULONG 0xffffffffffffffffull
#define ULONG1 ((ULONG) 1ull)
#else 
typedef unsigned long ULONG;
#define SIZEOF_ULONG 4
#define FULL_ULONG 0xfffffffful
#define ULONG1 1ul
#endif

// Shifts term to the right 8*shift, and takes the lower 8 bits (corresponds to the input
// representation with unsigned char).
#define MASK_TERM(trm, shift) (((trm)>>(8*(shift))) & 0xfful)

// Generates 8bit mask with the given offset (in terms of bytes)
#define GET_MASK(offset) (((ULONG) 0xfful) << (8*(offset)))

// Reads input (uchar) to the term on the given offset.
// Offset is in terms of bytes.
// Value produced is term as befure, but on the correct place input is substituted.
#define READ_TERM_1(trm, input, offset) (((trm) & (~(GET_MASK(offset)))) | (((ULONG)(input))<<((offset)*8)))

// Evaluate internal representation of a term, using correct size of an internal type.
#if defined(COMPILER_X64)
// Compilation for 64 bit platform.
#else 
// Determine if long is of length 4 B
#if defined(__SIZEOF_LONG__) && __SIZEOF_LONG__ == 4
// Assume we are compiling for 32 bit platform, if macro for size of long is defined,
// and is of size 4 bytes, use macro to evaluate 4 B type.
#endif
#endif 

// Using default implementation by function call. 
// System is neither x86_64 nor GCC having __SIZEOF_LONG__ set to 4.
//#ifndef 
//#define TERM_ITEM_EVAL_GENOME(trm, input)  
//#endif

// Fast ceiling function for integers.
#define OWN_CEIL(x)  (    (((int)(x)) < (x)) ? ((int)(x))+1 : ((int)(x))    )
#define OWN_FLOOR(x) (    (((int)(x)) < (x)) ? ((int)(x))-1 : ((int)(x))    )

// Maximal base size for FGb (number of polynomials).
#define FGb_MAXI_BASE 100000

/**
 * Reads 8-bit buffer to the 64 bit buffer.
 * iBuff has to be big enough to fit the input buffer.
 * 
 * @param input     input buffer to read.
 * @param size      size of the input buffer in bytes to read.
 * @param iBuff     destination buffer.
 */
void readUcharToUlong(const uchar * input, uint size, ULONG * iBuff);

/**
 * Reads 64 bit buffer to the 8 bit buffer.
 * 
 * @param output    output buffer to write.
 * @param size      size of the input buffer to read in bytes.
 * @param iBuff     input buffer to read (from LSB).
 */
void readUlongToUchar(uchar * output, uint size, const ULONG * iBuff);

/**
 * Fills given buffer with random values.
 * @param buffer
 * @param size
 */
void randomBuffer(uchar * buffer, uint size);

/**
 * Computes Hamming weight of the numeric type. 
 * Nifty parallel counting. 
 * 
 * For more inspiration take a look at http://bisqwit.iki.fi/source/misc/bitcounting/
 * @return 
 */
template<class TestType>
uint hamming_weight_fast(TestType k){
    TestType n = k;
    const unsigned TEST_BITS = sizeof(TestType) * 8;
    TestType m1 = (~(TestType)0) / 3;   // Binary 01010101...
    TestType m2 = (~(TestType)0) / 5;   // Binary 00110011...
    TestType m4 = (~(TestType)0) / 17;  // Binary 00001111...
    TestType h01 = (~(TestType)0) / 255; // Hex 0101...

    n = (n & m1) + ((n >> 1) & m1);
    n = (n & m2) + ((n >> 2) & m2);
    n = (n & m4) + ((n >> 4) & m4);

    return (n * h01) >> (TEST_BITS-8);
    // ^incidentally same as n % 255
}

/**
 * Computes Hamming weight of the numeric type. 
 * For more inspiration take a look at http://bisqwit.iki.fi/source/misc/bitcounting/
 * 
 * @param N
 * @return 
 */
template<class T>
uint hamming_weight(T n){
    uint result=0;
    while(n){
        result++;
        n &= n-1;   // Zero the lowest-order one-bit
    }
    return result;
}

/**
 * Computes hamming weight of the numeric array
 * @param n
 * @param size
 * @return 
 */
template<class T>
uint hamming_weight_array(const T * arr, uint size){
    uint result=0;
    for(uint i=0; i<size; i++){
        T n = arr[i];
        while(n){
            result++;
            n &= n-1;   // Zero the lowest-order one-bit
        }
    }
    return result;
}

template<class T>
uint hamming_weight_array(const std::vector<T> & inp, uint offset, uint size){
    uint result=0;
    for(uint i=0; i<size; i++){
        T n = inp[offset+i];
        while(n){
            result++;
            n &= n-1;   // Zero the lowest-order one-bit
        }
    }
    return result;
}

template<class T>
void dumpHex(std::ostream & c, const std::vector<T> & inp, unsigned int size, bool endl=1) {
    c << std::showbase // show the 0x prefix
      << std::internal // fill between the prefix and the number
      << std::setfill('0'); // fill with 0s
    
    for(unsigned int i = 0; i < size; i++){
        const ULONG toDisp = static_cast<const ULONG>(inp[i]);
        if (toDisp==0){
            c << std::hex << "0x" << std::setw(2 * sizeof(T)) << 0 << " ";
        } else {
            c << std::hex << std::setw(2 * sizeof(T)+2) << toDisp << " ";
        }
    }
    
    if (endl){
        c << std::endl;
    }
}

template<class T>
void dumpHex(std::ostream & c, const T * inp, unsigned int size, bool endl=1) {
    c << std::showbase // show the 0x prefix
      << std::internal // fill between the prefix and the number
      << std::setfill('0'); // fill with 0s
    
    for(unsigned int i = 0; i < size; i++){
        const ULONG toDisp = static_cast<const ULONG>(inp[i]);
        if (toDisp==0){
            c << std::hex << "0x" << std::setw(2 * sizeof(T)) << 0 << " ";
        } else {
            c << std::hex << std::setw(2 * sizeof(T)+2) << toDisp << " ";
        }
    }
    
    if (endl){
        c << std::endl;
    }
}

template<class T>
void dumpBin(std::ostream & c, const T * inp, unsigned int size, bool endl=1) {
    for(unsigned int i = 0; i < size; i++){
        const ULONG toDisp = static_cast<const ULONG>(inp[i]);
        std::bitset<sizeof(T)*8> x(toDisp);
        c << std::dec << std::setw(8 * sizeof(T)) << x << "b ";
    }
    
    if (endl){
        c << std::endl;
    }
}

template<class T>
void dumpBin(std::ostream & c, const std::vector<T> & inp, unsigned int size, bool endl=1) {    
    for(unsigned int i = 0; i < size; i++){
        const ULONG toDisp = static_cast<const ULONG>(inp[i]);
        std::bitset<sizeof(T)*8> x(toDisp);
        c << std::dec << std::setw(8 * sizeof(T)) << x << "b ";
    }
    
    if (endl){
        c << std::endl;
    }
}

inline void dumpUcharHex(std::ostream & c, const uchar* inp, unsigned int size, bool endl=1){ dumpHex(c, inp, size, endl); }
inline void dumpUlongHex(std::ostream & c, const ULONG* inp, unsigned int size, bool endl=1){ dumpHex(c, inp, size, endl); }
void dumpUchar   (std::ostream & c, const uchar * inp, unsigned int size, bool endl=1);

/**
 * Computes number of bits that matches.
 * 
 * @param a
 * @param b
 * @param bitPosStart
 * @param bitPosEnd
 * @return 
 */
uint numBitMatches(const uchar * a, const uchar * b, uint bitPosStart, uint bitPosEnd, uint offsetA=0, uint offsetB=0);

/**
 * Solves linear system. A does not have to be a square matrix.
 * M = [A|c].
 * 
 * After algorithm finishes, M contains the solution in the last column.
 * 
 * @param M     Output matrix, solution.
 * @param A     Input matrix to solve. Linear variables.
 * @param b     Constant vector for the linear equations.
 * @param w     Maximal rank to go to.
 * @return 
 */
long gaussPh4r05(NTL::mat_GF2& M, const NTL::mat_GF2& A, const NTL::vec_GF2 & b, long w);

#endif	/* BASE_H */
