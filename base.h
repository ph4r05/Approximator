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

template<class T>
void dumpHex(std::ostream & c, const T * inp, unsigned int size, bool endl=1) {
    c << std::showbase // show the 0x prefix
      << std::internal // fill between the prefix and the number
      << std::setfill('0'); // fill with 0s
    
    for(unsigned int i = 0; i < size; i++){
        c << std::hex << std::setw(sizeof(T)) << (static_cast<const ULONG>(inp[i])) << " ";
    }
    
    if (endl){
        c << std::endl;
    }
}

inline void dumpUcharHex(std::ostream & c, const uchar* inp, unsigned int size, bool endl=1){ dumpHex(c, inp, size, endl); }
inline void dumpUlongHex(std::ostream & c, const ULONG* inp, unsigned int size, bool endl=1){ dumpHex(c, inp, size, endl); }
void dumpUchar   (std::ostream & c, const uchar * inp, unsigned int size, bool endl=1);

#endif	/* BASE_H */
