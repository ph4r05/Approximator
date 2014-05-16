/* 
 * File:   base.h
 * Author: ph4r05
 *
 * Created on May 16, 2014, 2:06 PM
 */

#ifndef BASE_H
#define	BASE_H
#include <iostream>

typedef unsigned long long ULONG;
typedef unsigned char uchar;
typedef unsigned int  uint;

// Fast ceiling function for integers.
#define OWN_CEIL(x)  (    (((int)(x)) < (x)) ? ((int)(x))+1 : ((int)(x))    )
#define OWN_FLOOR(x) (    (((int)(x)) < (x)) ? ((int)(x))-1 : ((int)(x))    )

void dumpUcharHex(std::ostream & c, uchar* inp, unsigned int size);
void dumpUchar   (std::ostream & c, uchar * inp, unsigned int size);

#endif	/* BASE_H */

