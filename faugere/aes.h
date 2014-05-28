/* 
 * File:   aes.h
 * Author: ph4r05
 *
 * Created on May 27, 2014, 10:19 AM
 */

#ifndef AES_H
#define	AES_H

#include "../base.h"
#define AES_ROUNDS(keysize) ((keysize)==128 ? 10 : ((keysize)==192 ? 12 : 14));
void KeyExpansion(const uchar key[], uint w[], int keysize);
void aes_encrypt(const uchar in[], uchar out[], const uint key[], int keysize, int rounds);
void aes_decrypt(const uchar in[], uchar out[], const uint key[], int keysize, int rounds);

#endif	/* AES_H */

