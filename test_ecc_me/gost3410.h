#ifndef _GOST_3410_
#define _GOST_3410_
// signature algorithm - GOST 34.10-2012
// apply 256-bits elliptic curve id-GostR3410-2001-CryptoPro-A-ParamSet - see rfc4357
// apply 512-bits elliptic curve id-tc26-gost-3410-12-512-paramSetA - see rfc7836

// Comments: we used Big-Endian only for simple realization
//
#include <stdint.h>

#include "ecc_Gp.h"

PTBaseEccGp cryptoEccBase256Init();//id-GostR3410-2001-CryptoPro-A-ParamSet
PTBaseEccGp cryptoEccBase512Init();//id-tc26-gost-3410-12-512-paramSetA

void cryptoEccBaseDelete(PTBaseEccGp *ppEccBase);

// 1 - generated
int cryptoGenerateKeyPair(PTBaseEccGp pEccBase, uint8_t *pubKeyX, uint8_t *pubKeyY, uint8_t *privKey);
int cryptoGenerateKeyPairRnd(PTBaseEccGp pEccBase, const uint8_t *rnd, uint8_t *pubKeyX, uint8_t *pubKeyY, uint8_t *privKey);
// 1 - pubkey value is true
int cryptoCheckPubKeyValue(PTBaseEccGp pEccBase, const uint8_t *X, const uint8_t *Y);
//
int cryptoVkoGenerate(PTBaseEccGp pEccBase, const uint8_t *PrivKey,  //in
                      const uint8_t *UKM,                            //in
                      uint8_t *pubKeyX, uint8_t *pubKeyY);           //in|out 
//
int cryptoMakeSignature(PTBaseEccGp pEccBase, const uint8_t *hash, //in
                        const uint8_t *privKey,                    //in
                        uint8_t *sig_r, uint8_t *sig_s);           //out
//
// 1 - signature is verify
int cryptoVerifySignature(PTBaseEccGp pEccBase, const uint8_t *sig_r, const uint8_t *sig_s,
                          const uint8_t *Hash,
                          const uint8_t *pubKeyX, const uint8_t *pubKeyY);
//

//extern, must be released
//return iSize(>0) if successed
//return -1 if error
int getRNG(uint8_t *rnd, int iSize);
//
#endif//#ifndef _GOST_3410_
