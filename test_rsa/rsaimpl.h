#ifndef _RSA_IMPLEMENTATION_
#define _RSA_IMPLEMENTATION_

//RSA implementation
#include "Vlong_Gp.h"

//pV = pVx^pVe mod pBr->M
void br_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTBarrettReduction pBr);

//pV = pVx^pVe mod pQr->M
void qr_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTQuisquaterReduction pQr);


//FIPS PUB 186-4 C.3.1 Miller-Rabin Probabilistic Primality Test
int is_probable_prime_br(const PTBarrettReduction pBw);
int is_probable_prime_qr(const PTQuisquaterReduction pQw);
int is_probable_prime_me(const PTMontExp pMe);

//Miller-Rabin Probabilistic Primality Test for untrusted prime candidates
//rounds = [1....5] Error probability = [2^-20....2^-100]
int is_probable_prime(PTVlong pVprime_candidate, int rounds);

//if return 1, return init pVprime
int generate_prime(PTVlong pVprime, const unsigned int bit_len);

//------------------------------------------------------------------------------
typedef 
struct TRSAPrivateQuintuple{ //TRSAPrivateQuintuple
  TVlong P,Q;             //
  TVlong dP,dQ;           // 
  TVlong Qinv;            // 
} TRSAPrivateQuintuple, *PTRSAPrivateQuintuple;

typedef 
struct TRSAPublicKey{ //TRSAPublicKey
  TVlong e;               //
  TVlong N;               // 
} TRSAPublicKey, *PTRSAPublicKey;

void rsapq_init(PTRSAPrivateQuintuple pPQ);
void rsapq_clear(PTRSAPrivateQuintuple pPQ);
void rsapq_delete(PTRSAPrivateQuintuple pPQ);

int rsapq_calculate(PTRSAPrivateQuintuple pPQ, const PTVlong pP, const PTVlong pQ, const PTVlong p_e);

void rsapk_init(PTRSAPublicKey pPK);
void rsapk_clear(PTRSAPublicKey pPK);
void rsapk_delete(PTRSAPublicKey pPK);
//------------------------------------------------------------------------------

//PKCS#1: RSA Cryptography standart. RSAEP: c = m^e mod N
int rsapq_ep(PTVlong pC, const PTVlong pM, const PTRSAPublicKey pPK);
//PKCS#1: RSA Cryptography standart. RSADP: m = c^d mod N
int rsapq_dp(PTVlong pM, const PTVlong pC, const PTRSAPrivateQuintuple pPQ);

//PKCS#1: RSA Cryptography standart. RSAVP1: m = s^e mod N
int rsapq_vp1(PTVlong pM, const PTVlong pS, const PTRSAPublicKey pPK);
//PKCS#1: RSA Cryptography standart. RSASP1: s = m^d mod N
int rsapq_sp1(PTVlong pS, const PTVlong pM, const PTRSAPrivateQuintuple pPQ);

//generate RSA key-pair
//input: 
// - bit_len - length N in bits (512,1024,2048,4096)
// - p_e, if p_e is null, then e is 010001(2^16+1)
//pPK and pPQ must be init before rsa_generate_kp call
//output: pPK(e & N) and pPQ(P & Q only!)
int rsa_generate_kp(unsigned int bit_len, const PTVlong p_e, PTRSAPublicKey pPK, PTRSAPrivateQuintuple pPQ);
//------------------------------------------------------------------------------


#endif //_RSA_IMPLEMENTATION_