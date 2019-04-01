#ifndef _VLONG_GP_H
#define _VLONG_GP_H

#include "Vlong.h"

//--- Modular arithmetic implementation ---

int vl_gcd(PTVlong pV, const PTVlong pVx, const PTVlong pVy);//pV = GCD(pVx,pVy)
void vl_modinv_all(PTVlong pV, const PTVlong pVa, const PTVlong pVm);// modular inverse
void vl_modinv(PTVlong pV, const PTVlong pVa, const PTVlong pVm);// modular inverse // m must be a prime(!)

//--- BarrettReduction ---

//объект для выполнения быстрых операций a*b mod p
typedef
struct TBarrettReduction{ // struct for Barrett reduction - fast calc R = N mod M
  unsigned int uiN;       // bits for M
  unsigned int uiMu;      // bits for mu
  TVlong M;               // M - module
  TVlong mu;              // mu = 2^2n / M
  TVlong Q,T,res;         // buffer
} TBarrettReduction, *PTBarrettReduction;

void br_init(PTBarrettReduction pBr, const PTVlong pM);   //init pBr with module pM
void br_init_uchar_BE(PTBarrettReduction pBr,             //fast init with precalculated settings
                      const unsigned char *psM, unsigned int lenM,
                      const unsigned char *psMu, unsigned int lenMu);
void br_init_copy(PTBarrettReduction pBr, const PTBarrettReduction pBrIni);   //init pBr with pBrIni (copy)
void br_delete(PTBarrettReduction pBr);                   //delete pBr
void br_module(PTVlong pV, const PTBarrettReduction pBr); //pV = pV mod pBr->M

//pV = pVx*pVy mod pBr->M
void br_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTBarrettReduction pBr);

//--- end BarrettReduction ---


//--- Montgomery modular exponentiation ---
typedef
struct TMontExp{ // struct for Montgomery modular exponentiation - fast calc  x^e mod M
  unsigned int uiN;       // bits for R = 2^bits(M)
  TVlong M;               // M - module
  TVLongUnit n1short;     // n1s = -M^(-1) mod 2^bpu
  //TVlong n1;              // n1 = R - modinv_all(M,R) = -M^(-1) mod R
  TVlong r2;              // r2 = R*R mod M
  TVlong T,K;             // buffer
} TMontExp, *PTMontExp;

void me_init(PTMontExp pMe, const PTVlong pM);   //init pMe with module pM
void me_delete(PTMontExp pMe);
//pV = pVx^pVe mod pMe->M
void me_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTMontExp pMe);
//pV = 2^pVe mod pMe->M
void me_2modexp(PTVlong pV, const PTVlong pVe, const PTMontExp pMe);
//pV = pVx*pVy mod pBr->M
void me_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTMontExp pMe);
//
void me_module(PTVlong pV, const PTMontExp pMe);//pV = pV mod pMe->M
//pV = pV*R^-1 mod M
void me_mont(PTVlong pV, const PTMontExp pMe);
//pV = pV*pVy*R^-1 mod M
void me_mul(PTVlong pV, const PTVlong pVy, const PTMontExp pMe);
//--- end Montgomery modular exponentiation ---

//--- Quisquater reduction ---
typedef
struct TQuisquaterReduction{ // struct for Quisquater reduction - fast calc R = N mod M
  unsigned int uiN;       // bits for M
  TVlong M;               // M - module
  //TVlong del;             // del = |2^(n+bpu) / M|
  TVlong N;               // N = |2^(n+bpu) / M| * M
  TVlong Nsh;             // Nsh = 2^(n+bpu) - N
  TVlong res;//,q;        // buffer
} TQuisquaterReduction, *PTQuisquaterReduction;

void qr_init(PTQuisquaterReduction pQr, const PTVlong pM);   //init pQr with module pM
void qr_delete(PTQuisquaterReduction pQr);                   //delete pQr

//pV = pVx*pVy mod pQr->N(!) pVx,pVy may len = len(M) + BPU
void qr_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTQuisquaterReduction pQr);
//pV = (pV mod pQr->N) mod pQr->M
void qr_module(PTVlong pV, const PTQuisquaterReduction pQr);

//--- end Quisquater reduction ---

#endif//#ifndef _VLONG_GP_H
