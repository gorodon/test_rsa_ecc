#ifndef _VLONG_H
#define _VLONG_H

#include <stdint.h>

typedef unsigned int TVLongUnit;

typedef 
struct TVlongValue {
  TVLongUnit *pa;   // array of units
  unsigned int z;   // units allocated
  unsigned int n;   // used units
} TVlongValue, *PTVlongValue;

typedef 
struct TVlong {
  TVlongValue value;
  int negative;
} TVlong, *PTVlong;

// Macros for doing double precision multiply
#define BPU ( 8*sizeof(TVLongUnit) )       // Number of bits in an unit
#define lo(x) ( (x) & (((TVLongUnit)1<<(BPU/2))-1) ) // lower half of unit
#define hi(x) ( (x) >> (BPU/2) )         // upper half
#define lh(x) ( (x) << (BPU/2) )         // make upper half

//--- vlong implementation ---
void vl_init(PTVlong pV);
void vl_inite(PTVlong pV, unsigned int N);//N - size in TVlongValue
void vl_delete(PTVlong pV);
void vl_clear(PTVlong pV);
void vl_normal(PTVlong pV);
void vl_reserve(PTVlong pV, unsigned int x);

void vl_set(PTVlong pV, unsigned int i, TVLongUnit x);
TVLongUnit vl_get(PTVlong pV, unsigned int i);

unsigned int vl_bits(const PTVlong pV);
int vl_test(const PTVlong pV, unsigned int i);
int vl_is_negative(const PTVlong pV);

void vl_shr(PTVlong pV);
void vl_shrx(PTVlong pV, unsigned int x);
void vl_shl(PTVlong pV);
void vl_shlx(PTVlong pV, unsigned int x);

int vl_cf(const PTVlong pVl, const PTVlong pVr);
int vl_fast_compare(const PTVlong pVx, const PTVlong pVy);  //0 if x!=y; 1 if x==y

void vl_copy(PTVlong pVdst, const PTVlong pVsrc);

void vl_add(PTVlong pV, const PTVlong pVx);//pV = pV + pVx
void vl_word_uadd(PTVlong pV, TVLongUnit a);//pV = pV + a  (pV,a > 0)
void vl_inc(PTVlong pV);//pV+=1 (pV must >= 0)
void vl_dec(PTVlong pV);//pV-=1 (pV must >= 1)
void vl_sub(PTVlong pV, const PTVlong pVx);//pV = pV - pVx
void vl_usub(PTVlong pV, const PTVlong pVx);//unsigned pV = pV - pVx (pV must >= pVx)

void vl_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy); //pVres = pVx*pVy
void vl_mul(PTVlong pVx, const PTVlong pVy);//pVx*=pVy
void vl_fast_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy, unsigned int keep );//res = (x*y) mod 2^keep
void vl_fast_add_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy, unsigned int keep );//res += (x*y) mod 2^keep

TVLongUnit vl_mod_word(const PTVlong pVres, TVLongUnit w);//return pVres%w
void vl_dive(PTVlong pVres, PTVlong pVrem, const PTVlong pVx, const PTVlong pVy);//x = res*y + rem; res = x/y rem = x%y

//--- in-out ---
void vl_set_uchar_BE(PTVlong pV, const unsigned char *s, unsigned int len);
int vl_get_uchar_BE(unsigned char *s, unsigned int *plen, const PTVlong pV);
int vl_get_uchar_LE(unsigned char *s, unsigned int *plen, const PTVlong pV);
int vl_set_vlong(PTVlong pV, const TVLongUnit *puiVL, unsigned int uiLength);
int vl_get_vlong(TVLongUnit *puiVL, unsigned int *puiLength, const PTVlong pV);//if puiVL==0, puiLength - size of puiVL

//--- end vlong implementation ---

#endif//#ifndef _VLONG_H