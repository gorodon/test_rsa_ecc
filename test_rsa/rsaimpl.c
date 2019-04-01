#include <stdlib.h>
#include <stdio.h>
#include "rsaimpl.h"
#include "rng.h"

//test
void print_vlong(PTVlong pvl);
//test
inline
void vl_fast_inc(PTVlong pV){
	if(pV){
		TVLongUnit m = vl_get(pV, 0) + 1;
		if(m)
			vl_set(pV,0,m);
		else
			vl_inc(pV);
	}
}
//
inline
void vl_fast_dec(PTVlong pV){
	if(pV){
		TVLongUnit m = vl_get(pV, 0);
		if(m)
			vl_set(pV,0,m - 1);
		else
			vl_dec(pV);
	}
}
// ------------------------------------------------------------------------
//pV = pVx^pVe mod pQr->M
void qr_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTQuisquaterReduction pQr)
{
	unsigned int bits,i;
	TVlong t;
	//
	if(pV && pVx && pVe && pQr){
		vl_init(&t);
		vl_copy(&t,pVx);
		//
		vl_clear(pV);
		vl_set(pV,0,1);
		//
		i = 0;
		bits = vl_bits(pVe);
		while (1)
		{
			if (vl_test(pVe,i)){
				qr_modmul(pV,pV,&t,pQr);
			}
			i += 1;
			if ( i == bits ) break;
			qr_modmul(&t,&t,&t,pQr);
		}
		//
		vl_delete(&t);
		// !!!
		qr_module(pV,pQr);
	}
}
// ------------------------------------------------------------------------
//pV = pVx^pVe mod pBr->M
void br_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTBarrettReduction pBr)
{
	unsigned int bits,i;
	TVlong t;
	//
	if(pV && pVx && pVe && pBr){
		vl_init(&t);
		vl_copy(&t,pVx);
		//
		vl_clear(pV);
		vl_set(pV,0,1);
		//
		i = 0;
		bits = vl_bits(pVe);
		while (1)
		{
			if (vl_test(pVe,i)){
				br_modmul(pV,pV,&t,pBr);//me_mul(pV,&t,pMe);
			}
			i += 1;
			if ( i == bits ) break;
			br_modmul(&t,&t,&t,pBr);//me_mul(&t,&t,pMe);
		}
		//
		vl_delete(&t);
	}
}
// ------------------------------------------------------------------------
static int fpseudornd_init(int *pseed);
static int fpseudornd(int *pseed, uint8_t rng[32]);
// ------------------------------------------------------------------------
//generate random Vlong pV < pVmax (may use pseudorandom)
static int generate_random_Vlong(PTVlong pV, const PTVlong pVmax, int *prndseed)
{
	int res = 0,i;
	TVLongUnit temp[32/sizeof(TVLongUnit)];
	//
	vl_clear(pV);
	vl_reserve(pV,pVmax->value.n + 1);
	while(1){
		if(0 != fpseudornd(prndseed,(uint8_t*)temp))
			break;
		for(i = 0; i < sizeof(temp)/sizeof(TVLongUnit);++i){
			vl_shlx(pV,BPU);
			vl_set(pV, 0, temp[i]);
			if(pV->value.n >= pVmax->value.n){
				res = 1;
				break;
			}
		}
		if(res) break;
	}
	if(res){
		while(vl_cf(pV,pVmax) > 0)//pV > pVmax
			vl_shr(pV);//pV = pV >> 1
	}
	//
	memset(temp, 0 , sizeof(temp));
	//
	if((pV->value.n == 0) || ((pV->value.n == 1)&&(pV->value.pa[0] == (TVLongUnit)1)))
		res = 0;
	//
	return res;
}
//
//generate random Vlong - must use true random!
static int generate_random_Vlong_byLen(PTVlong pV, const unsigned int bytes_len){
	int res = 0,i;
	TVLongUnit temp[32/sizeof(TVLongUnit)];
	const unsigned int unit_len = bytes_len/sizeof(TVLongUnit);
	//
	vl_clear(pV);
	vl_reserve(pV,unit_len + 1);
	while(1){
		if(NO_RNG_ERROR != rngGet((uint8_t*)temp))
			break;
		for(i = 0; i < sizeof(temp)/sizeof(TVLongUnit);++i){
			vl_shlx(pV,BPU);
			vl_set(pV, 0, temp[i]);
			if(pV->value.n >= unit_len){
				res = 1;
				break;
			}
		}
		if(res) break;
	}
	//
	memset(temp, 0 , sizeof(temp));
	//
	if((pV->value.n == 0) || ((pV->value.n == 1)&&(pV->value.pa[0] == (TVLongUnit)1)))
		res = 0;
	//
	return res;
}
// ------------------------------------------------------------------------
//FIPS PUB 186-4 C.3.1 Miller-Rabin Probabilistic Primality Test
int is_probable_prime_me(const PTMontExp pMe)
{
	//Set w = the integer to be tested, w = 1 + m*2^a, where m is odd and 2^a is the largest
	//power of 2 dividing (w - 1)
	int result = 0;
	int rndseed;
	TVlong m,b,z,w_1;
	unsigned int a = 0, i, j, ucontinue, rounds;
	//
	if(!vl_test(&pMe->M,0))//w-четное
		return result;
	if(0 != fpseudornd_init(&rndseed))
		return result;
	//
	vl_inite(&m,pMe->M.value.n+1);vl_copy(&m,&pMe->M);//m = w
	vl_inite(&b,pMe->M.value.n+1);
	vl_inite(&z,pMe->M.value.n+1);
	vl_inite(&w_1,pMe->M.value.n+1);
	//
	vl_fast_dec(&m);//m = w - 1
	vl_copy(&w_1,&m);
	//
	while(!vl_test(&m,0)){//m-четное
		++a;
		vl_shr(&m);
	}//m = (w - 1) / 2^a and a > 0
	//
	//number of rounds writes in FIPS PUB 186-4 Tables C2,C3
	rounds = 27;//256 <= length <= 512 bits
	if(vl_bits(&pMe->M) >= 512)
		rounds = 9 + 1;//(may be 7)for length 512 bits 
	if(vl_bits(&pMe->M) >= 1024)
		rounds = 4 + 1;//4 for length 1024 bits
	//
	for(i = 0; i < rounds; ++i){
		//Generate a random vlong b in the range 1 < b < w - 1
		if(0 == generate_random_Vlong(&b,&w_1,&rndseed))
			break;
		//Set z = b^m mod w.
		me_modexp(&z,&b,&m,pMe);
		//If z == 1, or if z == w - 1, continue
		if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)) continue;
		if(vl_cf(&z,&w_1) == 0) continue;
		//Step 4.5
		ucontinue = 0;
		for(j = 1; j < a; ++j){
			me_modmul(&z,&z,&z,pMe);//z = z^2 mod w
			if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)){
				break;//w is COMPOSITE
			}
			if(vl_cf(&z,&w_1) == 0){//4.5.2 z == w-1
				ucontinue = 1;
				break;//continue test
			}
		}
		//
		if(ucontinue)
			continue;
		//
		break;//w is COMPOSITE
	}
	//
	if(rounds == i) result = 1;
		else	printf(" \n is_probable_prime_me(i=%d) ", i);//test
	//
	vl_delete(&m);
	vl_delete(&b);
	vl_delete(&z);
	vl_delete(&w_1);
	//
	rndseed = 0;
	//
	return result;
}

int is_probable_prime_br(const PTBarrettReduction pBw)
{
	//Set w = the integer to be tested, w = 1 + m*2^a, where m is odd and 2^a is the largest
	//power of 2 dividing (w - 1)
	int result = 0;
	int rndseed;
	TVlong m,b,z,w_1;
	unsigned int a = 0, i, j, ucontinue, rounds;
	//
	if(!vl_test(&pBw->M,0))//w-четное
		return result;
	if(0 != fpseudornd_init(&rndseed))
		return result;
	//
	vl_inite(&m,pBw->M.value.n+1);vl_copy(&m,&pBw->M);//m = w
	vl_inite(&b,pBw->M.value.n+1);
	vl_inite(&z,pBw->M.value.n+1);
	vl_inite(&w_1,pBw->M.value.n+1);
	//
	vl_fast_dec(&m);//m = w - 1
	vl_copy(&w_1,&m);
	//
	while(!vl_test(&m,0)){//m-четное
		++a;
		vl_shr(&m);
	}//m = (w - 1) / 2^a and a > 0
	//
	//number of rounds writes in FIPS PUB 186-4 Tables C2,C3
	rounds = 27;//256 <= length <= 512 bits
	if(vl_bits(&pBw->M) >= 512)
		rounds = 9 + 1;//(may be 7)for length 512 bits 
	if(vl_bits(&pBw->M) >= 1024)//for length 1024 bits
		rounds = 4 + 1;//
	//
	for(i = 0; i < rounds; ++i){
		//Generate a random vlong b in the range 1 < b < w - 1
		if(0 == generate_random_Vlong(&b,&w_1,&rndseed))
			break;
		//Set z = b^m mod w.
		br_modexp(&z,&b,&m,pBw);
		//If z == 1, or if z == w - 1, continue
		if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)) continue;
		if(vl_cf(&z,&w_1) == 0) continue;
		//Step 4.5
		ucontinue = 0;
		for(j = 1; j < a; ++j){
			br_modmul(&z,&z,&z,pBw);//z = z^2 mod w
			if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)){
				break;//w is COMPOSITE
			}
			if(vl_cf(&z,&w_1) == 0){//4.5.2 z == w-1
				ucontinue = 1;
				break;//continue test
			}
		}
		//
		if(ucontinue)
			continue;
		//
		break;//w is COMPOSITE
	}
	//
	if(rounds == i) result = 1;
		else	printf(" \n is_probable_prime_br(i=%d) ", i);//test
	//
	vl_delete(&m);
	vl_delete(&b);
	vl_delete(&z);
	vl_delete(&w_1);
	//
	rndseed = 0;
	//
	return result;
}
/*
int is_probable_prime_qr(const PTQuisquaterReduction pQw)
{
	//Set w = the integer to be tested, w = 1 + m*2^a, where m is odd and 2^a is the largest
	//power of 2 dividing (w - 1)
	int result = 0;
	int rndseed;
	TVlong m,b,z,w_1;
	unsigned int a = 0, i, j, ucontinue, rounds;
	//
	if(!vl_test(&pQw->M,0))//w-четное
		return result;
	if(0 != fpseudornd_init(&rndseed))
		return result;
	//
	vl_inite(&m,pQw->M.value.n+1);vl_copy(&m,&pQw->M);//m = w
	vl_inite(&b,pQw->M.value.n+1);
	vl_inite(&z,pQw->M.value.n+1);
	vl_inite(&w_1,pQw->M.value.n+1);
	//
	vl_fast_dec(&m);//m = w - 1
	vl_copy(&w_1,&m);
	//
	while(!vl_test(&m,0)){//m-четное
		++a;
		vl_shr(&m);
	}//m = (w - 1) / 2^a and a > 0
	//
	rounds = 9 + 1;//(may be 7)for length 512 bits //FIPS PUB 186-4 Tables C2,C3
	if(vl_bits(&pQw->M) >= 1024)//for length 1024 bits
		rounds = 4 + 1;//
	//
	for(i = 0; i < rounds; ++i){
		//Generate a random vlong b in the range 1 < b < w - 1
		if(0 == generate_random_Vlong(&b,&w_1,&rndseed))
			break;
		//Set z = b^m mod w.
		qr_modexp(&z,&b,&m,pQw);
		//If z == 1, or if z == w - 1, continue
		if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)) continue;
		if(vl_cf(&z,&w_1) == 0) continue;
		//Step 4.5
		ucontinue = 0;
		for(j = 1; j < a; ++j){
			qr_modmul(&z,&z,&z,pQw);qr_module(&z,pQw);//z = z^2 mod w
			if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)1)){
				break;//w is COMPOSITE
			}
			if(vl_cf(&z,&w_1) == 0){//4.5.2 z == w-1
				ucontinue = 1;
				break;//continue test
			}
		}
		//
		if(ucontinue)
			continue;
		//
		break;//w is COMPOSITE
	}
	//
	if(rounds == i) result = 1;
		else	printf(" \n is_probable_prime_qr(i=%d) ", i);//test
	//
	vl_delete(&m);
	vl_delete(&b);
	vl_delete(&z);
	vl_delete(&w_1);
	//
	rndseed = 0;
	//
	return result;
}
//*/

//Miller-Rabin Probabilistic Primality Test for untrusted prime candidates
int is_probable_prime(PTVlong pVprime_candidate, int rounds)
{
	TMontExp me;
	int res = 0;
	if(pVprime_candidate && (rounds > 0)){
		me_init(&me,pVprime_candidate);
		for(rounds; rounds > 0; --rounds){
			if(1 != (res = is_probable_prime_me(&me))){
				break;
			}
		}
		me_delete(&me);
	}
	return res;
}

//Fermat test a^(w-1) = 1 mod w; a^w = a mod w
/*
static
int prime_candidate_me(const PTMontExp pMe)
{
	int res = 0;
	TVlong z;
	vl_inite(&z,pMe->M.value.n+1);
	vl_set(&z,0,(TVLongUnit)3);
	me_modexp(&z,&z,&pMe->M,pMe);
	if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)3))// z == 3
		res = 1;
	vl_delete(&z);
	return res;
}
//*/
static
int prime_candidate_me(const PTMontExp pMe)
{
	int res = 0;
	TVlong z;
	vl_inite(&z,pMe->M.value.n+1);
	me_2modexp(&z,&pMe->M,pMe);//z = 2^M mod M
	if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)2))// z == 2
		res = 1;
	vl_delete(&z);
	return res;
}


static
int prime_candidate_br(const PTBarrettReduction pBw)
{
	int res = 0;
	TVlong z;
	vl_inite(&z,pBw->M.value.n+1);
	vl_set(&z,0,(TVLongUnit)3);
	br_modexp(&z,&z,&pBw->M,pBw);
	if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)3))// z == 3
		res = 1;
	vl_delete(&z);
	return res;
}
static
int prime_candidate_qr(const PTQuisquaterReduction pQw)
{
	int res = 0;
	TVlong z;
	vl_inite(&z,pQw->M.value.n+1);
	vl_set(&z,0,(TVLongUnit)3);
	qr_modexp(&z,&z,&pQw->M,pQw);
	if((z.value.n == 1)&&(z.value.pa[0] == (TVLongUnit)3))// z == 3
		res = 1;
	vl_delete(&z);
	return res;
}

//#include <windows.h>//for tests only
extern int GetTickCountMy();//time in ms
//*
//if return 1, init P_prime
int generate_prime(PTVlong pVprime, const unsigned int bit_len){
	TMontExp me;
	unsigned int i,j;
	int res = 0;
	// const uint16_t primelist[80] = {
	// 3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,
	// 59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,
	// 137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
	// 227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
	// 313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419
	// };//primelist
	const uint16_t primelist[] = {
	3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
	101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
	191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
	281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
	389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
	491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
	607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
	719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
	829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
	953,967,971,977,983,991,997,
	1009,1013,1019,1021,1031,1033,1039,1049,
	1051,1061,1063,1069,1087,1091,1093,1097,1103,1109,1117,1123,1129,1151,
	1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229
	};
	//
	int iStart;//test
	//
	uint16_t *prime_modules;//[sizeof(primelist)/sizeof(primelist[0])];
	//
	if(!pVprime || (bit_len < 128) || (0 != bit_len % 32))
		return res;
	//
	prime_modules = (uint16_t*)malloc(sizeof(primelist));
	//
	if(1 == generate_random_Vlong_byLen(pVprime,bit_len >> 3)){
		//
		iStart = GetTickCountMy();//test 1 step
		//set first m to odd
		pVprime->value.pa[0] |= (TVLongUnit)0x01;
		//set m-s (bit_len-1) & (bit_len-2) bits in 1
		i = bit_len - 2;
		vl_set(pVprime,i/BPU,vl_get(pVprime,i/BPU)|((TVLongUnit)3 << (i%BPU)));
		//
		//compute primelist modules
		for (i=0;i<sizeof(primelist)/sizeof(primelist[0]);i+=1)
		{
			unsigned int p = primelist[i];
			unsigned int r = (unsigned int)vl_mod_word(pVprime,p);
			if(r && (r < p)) r = p - r;
			prime_modules[i] = r;
		}
		//test start
		printf("\ngenerate_primes: 1 step = %d ms", GetTickCountMy() - iStart);
		unsigned num = 0;
		iStart = GetTickCountMy();
		//test end
		//find prime...
		for(i = 0; i < 1024*1024; i+=2, vl_word_uadd(pVprime,2))
		{
			unsigned int delta;
			//
			for(delta = 0; delta < 1024*1024; delta+=2){
				int composite = 0;
				unsigned int mod;
				//
				for (j=0;j<sizeof(primelist)/sizeof(primelist[0]);j+=1){
					
					mod = (i+delta)%(unsigned int)primelist[j];
					//check that gcd((pVprime+delta-1,primelist) == 1
					if(mod - prime_modules[j] <= 1){
						composite = 1;
						break;
					}
					// if(prime_modules[j] == mod){
						// composite = 1;
						// break;
					// }
				}
				//
				if(!composite)//found prime candidate with delta shift
					break;
			}
			//found prime candidate with delta shift
			vl_word_uadd(pVprime,delta);
			i+=delta;
			//
			num+=1;//test str
			me_init(&me,pVprime);
			if(1 == prime_candidate_me(&me))
			{
				if(1 == (res = is_probable_prime_me(&me))){
					me_delete(&me);
					break;
				}
			}
			me_delete(&me);
		}
		printf("\ngenerate_primes: 2 step = %d ms", GetTickCountMy() - iStart);//test
		printf("\ngenerate_primes: num test = %d, difference = %d", num, i);//test
	}
	//
	free(prime_modules);
	//
	return res;
}
//*/

/*
//--- bitmap manage functions
inline int bs_get_bit(TVLongUnit* pmap,unsigned int x){
	return (pmap[x/BPU] >> (x%BPU))&0x01;
}

inline void bs_set_bit(TVLongUnit* pmap,unsigned int x,int mbit){
	if(mbit)
		pmap[x/BPU] |=  (TVLongUnit)1 << (x%BPU);
	else
		pmap[x/BPU] &=  ~((TVLongUnit)1 << (x%BPU));
}
//--- bitmap functions
//*/

//------------------------------------------------------------------------------
void rsapq_init(PTRSAPrivateQuintuple pPQ){
	if(pPQ){
		vl_init(&pPQ->P);
		vl_init(&pPQ->Q);
		vl_init(&pPQ->dP);
		vl_init(&pPQ->dQ);
		vl_init(&pPQ->Qinv);
	}
};

void rsapq_clear(PTRSAPrivateQuintuple pPQ){
	if(pPQ){
		vl_clear(&pPQ->P);
		vl_clear(&pPQ->Q);
		vl_clear(&pPQ->dP);
		vl_clear(&pPQ->dQ);
		vl_clear(&pPQ->Qinv);
	}
};

void rsapq_delete(PTRSAPrivateQuintuple pPQ){
	if(pPQ){
		vl_delete(&pPQ->P);
		vl_delete(&pPQ->Q);
		vl_delete(&pPQ->dP);
		vl_delete(&pPQ->dQ);
		vl_delete(&pPQ->Qinv);
	}
};

int rsapq_calculate(PTRSAPrivateQuintuple pPQ, const PTVlong pP, const PTVlong pQ, const PTVlong p_e)
{
	int res = 0;
	if(pPQ && pP && pQ && p_e){
		//
		//if(&pPQ->P != pP)
			vl_copy(&pPQ->P,pP);
		//if(&pPQ->Q != pQ)
			vl_copy(&pPQ->Q,pQ);
		//
		vl_dec(&pPQ->P);//P = P - 1
		vl_modinv_all(&pPQ->dP,p_e,&pPQ->P);//dP = e^(-1) mod (P-1)
		vl_inc(&pPQ->P);//restore P

		vl_dec(&pPQ->Q);//Q = Q - 1
		vl_modinv_all(&pPQ->dQ,p_e,&pPQ->Q);//dQ = e^(-1) mod (Q-1)
		vl_inc(&pPQ->Q);//restore Q

		if(vl_cf(&pPQ->P,&pPQ->Q) < 0){//P < Q
			vl_usub(&pPQ->Q,&pPQ->P);// Q = Q mod P = Q - P;
			vl_modinv(&pPQ->Qinv,&pPQ->Q,&pPQ->P);//Qinv = Q^(-1) mod P
			vl_add(&pPQ->Q,&pPQ->P);//restore Q
		}else{
			vl_modinv(&pPQ->Qinv,&pPQ->Q,&pPQ->P);//Qinv = Q^(-1) mod P
		}
		//
		res = 1;
	}
	return res;
};

void rsapk_init(PTRSAPublicKey pPK)
{
	if(pPK){
		vl_init(&pPK->e);
		vl_init(&pPK->N);
	}
};

void rsapk_clear(PTRSAPublicKey pPK)
{
	if(pPK){
		vl_clear(&pPK->e);
		vl_clear(&pPK->N);
	}
};

void rsapk_delete(PTRSAPublicKey pPK)
{
	if(pPK){
		vl_delete(&pPK->e);
		vl_delete(&pPK->N);
	}
};

//*
//PKCS#1: RSA Cryptography standart. RSADP:
int rsapq_dp(PTVlong pM, const PTVlong pC, const PTRSAPrivateQuintuple pPQ)
{
	int res = 0;
	TMontExp me;
	TVlong m1,m2;
	//
	if(pM && pC && pPQ){
		vl_init(&m1);
		vl_init(&m2);
		//
		vl_mule(&m1,&pPQ->Q,&pPQ->P);
		//C < P*Q
		if(vl_cf(pC,&m1) < 0)//if(vl_bits(pC) <= vl_bits(&pPQ->Q)+vl_bits(&pPQ->P))
		{
			vl_clear(pM);
			vl_clear(&m1);
			//hereinafter: c = pM
			//m2 = c^dQ mod Q
			me_init(&me,&pPQ->Q);//me - mod Q
			vl_copy(pM,pC);
			me_module(pM,&me);//c = pC mod Q
			me_modexp(&m2,pM,&pPQ->dQ,&me);
			me_delete(&me);
			//m1 = c^dP mod P
			me_init(&me,&pPQ->P);//me - mod P
			vl_copy(pM,pC);
			me_module(pM,&me);//c = pC mod P
			me_modexp(&m1,pM,&pPQ->dP,&me);
			// h = (m1-m2)*Qinv mod P
			while(vl_cf(&m1,&m2) < 0){//m1 < m2
				vl_add(&m1,&pPQ->P); //m1 += P
			}
			//hereinafter: h = pM
			vl_copy(pM,&m1);
			vl_usub(pM,&m2);//h = (m1-m2) mod P
			me_modmul(pM,pM,&pPQ->Qinv,&me);
			me_delete(&me);
			//m = Q*h + m2
			//vl_mul(pM,&pPQ->Q);//h = h*Q
			//vl_add(pM,&m2);
			vl_mule(&m1,pM,&pPQ->Q);//m1 = h*Q
			vl_add(&m1,&m2);//m1 = h*Q + m2
			vl_copy(pM,&m1);
			//
			res = 1;
		}
		//
		vl_delete(&m1);
		vl_delete(&m2);
	}
	//
	return res;
};
//*/
/*
//PKCS#1: RSA Cryptography standart. RSADP:
int rsapq_dp(PTVlong pM, const PTVlong pC, const PTRSAPrivateQuintuple pPQ)
{
	int res = 0;
	TBarrettReduction br;
	TVlong m1,m2;
	//
	if(pM && pC && pPQ){
		vl_init(&m1);
		vl_init(&m2);
		//
		vl_mule(&m1,&pPQ->Q,&pPQ->P);
		//C < P*Q
		if(vl_cf(pC,&m1) < 0)
		{
			//printf("\n-rsapq_dp: pC = ");print_vlong(pC);

			vl_clear(&m1);
			//hereinafter: c = pM
			//m2 = c^dQ mod Q
			br_init(&br,&pPQ->Q);//br - mod Q
			vl_clear(pM);vl_copy(pM,pC);br_module(pM,&br);//c = pC mod Q

			//printf("\n-rsapq_dp: pC mod Q = ");print_vlong(pM);

			//if(pM->value.n == 0)
			//	printf("\n-rsapq_dp: pC mod Q = 0!!!");
			br_modexp(&m2,pM,&pPQ->dQ,&br);
			br_delete(&br);

			//printf("\n-rsapq_dp: c^dQ mod Q = ");print_vlong(&m2);

			//m1 = c^dP mod P
			br_init(&br,&pPQ->P);//br - mod P
			vl_clear(pM);vl_copy(pM,pC);br_module(pM,&br);//c = pC mod P

			//printf("\n-rsapq_dp: pC mod P = ");print_vlong(pM);

			//if(pM->value.n == 0)
			//	printf("\n-rsapq_dp: pC mod P = 0!!!");
			br_modexp(&m1,pM,&pPQ->dP,&br);

			//printf("\n-rsapq_dp: c^dP mod P = ");print_vlong(&m1);

			// h = (m1-m2)*Qinv mod P
			while(vl_cf(&m1,&m2) < 0){//m1 < m2
				vl_add(&m1,&pPQ->P); //m1 += P
			}
			//hereinafter: h = pM
			vl_clear(pM);vl_copy(pM,&m1);//
			vl_usub(pM,&m2);//h = (m1-m2) mod P
			br_modmul(pM,pM,&pPQ->Qinv,&br);
			br_delete(&br);
			//m = Q*h + m2
			//vl_mul(pM,&pPQ->Q);//h = h*Q
			//vl_add(pM,&m2);
			vl_mule(&m1,pM,&pPQ->Q);//m1 = h*Q
			vl_add(&m1,&m2);//m1 = h*Q + m2
			vl_copy(pM,&m1);
			//
			res = 1;
		}
		//
		vl_delete(&m1);
		vl_delete(&m2);
	}
	//
	return res;
};
//*/
/*
//PKCS#1: RSA Cryptography standart. RSADP:
int rsapq_dp(PTVlong pM, const PTVlong pC, const PTRSAPrivateQuintuple pPQ)
{
	int res = 0;
	TQuisquaterReduction qr;
	TVlong m1,m2;
	//
	if(pM && pC && pPQ){
		vl_init(&m1);
		vl_init(&m2);
		//
		vl_mule(&m1,&pPQ->Q,&pPQ->P);
		//C < P*Q
		if(vl_cf(pC,&m1) < 0)
		{
			//printf("\n-rsapq_dp: pC = ");print_vlong(pC);
			//
			vl_clear(&m1);
			//hereinafter: c = pM
			//m2 = c^dQ mod Q
			qr_init(&qr,&pPQ->Q);//qr - mod Q
			vl_clear(pM);
			//vl_dive(0,pM,pC,&pPQ->Q);//c = pC mod Q
			if(pC->value.n < pPQ->Q.value.n){
				vl_copy(pM,pC);
			}else{
				vl_set(pM,0,1);qr_modmul(pM,pC,pM,&qr);qr_module(pM,&qr);//c = pC mod Q
			}
			qr_modexp(&m2,pM,&pPQ->dQ,&qr);
			qr_delete(&qr);
			//printf("\n-rsapq_dp: c^dQ mod Q = ");print_vlong(&m2);

			//m1 = c^dP mod P
			qr_init(&qr,&pPQ->P);//qr - mod P
			vl_clear(pM);
			//vl_dive(0,pM,pC,&pPQ->P);//c = pC mod P
			if(pC->value.n < pPQ->P.value.n){
				vl_copy(pM,pC);
			}else{
				vl_set(pM,0,1);qr_modmul(pM,pC,pM,&qr);qr_module(pM,&qr);//c = pC mod P
			}
			//printf("\n-rsapq_dp: pC mod P = ");print_vlong(pM);
			qr_modexp(&m1,pM,&pPQ->dP,&qr);
			//printf("\n-rsapq_dp: c^dP mod P = ");print_vlong(&m1);

			// h = (m1-m2)*Qinv mod P
			while(vl_cf(&m1,&m2) < 0){//m1 < m2
				vl_add(&m1,&pPQ->P); //m1 += P
			}
			//hereinafter: h = pM
			vl_clear(pM);vl_copy(pM,&m1);//
			vl_usub(pM,&m2);//h = (m1-m2) mod P
			qr_modmul(pM,pM,&pPQ->Qinv,&qr);qr_module(pM,&qr);
			qr_delete(&qr);
			//m = Q*h + m2
			//vl_mul(pM,&pPQ->Q);//h = h*Q
			//vl_add(pM,&m2);
			vl_mule(&m1,pM,&pPQ->Q);//m1 = h*Q
			vl_add(&m1,&m2);//m1 = h*Q + m2
			vl_copy(pM,&m1);
			//
			res = 1;
		}
		//
		vl_delete(&m1);
		vl_delete(&m2);
	}
	//
	return res;
};
//*/

//PKCS#1: RSA Cryptography standart. RSAEP: c = m^e mod N
static
int rsapq_ep_me(PTVlong pC, const PTVlong pM, const PTRSAPublicKey pPK)
{
	int res = 0;
	TMontExp me;
	if(pC && pM && pPK){
		if(vl_cf(pM,&pPK->N) < 0){//M < N
			me_init(&me,&pPK->N);
			me_modexp(pC,pM,&pPK->e,&me);
			me_delete(&me);
			res = 1;
		}
	}
	return res;
};

//PKCS#1: RSA Cryptography standart. RSAEP: c = m^e mod N
static
int rsapq_ep_qr(PTVlong pC, const PTVlong pM, const PTRSAPublicKey pPK)
{
	int res = 0;
	TQuisquaterReduction qrN;
	if(pC && pM && pPK){
		if(vl_cf(pM,&pPK->N) < 0){//M < N
			qr_init(&qrN,&pPK->N);
			qr_modexp(pC,pM,&pPK->e,&qrN);
			qr_delete(&qrN);
			res = 1;
		}
	}
	return res;
};

/*
int rsapq_ep_br(PTVlong pC, const PTVlong pM, const PTRSAPublicKey pPK)
{
	int res = 0;
	TBarrettReduction brN;
	if(pC && pM && pPK){
		if(vl_cf(pM,&pPK->N) < 0){//M < N
			br_init(&brN,&pPK->N);
			br_modexp(pC,pM,&pPK->e,&brN);
			br_delete(&brN);
			res = 1;
		}
	}
	return res;
};
//*/

//PKCS#1: RSA Cryptography standart. RSAEP: c = m^e mod N
int rsapq_ep(PTVlong pC, const PTVlong pM, const PTRSAPublicKey pPK){
	int res = 0;
	if(pC && pM && pPK){
		if(pPK->e.value.n < 2){
			res = rsapq_ep_qr(pC,pM,pPK);//e is short
		}else{
			res = rsapq_ep_me(pC,pM,pPK);//e is long
		}
	}
	return res;
}

//PKCS#1: RSA Cryptography standart. RSADP: m = c^d mod N, by Montgomery me, for pPK->e=d
int rsapq_dp_me(PTVlong pM, const PTVlong pC, const PTRSAPublicKey pPK)
{
	int res = 0;
	TMontExp me;
	if(pC && pM && pPK){
		if(vl_cf(pC,&pPK->N) < 0){//C < N
			me_init(&me,&pPK->N);
			me_modexp(pM,pC,&pPK->e,&me);
			me_delete(&me);
			res = 1;
		}
	}
	return res;
};

//PKCS#1: RSA Cryptography standart. RSASP1: s = m^d mod N
int rsapq_sp1(PTVlong pS, const PTVlong pM, const PTRSAPrivateQuintuple pPQ)
{
	return rsapq_dp(pS,pM,pPQ);
};

//PKCS#1: RSA Cryptography standart. RSAVP1: m = s^e mod N
int rsapq_vp1(PTVlong pM, const PTVlong pS, const PTRSAPublicKey pPK)
{
	return rsapq_ep(pM,pS,pPK);
};

//generate RSA key-pair
// - bit_len - length N in bits (1024,2048)
// - p_e, if p_e is null, then e is 010001(2^16+1)
//pPK and pPQ must be init before rsa_generate_kp call
//output: pPK(e & N) and pPQ(P & Q only!)
int rsa_generate_kp(unsigned int bit_len, const PTVlong p_e, PTRSAPublicKey pPK, PTRSAPrivateQuintuple pPQ)
{
	int res = 0;
	TVlong m;
	vl_init(&m);
	//
	while(pPK && pPQ && (bit_len >= 512)){
		//init e
		if(p_e){
			vl_copy(&pPK->e,p_e);
		}else{
			unsigned char str_e[] = {	0x01, 0x00, 0x01 };
			vl_set_uchar_BE(&pPK->e,str_e,sizeof(str_e));
		}
		//Generate P
		if(1 != generate_prime(&pPQ->P,bit_len >> 1)){
			res = -1;
			break;
		}
		//check GCD(P-1,e) == 1
		vl_dec(&pPQ->P);//P-=1
		if(1 == vl_gcd(&m,&pPQ->P,&pPK->e)){
			if((m.value.n == 1)&&(m.value.pa[0] == (TVLongUnit)1))// m == 1
				vl_inc(&pPQ->P);//restore P
			else{
				//res = -2;
				//break;//error
				continue;
			}
		}else{
			res = -3;
			break;
		}
		//
		//Generate Q
		if(1 != generate_prime(&pPQ->Q,bit_len >> 1)){
			res = -11;
			break;
		}
		//check GCD(Q-1,e) == 1
		vl_dec(&pPQ->Q);//Q-=1
		if(1 == vl_gcd(&m,&pPQ->Q,&pPK->e)){
			if((m.value.n == 1)&&(m.value.pa[0] == (TVLongUnit)1))// m == 1
				vl_inc(&pPQ->Q);//restore Q
			else{
				//res = -12;
				//break;//error
				continue;
			}
		}else{
			res = -13;
			break;
		}
		//
		//calculate N, check bit_len
		vl_mule(&pPK->N,&pPQ->Q,&pPQ->P);
		if(bit_len != vl_bits(&pPK->N)){
			res = -20;
			break;
		}
		//
		res = 1;
		break;
	}
	//
	vl_delete(&m);
	//
	if(1 != res){
		rsapk_clear(pPK);
		rsapq_clear(pPQ);
	}
	//
	return res;
};
// ------------------------------------------------------------------------
//light pseudorandom:
static inline int pseudornd(int *pnext){
  int hi, lo, x;
  // Can't be initialized with 0, so use another value.
  if (*pnext == 0)
    *pnext = 123459876;
  hi = *pnext / 127773;
  lo = *pnext % 127773;
  x = 16807 * lo - 2836 * hi;
  if (x < 0)
    x += 0x7fffffff;
  *pnext = x;
  return (x & 0x7fffffff);
}

static int fpseudornd(int *pseed, uint8_t rng[32]){
	int i,ret;
	ret = 1;
	if(pseed && rng){
		for(i = 0; i<32; i++){
		  rng[i] = (uint8_t)pseudornd(pseed);
		}
		ret = 0;//ok
	}else{
		ret = -1;
	}
	return ret;
}

static int fpseudornd_init(int *pseed){
	int ret = 1;
	int temp[32/sizeof(int)];
	if(pseed){
		if(NO_RNG_ERROR != rngGet((uint8_t*)temp))
			ret = -1;
		else{
			*pseed = temp[0];
			memset(temp, 0 , sizeof(temp));
			ret = 0;//ok
		}
	}else{
		ret = -2;
	}
	return ret;
}
//------------------------------------------------------------------------------
