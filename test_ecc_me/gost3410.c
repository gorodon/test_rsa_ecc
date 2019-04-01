#include <string.h>
#include <stdlib.h>

#include "gost3410.h"
#include "ecc_Gp.h"
//
//Comments: we use Big-Endian only for simple realization
//
PTBaseEccGp cryptoEccBase256Init(){

	//id-GostR3410-2001-CryptoPro-A-ParamSet

	TPointEcc ecc_G;
	TVlong t;

	PTBaseEccGp pEccBase = NULL;

	const unsigned char str_p[] = {
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0x97};
	const unsigned char str_q[] = {
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0x6C,0x61,0x10,0x70,0x99,0x5A,0xD1,0x00,0x45,0x84,0x1B,0x09,0xB7,0x61,0xB8,0x93};
	const unsigned char str_a[] = {//-3
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0x94};
	const unsigned char str_b[] = {0xa6};
	const unsigned char str_Gx[]= {0x01};
	const unsigned char str_Gy[]= {
		0x8D,0x91,0xE4,0x71,0xE0,0x98,0x9C,0xDA,0x27,0xDF,0x50,0x5A,0x45,0x3F,0x2B,0x76,
		0x35,0x29,0x4F,0x2D,0xDF,0x23,0xE3,0xB1,0x22,0xAC,0xC9,0x9C,0x9E,0x9F,0x1E,0x14};

	//int iRet = 0;
	pEccBase = (PTBaseEccGp)malloc(sizeof(TBaseEccGp));
	//
	if(pEccBase){
		vl_init(&t);
		//
		pEccBase->iModSize = 32;//256 bits
		//
		vl_set_uchar_BE(&t,str_p,sizeof(str_p));
		me_init(&pEccBase->mMe_p,&t);
		//
		vl_set_uchar_BE(&t,str_q,sizeof(str_q));
		me_init(&pEccBase->mMe_q,&t);
		vl_delete(&t);
		//
		vl_init(&pEccBase->va_r);
		vl_init(&pEccBase->vb);
		vl_set_uchar_BE(&pEccBase->va_r,str_a,sizeof(str_a));
		me_mul(&pEccBase->va_r,&pEccBase->mMe_p.r2,&pEccBase->mMe_p);
		vl_set_uchar_BE(&pEccBase->vb,str_b,sizeof(str_b));
		//
		fEccPoint_init(&ecc_G);
		fEccJacPoint_init(&pEccBase->mG);
		vl_set_uchar_BE(&ecc_G.vx,str_Gx,sizeof(str_Gx));
		vl_set_uchar_BE(&ecc_G.vy,str_Gy,sizeof(str_Gy));
		fEcc_P2Jac(&pEccBase->mMe_p,&ecc_G,&pEccBase->mG);
		fEccPoint_delete(&ecc_G);
	}
	//
	return pEccBase;
}

PTBaseEccGp cryptoEccBase512Init(){

	//id-tc26-gost-3410-12-512-paramSetA

	TPointEcc ecc_G;
	TVlong t;

	PTBaseEccGp pEccBase = NULL;

	const unsigned char str_p[] = {
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0xC7};
	const unsigned char str_q[] = {
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0x27,0xE6,0x95,0x32,0xF4,0x8D,0x89,0x11,0x6F,0xF2,0x2B,0x8D,0x4E,0x05,0x60,0x60,
		0x9B,0x4B,0x38,0xAB,0xFA,0xD2,0xB8,0x5D,0xCA,0xCD,0xB1,0x41,0x1F,0x10,0xB2,0x75};
	const unsigned char str_a[] = {//-3
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
		0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0xC4};
	const unsigned char str_b[] = {
		0xE8,0xC2,0x50,0x5D,0xED,0xFC,0x86,0xDD,0xC1,0xBD,0x0B,0x2B,0x66,0x67,0xF1,0xDA,
		0x34,0xB8,0x25,0x74,0x76,0x1C,0xB0,0xE8,0x79,0xBD,0x08,0x1C,0xFD,0x0B,0x62,0x65,
		0xEE,0x3C,0xB0,0x90,0xF3,0x0D,0x27,0x61,0x4C,0xB4,0x57,0x40,0x10,0xDA,0x90,0xDD,
		0x86,0x2E,0xF9,0xD4,0xEB,0xEE,0x47,0x61,0x50,0x31,0x90,0x78,0x5A,0x71,0xC7,0x60};
	const unsigned char str_Gx[]= {0x03};
	const unsigned char str_Gy[]= {
		0x75,0x03,0xCF,0xE8,0x7A,0x83,0x6A,0xE3,0xA6,0x1B,0x88,0x16,0xE2,0x54,0x50,0xE6,
		0xCE,0x5E,0x1C,0x93,0xAC,0xF1,0xAB,0xC1,0x77,0x80,0x64,0xFD,0xCB,0xEF,0xA9,0x21,
		0xDF,0x16,0x26,0xBE,0x4F,0xD0,0x36,0xE9,0x3D,0x75,0xE6,0xA5,0x0E,0x3A,0x41,0xE9,
		0x80,0x28,0xFE,0x5F,0xC2,0x35,0xF5,0xB8,0x89,0xA5,0x89,0xCB,0x52,0x15,0xF2,0xA4};

	//int iRet = 0;
	pEccBase = (PTBaseEccGp)malloc(sizeof(TBaseEccGp));
	//
	if(pEccBase){
		vl_init(&t);
		//
		pEccBase->iModSize = 64;//512 bits
		//
		vl_set_uchar_BE(&t,str_p,sizeof(str_p));
		me_init(&pEccBase->mMe_p,&t);
		//
		vl_set_uchar_BE(&t,str_q,sizeof(str_q));
		me_init(&pEccBase->mMe_q,&t);
		vl_delete(&t);
		//
		vl_init(&pEccBase->va_r);
		vl_init(&pEccBase->vb);
		vl_set_uchar_BE(&pEccBase->va_r,str_a,sizeof(str_a));
		me_mul(&pEccBase->va_r,&pEccBase->mMe_p.r2,&pEccBase->mMe_p);
		vl_set_uchar_BE(&pEccBase->vb,str_b,sizeof(str_b));
		//
		fEccPoint_init(&ecc_G);
		fEccJacPoint_init(&pEccBase->mG);
		vl_set_uchar_BE(&ecc_G.vx,str_Gx,sizeof(str_Gx));
		vl_set_uchar_BE(&ecc_G.vy,str_Gy,sizeof(str_Gy));
		fEcc_P2Jac(&pEccBase->mMe_p,&ecc_G,&pEccBase->mG);
		fEccPoint_delete(&ecc_G);
	}
	//
	return pEccBase;
}

void cryptoEccBaseDelete(PTBaseEccGp *ppEccBase){
	if(ppEccBase && (*ppEccBase)){
		PTBaseEccGp pEccBase = *ppEccBase;
		//
		me_delete(&pEccBase->mMe_p);
		me_delete(&pEccBase->mMe_q);
		//
		vl_delete(&pEccBase->va_r);
		vl_delete(&pEccBase->vb);
		//
		fEccJacPoint_delete(&pEccBase->mG);
		//
		//
		free(pEccBase);
		//
		*ppEccBase = NULL;
	}
}

int cryptoCheckPubKeyValue(PTBaseEccGp pEccBase, const uint8_t *X, const uint8_t *Y)
{
	PTBaseEccGp p_ecc_base;
	TPointEcc ecc_point;
	int result = 0;

	p_ecc_base = pEccBase;

	if(p_ecc_base){
		fEccPoint_init(&ecc_point);
		//
		vl_set_uchar_BE(&ecc_point.vx, X, p_ecc_base->iModSize);
		vl_set_uchar_BE(&ecc_point.vy, Y, p_ecc_base->iModSize);
		//
		result = fEcc_CheckPoint(p_ecc_base, &ecc_point);
		//
		fEccPoint_delete(&ecc_point);
	}
	return result;
}

int cryptoGenerateKeyPair(PTBaseEccGp pEccBase, uint8_t *pubKeyX, uint8_t *pubKeyY, uint8_t *privKey)
{
	uint8_t temp[64];//max Mod size
	int iRet = 0;
	if(!pEccBase)
	{
		return iRet;
	}
	
	if(getRNG(temp,pEccBase->iModSize) != pEccBase->iModSize)
	{
		memset(temp, 0, sizeof(temp));
		return 0;
	}
	//
	iRet = cryptoGenerateKeyPairRnd(pEccBase, temp, pubKeyX, pubKeyY, privKey);
	memset(temp, 0, sizeof(temp));
	//
	return iRet;
}

int cryptoGenerateKeyPairRnd(PTBaseEccGp pEccBase, const uint8_t *rnd, uint8_t *pubKeyX, uint8_t *pubKeyY, uint8_t *privKey)
{
	int i;
	int iRet = 0;
	uint8_t temp[64];//max Mod size

	PTBaseEccGp p_ecc_base;
	TJacPointEcc jacPoint;
	TPointEcc publicKey;
	TVlong privateKey;

	p_ecc_base = pEccBase;

	if (p_ecc_base && rnd && pubKeyX && pubKeyY && privKey)
	{
		fEccJacPoint_init(&jacPoint);
		fEccPoint_init(&publicKey);
		vl_init(&privateKey);
		//
		while (1)
		{
			vl_set_uchar_BE(&privateKey, rnd, p_ecc_base->iModSize);
			while (vl_cf(&privateKey, &p_ecc_base->mMe_q.M) >= 0)//privateKey >= q
				vl_usub(&privateKey, &p_ecc_base->mMe_q.M);      //privateKey-=q
			//
			if (privateKey.value.n == 0)//privateKey == 0
				break;//err
			//p_ecc_base->mG -> jacPoint
			fEcc_Jac2Jac(&jacPoint, &p_ecc_base->mG);
			//publicKey = [privateKey]*G Ep(a,b)
			if (fEcc_MulJacP(&p_ecc_base->mMe_p, &p_ecc_base->va_r, &privateKey, &jacPoint) != 1)
				break;//err
			//jacPoint -> affine publicKey
			if(fEcc_Jac2P(&p_ecc_base->mMe_p, &jacPoint, &publicKey) != 1)
				break;//err
			//
			memset(pubKeyX, 0, p_ecc_base->iModSize);
			memset(pubKeyY, 0, p_ecc_base->iModSize);
			memset(privKey, 0, p_ecc_base->iModSize);
			//
			i = p_ecc_base->iModSize;
			iRet = vl_get_uchar_BE(temp, &i, &publicKey.vx);
			if (1 != iRet)
				break;
			memcpy(pubKeyX + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
			//
			i = p_ecc_base->iModSize;
			iRet = vl_get_uchar_BE(temp, &i, &publicKey.vy);
			if (1 != iRet)
				break;
			memcpy(pubKeyY + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
			//
			i = p_ecc_base->iModSize;
			iRet = vl_get_uchar_BE(temp, &i, &privateKey);
			if (1 != iRet)
				break;
			memcpy(privKey + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
			//
			break;
		}
		//
		fEccJacPoint_delete(&jacPoint);
		fEccPoint_delete(&publicKey);
		vl_delete(&privateKey);
		memset(temp, 0, sizeof(temp));
	}
	//
	return iRet;
}
//
int cryptoVkoGenerate(PTBaseEccGp pEccBase, const uint8_t *PrivKey,  //in
                      const uint8_t *UKM,                            //in
                      uint8_t *pubKeyX, uint8_t *pubKeyY)            //in|out 
{
	uint8_t temp[64];//max Mod size
	int i;
	int iRet = 0;
	//
	PTBaseEccGp p_ecc_base;
	TJacPointEcc jacPoint;
	TPointEcc publicKey;
	TVlong privateKey,vukm,r;
	//
	p_ecc_base = pEccBase;
	//
	if (p_ecc_base && pubKeyX && pubKeyY && PrivKey && UKM)
	{
		fEccJacPoint_init(&jacPoint);
		fEccPoint_init(&publicKey);
		vl_init(&privateKey); vl_init(&vukm); vl_init(&r);
		//
		while (1)
		{
			vl_set_uchar_BE(&privateKey, PrivKey, p_ecc_base->iModSize);
			while (vl_cf(&privateKey, &p_ecc_base->mMe_q.M) >= 0) //privateKey >= q
				vl_usub(&privateKey, &p_ecc_base->mMe_q.M);       //privateKey-=q
			if (privateKey.value.n == 0)
				break;//err
			//
			vl_set_uchar_BE(&vukm, UKM, p_ecc_base->iModSize);
			while (vl_cf(&vukm, &p_ecc_base->mMe_q.M) >= 0) //vukm >= q
				vl_usub(&vukm, &p_ecc_base->mMe_q.M);       //vukm-=q
			if (vukm.value.n == 0)
				break;//err
			//
			vl_set_uchar_BE(&publicKey.vx, pubKeyX, p_ecc_base->iModSize);
			vl_set_uchar_BE(&publicKey.vy, pubKeyY, p_ecc_base->iModSize);
			//
			if (fEcc_CheckPoint(p_ecc_base, &publicKey) != 1)
				break;//err
			//
			fEcc_P2Jac(&p_ecc_base->mMe_p, &publicKey, &jacPoint);//publicKey -> jacPoint
			if (vl_fast_compare(&jacPoint.vX, &p_ecc_base->mG.vX)&&
				vl_fast_compare(&jacPoint.vY, &p_ecc_base->mG.vY)&&
				vl_fast_compare(&jacPoint.vZ, &p_ecc_base->mG.vZ))
				break;//err
			//
			me_modmul(&r, &vukm, &privateKey, &p_ecc_base->mMe_q);//r = ukm*privKey mod q
			//
			if (fEcc_MulJacP(&p_ecc_base->mMe_p, &p_ecc_base->va_r, &r, &jacPoint) != 1)
				break;//err
			//jacPoint -> affine publicKey
			if (fEcc_Jac2P(&p_ecc_base->mMe_p, &jacPoint, &publicKey) != 1)
				break;//err
			//
			memset(pubKeyX, 0, p_ecc_base->iModSize);
			memset(pubKeyY, 0, p_ecc_base->iModSize);
			//
			i = p_ecc_base->iModSize;
			iRet = vl_get_uchar_BE(temp, &i, &publicKey.vx);
			if (1 != iRet)
				break;
			memcpy(pubKeyX + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
			//
			i = p_ecc_base->iModSize;
			iRet = vl_get_uchar_BE(temp, &i, &publicKey.vy);
			if (1 != iRet)
				break;
			memcpy(pubKeyY + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
			//
			break;
		}
		//
		fEccJacPoint_delete(&jacPoint);
		fEccPoint_delete(&publicKey);
		vl_delete(&privateKey); vl_delete(&vukm); vl_delete(&r);
		memset(temp, 0, sizeof(temp));
		//
	}
	//
	return iRet;
}

int cryptoMakeSignature(PTBaseEccGp pEccBase, const uint8_t *hash, //in
                        const uint8_t *privKey,                    //in
                        uint8_t *sig_r, uint8_t *sig_s)            //out
{
	//GOST34.10-2012
	PTBaseEccGp p_ecc_base;
	TJacPointEcc jacPoint;
	TPointEcc pointC;
	TVlong r,s,h,k;
	uint8_t temp[64];//max Mod size
	int i;
	int iRet = 0;
	//
	p_ecc_base = pEccBase;
	//
	if(p_ecc_base && hash && privKey && sig_r && sig_s){
	//
	fEccJacPoint_init(&jacPoint);
	fEccPoint_init(&pointC);
	vl_init(&r);vl_init(&s);vl_init(&h);vl_init(&k);
	//
	vl_set_uchar_BE(&h,hash,p_ecc_base->iModSize);//h = hash(h)
	while(vl_cf(&h,&p_ecc_base->mMe_q.M) >= 0)//h >= q
		vl_usub(&h,&p_ecc_base->mMe_q.M);//h-=q
	if(h.value.n == 0)//h == 0
		vl_set(&h,0,1); //h = 1
	//
	while(1)
	{
		if(getRNG(temp,pEccBase->iModSize) != pEccBase->iModSize)
			break;//err
		//
		vl_set_uchar_BE(&k,temp,pEccBase->iModSize);//k = rand
		while(vl_cf(&k,&p_ecc_base->mMe_q.M) >= 0)//k >= q
			vl_usub(&k,&p_ecc_base->mMe_q.M);//k-=q
		if(k.value.n == 0)//k == 0
			continue;//err - repeat!
		//
		// вычисление точки эллиптической кривой C = [k]*G
		fEcc_Jac2Jac(&jacPoint, &p_ecc_base->mG);//p_ecc_base->mG -> jacPoint
		//
		if (fEcc_MulJacP(&p_ecc_base->mMe_p, &p_ecc_base->va_r, &k, &jacPoint) != 1)
			break;//err
		//jacPoint -> affine pointC
		if(fEcc_Jac2P(&p_ecc_base->mMe_p,&jacPoint,&pointC) != 1)
			break;
		//r = Cx mod q
		vl_copy(&r,&pointC.vx);//r = Cx
		while(vl_cf(&r,&p_ecc_base->mMe_q.M) >= 0)//r >= q
			vl_usub(&r,&p_ecc_base->mMe_q.M);//r-=q
		//
		if(r.value.n == 0)//r == 0
			continue;//err - repeat!
		//
		//s = (k*h + r*d) mod q //d - privKey
		me_modmul(&s,&h,&k,&p_ecc_base->mMe_q);//s = k*h mod q
		//
		vl_set_uchar_BE(&k,privKey,pEccBase->iModSize);//k = privKey
		me_modmul(&h,&r,&k,&p_ecc_base->mMe_q);//h = r*d mod q
		//
		vl_add(&s,&h);//s = (k*h + r*d)
		while(vl_cf(&s,&p_ecc_base->mMe_q.M) >= 0)//s >= q //s = (k*h + r*d) mod q
			vl_usub(&s,&p_ecc_base->mMe_q.M);//s-=q
		//
		if(s.value.n == 0)//s == 0
			continue;//err - repeat!
		//
		memset(sig_r, 0, pEccBase->iModSize);
		memset(sig_s, 0, pEccBase->iModSize);
		//
		i = p_ecc_base->iModSize;
		iRet = vl_get_uchar_BE(temp, &i, &r);
		if (1 != iRet)
			break;
		memcpy(sig_r + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
		//
		i = p_ecc_base->iModSize;
		iRet = vl_get_uchar_BE(temp, &i, &s);
		if (1 != iRet)
			break;
		memcpy(sig_s + (p_ecc_base->iModSize - i), temp, i);//left-side <null> padding
		//
		break;
	}
	//
	fEccJacPoint_delete(&jacPoint);
	fEccPoint_delete(&pointC);
	vl_delete(&r);vl_delete(&s);vl_delete(&h);vl_delete(&k);
	memset(temp, 0, sizeof(temp));
	//
	}
	//
	return iRet;
}

int cryptoVerifySignature(PTBaseEccGp pEccBase, const uint8_t *sig_r, const uint8_t *sig_s,
                          const uint8_t *Hash,
                          const uint8_t *pubKeyX, const uint8_t *pubKeyY)
{
	//GOST34.10-2012
	PTBaseEccGp p_ecc_base;
	TJacPointEcc jacPointP,jacPointQ;
	TPointEcc publicKey;
	TVlong r,s,z1,z2,h;
	//
	int iRet = 0;
	//
	p_ecc_base = pEccBase;
	//
	if(p_ecc_base && sig_r && sig_s && Hash && pubKeyX && pubKeyY){
	//
	fEccJacPoint_init(&jacPointP);fEccJacPoint_init(&jacPointQ);
	fEccPoint_init(&publicKey);
	vl_init(&r);vl_init(&s);vl_init(&z1);vl_init(&z2);vl_init(&h);
	//
	while(1)
	{
		//
		vl_set_uchar_BE(&r,sig_r,p_ecc_base->iModSize);
		vl_set_uchar_BE(&s,sig_s,p_ecc_base->iModSize);
		//r == 0 || r >= q
		if((r.value.n == 0)||(vl_cf(&r,&p_ecc_base->mMe_q.M) >= 0))
			break;//err
		//s == 0 || s >= q
		if((s.value.n == 0)||(vl_cf(&s,&p_ecc_base->mMe_q.M) >= 0))
			break;//err
		//
		vl_set_uchar_BE(&z1,Hash,p_ecc_base->iModSize);//z1 = hash(h)
		//
		if(vl_cf(&z1,&p_ecc_base->mMe_q.M) >= 0)//z1 >= q
			vl_usub(&z1,&p_ecc_base->mMe_q.M);//z1-=q
		if(z1.value.n == 0)//z1 == 0
			vl_set(&z1,0,1); //z1 = 1
		//
		vl_set_uchar_BE(&publicKey.vx,pubKeyX,p_ecc_base->iModSize);
		vl_set_uchar_BE(&publicKey.vy,pubKeyY,p_ecc_base->iModSize);
		if(fEcc_CheckPoint(p_ecc_base, &publicKey) != 1)
			break;//err
		//
		vl_modinv(&h,&z1,&p_ecc_base->mMe_q.M);//h = z1^-1 mod q
		//
		me_modmul(&z1,&h,&s,&p_ecc_base->mMe_q);//z1 = s*h^-1 mod q
		//
		vl_copy(&s,&p_ecc_base->mMe_q.M);//s = q
		vl_sub(&s,&r);//s = -r mod q
		me_modmul(&z2,&h,&s,&p_ecc_base->mMe_q);//z2 = -r*h^-1 mod q
		//
		fEcc_Jac2Jac(&jacPointP, &p_ecc_base->mG);//p_ecc_base->mG -> jacPointP
		//
		fEcc_P2Jac(&p_ecc_base->mMe_p, &publicKey, &jacPointQ);//publicKey -> jacPointQ
		//C = z1*P(G)+z2*Q(pubKey)
		if(fEcc_MulJacPQ(&p_ecc_base->mMe_p,&p_ecc_base->va_r,&z1,&jacPointP,&z2,&jacPointQ) != 1)
			break;//err
		//
		if(fEcc_Jac2P(&p_ecc_base->mMe_p,&jacPointQ,&publicKey) != 1)//C -> publicKey
			break;//err
		//s = Cx mod q
		vl_copy(&s,&publicKey.vx);//s = Cx
		while(vl_cf(&s,&p_ecc_base->mMe_q.M) >= 0)//s >= q
			vl_usub(&s,&p_ecc_base->mMe_q.M);//s-=q
		//
		iRet = vl_fast_compare(&r,&s);//iRet = 1 if r == Cx mod q
		//
		break;
	}
	//
	fEccJacPoint_delete(&jacPointP);fEccJacPoint_delete(&jacPointQ);
	fEccPoint_delete(&publicKey);
	vl_delete(&r);vl_delete(&s);vl_delete(&z1);vl_delete(&z2);vl_delete(&h);
	//
	}
	//
	return iRet;
}

