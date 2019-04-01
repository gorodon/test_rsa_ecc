#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#include <windows.h>
#endif

#include "Vlong_Gp.h"
#include "ecc_Gp.h"
#include "gost3410.h"

void print_hex_str(const unsigned char *strHex, unsigned int uLenHex, char *strRes);
void print_hex(const unsigned char *strHex, unsigned int uLenHex);
void print_vlong(PTVlong pvl);
//
void test_ecc_raw();
void test_ecc_raw_2();
void test_flow_gost3410_256();
void test_flow_gost3410_512();
void test_cryptoVkoGenerate(PTBaseEccGp pEccBase);
//
int GetTickCountMy();//time in ms
void rngInit();

int main(int argc, char *argv[])
{
  rngInit();
  //Stage 1
  //test_ecc_raw();
  //Stage 2
  //test_ecc_raw_2();
  //
  test_flow_gost3410_256();
  //
  test_flow_gost3410_512();
  //
  return 0;
}

void test_ecc_raw(){
  int iRet;
  TVlong p,a,d;
  TPointEcc eccPointG;
  TJacPointEcc jacPointP,jacPointQ,jacPointR;
  TMontExp mMe_p;
  TVlong pvl_buff[9];
  int i;
//for p
  unsigned char str_p[] = {
    0x80,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x04,0x31};
  unsigned char str_a[] = {0x07};
  unsigned char str_d[] = {
    0x7A,0x92,0x9A,0xDE,0x78,0x9B,0xB9,0xBE,0x10,0xED,0x35,0x9D,0xD3,0x9A,0x72,0xC1,
    0x1B,0x60,0x96,0x1F,0x49,0x39,0x7E,0xEE,0x1D,0x19,0xCE,0x98,0x91,0xEC,0x3B,0x28};
  unsigned char str_Gx[]= {0x02};
  unsigned char str_Gy[]= {
    0x08,0xE2,0xA8,0xA0,0xE6,0x51,0x47,0xD4,0xBD,0x63,0x16,0x03,0x0E,0x16,0xD1,0x9C,
    0x85,0xC9,0x7F,0x0A,0x9C,0xA2,0x67,0x12,0x2B,0x96,0xAB,0xBC,0xEA,0x7E,0x8F,0xC8};
  //
  vl_init(&p);
  vl_init(&a);
  vl_init(&d);
  //
  vl_set_uchar_BE(&p,str_p,sizeof(str_p));
  vl_set_uchar_BE(&a,str_a,sizeof(str_a));
  vl_set_uchar_BE(&d,str_d,sizeof(str_d));

  me_init(&mMe_p,&p); 
  fEccPoint_init(&eccPointG);

  for(i=0; i<9; i++)
    vl_init(&pvl_buff[i]);

  vl_set_uchar_BE(&eccPointG.vx,str_Gx,sizeof(str_Gx));
  vl_set_uchar_BE(&eccPointG.vy,str_Gy,sizeof(str_Gy));
  //
  printf("\np  = ");
  print_vlong(&p);
  printf("\na  = ");
  print_vlong(&a);
  
  printf("\na*R= ");
  me_mul(&a,&mMe_p.r2,&mMe_p);
  print_vlong(&a);

  printf("\nGx = ");
  print_vlong(&eccPointG.vx);
  printf("\nGy = ");
  print_vlong(&eccPointG.vy);
  //
  fEccJacPoint_init(&jacPointP);
  fEccJacPoint_init(&jacPointQ);
  fEccJacPoint_init(&jacPointR);
  //
  fEcc_P2Jac(&mMe_p,&eccPointG,&jacPointP);
  fEcc_P2Jac(&mMe_p,&eccPointG,&jacPointR);
  
  printf("\nR_X = ");
  print_vlong(&jacPointR.vX);
  printf("\nR_Y = ");
  print_vlong(&jacPointR.vY);
  printf("\nR_Z = ");
  print_vlong(&jacPointR.vZ);
  

  if(fEcc_Jac2P(&mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nfEcc_Jac2P=G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
  }
  //
  //iRet = fEcc_DubJacP(&mMe_p,&a,&jacPointR,pvl_buff);
  iRet = fEcc_AddJacP(&mMe_p,&a,&jacPointP,&jacPointR,pvl_buff);
  
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=2G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //Q=2G
  fEcc_P2Jac(&mMe_p,&eccPointG,&jacPointQ);
  //
  iRet = fEcc_AddJacP(&mMe_p,&a,&jacPointP,&jacPointR,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=2G+G=3G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //
  iRet = fEcc_AddJacP(&mMe_p,&a,&jacPointP,&jacPointR,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=3G+G=4G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //
  iRet = fEcc_DubJacP(&mMe_p,&a,&jacPointQ,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointQ,&eccPointG) == 1){
      printf("\nQ=2*2G=4G\nQx = ");
      print_vlong(&eccPointG.vx);
      printf("\nQy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //
  /*
  d = 7A929ADE789BB9BE10ED359DD39A72C11B60961F49397EEE1D19CE9891EC3B28
  Q = [d]*G
  Qx= 7F2B49E270DB6D90D8595BEC458B50C58585BA1D4E9B788F6689DBD8E56FD80B
  Qy= 26F1B489D6701DD185C8413A977B3CBBAF64D1C593D26627DFFB101A87FF77DA
  */
  printf("\n\nGOST example:\n d = ");
  print_vlong(&d);

  iRet = fEcc_MulJacP(&mMe_p,&a,&d,&jacPointP);
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointP,&eccPointG) == 1){
      printf("\nP=[d]*G\nPx = ");
      print_vlong(&eccPointG.vx);
      printf("\nPy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //[d]*G+[l]*G
  vl_shr(&d);
  printf("\n\nGOST example:\nd/2= ");
  print_vlong(&d);

  vl_set_uchar_BE(&eccPointG.vx,str_Gx,sizeof(str_Gx));
  vl_set_uchar_BE(&eccPointG.vy,str_Gy,sizeof(str_Gy));
  fEcc_P2Jac(&mMe_p,&eccPointG,&jacPointP);
  fEcc_P2Jac(&mMe_p,&eccPointG,&jacPointQ);

  iRet = fEcc_MulJacPQ(&mMe_p,&a,&d,&jacPointQ,&d,&jacPointP);
  if(iRet == 1){
    if(fEcc_Jac2P(&mMe_p,&jacPointP,&eccPointG) == 1){
      printf("\nP=[d]*G+[l]*G\nPx = ");
      print_vlong(&eccPointG.vx);
      printf("\nPy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //
  vl_delete(&p);
  vl_delete(&a);
  vl_delete(&d);
  me_delete(&mMe_p);
  fEccPoint_delete(&eccPointG);
  fEccJacPoint_delete(&jacPointP);
  fEccJacPoint_delete(&jacPointQ);
  fEccJacPoint_delete(&jacPointR);
  
  for(i=0; i<9; i++)
    vl_delete(&pvl_buff[i]);
  
  
}

void ecc_base_init(PTBaseEccGp pEccBase){
  TVlong t;
  unsigned char str_p[] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0x97};
  unsigned char str_q[] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0x6C,0x61,0x10,0x70,0x99,0x5A,0xD1,0x00,0x45,0x84,0x1B,0x09,0xB7,0x61,0xB8,0x93};
  unsigned char str_a[] = {
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,
    0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFD,0x94};
  unsigned char str_b[] = {0xa6};
  unsigned char str_Gx[]= {0x01};
  unsigned char str_Gy[]= {
    0x8D,0x91,0xE4,0x71,0xE0,0x98,0x9C,0xDA,0x27,0xDF,0x50,0x5A,0x45,0x3F,0x2B,0x76,
    0x35,0x29,0x4F,0x2D,0xDF,0x23,0xE3,0xB1,0x22,0xAC,0xC9,0x9C,0x9E,0x9F,0x1E,0x14};

  if(pEccBase){
    vl_init(&t);
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
    //
	TPointEcc ecc_G;
	fEccPoint_init(&ecc_G);
    fEccJacPoint_init(&pEccBase->mG);
	vl_set_uchar_BE(&ecc_G.vx,str_Gx,sizeof(str_Gx));
	vl_set_uchar_BE(&ecc_G.vy,str_Gy,sizeof(str_Gy));
	fEcc_P2Jac(&pEccBase->mMe_p,&ecc_G,&pEccBase->mG);
	fEccPoint_delete(&ecc_G);
  }
}

void ecc_base_delete(PTBaseEccGp pEccBase){
  if(pEccBase){
    me_delete(&pEccBase->mMe_p);
    me_delete(&pEccBase->mMe_q);
    //
    vl_delete(&pEccBase->vb);
	vl_delete(&pEccBase->va_r);
    //
    fEccJacPoint_delete(&pEccBase->mG);
  }
}

void test_ecc_raw_2(){
  int iRet,i;
  TVlong d;//p,a
  TBaseEccGp ecc_base;
  TPointEcc eccPointG;
  TJacPointEcc jacPointP,jacPointQ,jacPointR;
  TVlong pvl_buff[9];
  
  unsigned char str_d[] = {
    0x7A,0x92,0x9A,0xDE,0x78,0x9B,0xB9,0xBE,0x10,0xED,0x35,0x9D,0xD3,0x9A,0x72,0xC1,
    0x1B,0x60,0x96,0x1F,0x49,0x39,0x7E,0xEE,0x1D,0x19,0xCE,0x98,0x91,0xEC,0x3B,0x28};
  //
  vl_init(&d);
  vl_set_uchar_BE(&d,str_d,sizeof(str_d));
  //
  ecc_base_init(&ecc_base);
  //
  fEccPoint_init(&eccPointG);

  for(i=0; i<9; i++)
    vl_init(&pvl_buff[i]);

  //
  printf("\n\nTest2, CryptoPro-A:\nBr_p.M  = ");
  print_vlong(&ecc_base.mMe_p.M);
  //printf("\nBr_p.mu = ");
  //print_vlong(&ecc_base.mMe_p.mu);
  printf("\nBr_q.M  = ");
  print_vlong(&ecc_base.mMe_q.M);
  //printf("\nBr_q.mu = ");
  //print_vlong(&ecc_base.mBr_q.mu);
  //printf("\na  = ");
  //print_vlong(&ecc_base.va);
  printf("\na_r= ");
  print_vlong(&ecc_base.va_r);
  printf("\nb  = ");
  print_vlong(&ecc_base.vb);
  printf("\nGx = ");
  print_vlong(&ecc_base.mG.vX);
  printf("\nGy = ");
  print_vlong(&ecc_base.mG.vY);
  printf("\nGz = ");
  print_vlong(&ecc_base.mG.vZ);
  //
  fEccJacPoint_init(&jacPointP);
  fEccJacPoint_init(&jacPointQ);
  fEccJacPoint_init(&jacPointR);
  //
  fEcc_Jac2Jac(&jacPointP,&ecc_base.mG);
  fEcc_Jac2Jac(&jacPointR,&ecc_base.mG);
  //
  if(fEcc_Jac2P(&ecc_base.mMe_p,&ecc_base.mG,&eccPointG) == 1){
    printf("\nR = G\nRx = ");
    print_vlong(&eccPointG.vx);
    printf("\nRy = ");
    print_vlong(&eccPointG.vy);
  }
  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point G = %d",iRet);
  
  iRet = fEcc_DubJacP(&ecc_base.mMe_p,&ecc_base.va_r,&jacPointR,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=2G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }

  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point R = %d",iRet);

  //Q=2G
  //fEcc_P2Jac(&ecc_base.mMe_p,&eccPointG,&jacPointQ);
  fEcc_Jac2Jac(&jacPointQ,&jacPointR);
  
  //
  iRet = fEcc_AddJacP(&ecc_base.mMe_p,&ecc_base.va_r,&jacPointP,&jacPointR,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=2G+G=3G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }

  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point R = %d",iRet);
  //
  iRet = fEcc_AddJacP(&ecc_base.mMe_p,&ecc_base.va_r,&jacPointP,&jacPointR,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointR,&eccPointG) == 1){
      printf("\nR=3G+G=4G\nRx = ");
      print_vlong(&eccPointG.vx);
      printf("\nRy = ");
      print_vlong(&eccPointG.vy);
    }
  }

  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point R = %d",iRet);

  //
  iRet = fEcc_DubJacP(&ecc_base.mMe_p,&ecc_base.va_r,&jacPointQ,pvl_buff);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointQ,&eccPointG) == 1){
      printf("\nQ=2*2G=4G\nQx = ");
      print_vlong(&eccPointG.vx);
      printf("\nQy = ");
      print_vlong(&eccPointG.vy);
    }
  }

  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point Q = %d",iRet);
  //
  printf("\n\nCryptoPro-A:");
  //q*P = 0
  iRet = fEcc_MulJacP(&ecc_base.mMe_p,&ecc_base.va_r,&ecc_base.mMe_q.M,&jacPointQ);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointQ,&eccPointG) == 1){
      printf("\nQ=[q]*Q(=0)\nQx = ");
      print_vlong(&eccPointG.vx);
      printf("\nQy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  //
  printf("\n\nd  = ");
  print_vlong(&d);

  iRet = fEcc_MulJacP(&ecc_base.mMe_p,&ecc_base.va_r,&d,&jacPointP);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointP,&eccPointG) == 1){
      printf("\nP=[d]*G\nPx = ");
      print_vlong(&eccPointG.vx);
      printf("\nPy = ");
      print_vlong(&eccPointG.vy);
    }
  }

  iRet = fEcc_CheckPoint(&ecc_base,&eccPointG);
  printf("\nCheck point P = %d",iRet);
  //
  //[d]*G+[l]*G
  vl_shr(&d);
  printf("\n\nCryptoPro-A:\nd/2= ");
  print_vlong(&d);

  fEcc_Jac2Jac(&jacPointP,&ecc_base.mG);
  fEcc_Jac2Jac(&jacPointQ,&ecc_base.mG);

  iRet = fEcc_MulJacPQ(&ecc_base.mMe_p,&ecc_base.va_r,&d,&jacPointQ,&d,&jacPointP);
  if(iRet == 1){
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointP,&eccPointG) == 1){
      printf("\nP=[d]*G+[l]*G\nPx = ");
      print_vlong(&eccPointG.vx);
      printf("\nPy = ");
      print_vlong(&eccPointG.vy);
    }
  }
  
  int iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = fEcc_MulJacPQ(&ecc_base.mMe_p,&ecc_base.va_r,&d,&jacPointQ,&d,&jacPointP);
    if(iRet != 1)
      printf("\n Error fEcc_MulJacPQ!");
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointP,&eccPointG) != 1){
		printf("\n Error fEcc_Jac2P!");
	}
  }
  printf("\n\n %d*(fEcc_MulJacPQ+fEcc_Jac2P) = %d ms", i, GetTickCountMy() - iTstart);
  
  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = fEcc_MulJacP(&ecc_base.mMe_p,&ecc_base.va_r,&d,&jacPointP);
    if(iRet != 1)
      printf("\n Error fEcc_MulJacP!");
    if(fEcc_Jac2P(&ecc_base.mMe_p,&jacPointP,&eccPointG) != 1){
		printf("\n Error fEcc_Jac2P!");
	}
  }
  printf("\n\n %d*(fEcc_MulJacP+fEcc_Jac2P) = %d ms", i, GetTickCountMy() - iTstart);


  printf("\n");
  //
  for(i=0; i<9; i++)
    vl_delete(&pvl_buff[i]);

  vl_delete(&d);
  fEccPoint_delete(&eccPointG);
  fEccJacPoint_delete(&jacPointP);
  fEccJacPoint_delete(&jacPointQ);
  fEccJacPoint_delete(&jacPointR);
  ecc_base_delete(&ecc_base);
}
//------------------------------------------------------------------------------
void test_flow_gost3410_256(){
  int iRet,i,iTstart;
  uint8_t pubKey[2*32],privKey[32],hash[32],r[32],s[32],ukm[32];

  printf("\n\n--- Flow test id-GostR3410-2001-CryptoPro-A-ParamSet ---\n");
  
  PTBaseEccGp pEccBase = cryptoEccBase256Init();
  if(!pEccBase){
    printf("\n Error cryptoEccBase256Init!");
    return;
  }
  //*
  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoGenerateKeyPair(pEccBase,&pubKey[0],&pubKey[32],privKey);
    if(iRet != 1)
      printf("\n Error cryptoGenerateKeyPair!");
  }
  printf("\n\n 1000*cryptoGenerateKeyPair = %d ms", GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<100000; i++){
  iRet = cryptoCheckPubKeyValue(pEccBase,&pubKey[0],&pubKey[32]);
    if(iRet != 1)
      printf("\n Error cryptoCheckPubKeyValue!");
  }
  printf("\n\n %d*cryptoCheckPubKeyValue = %d ms", i, GetTickCountMy() - iTstart);
  //*/
  getRNG(hash,sizeof(hash));
  
  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoMakeSignature(pEccBase,hash,privKey,r,s);
    if(iRet != 1)
      printf("\n Error cryptoMakeSignature!");
  }
  printf("\n\n %d*cryptoMakeSignature = %d ms", i, GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoVerifySignature(pEccBase,r,s,hash,&pubKey[0],&pubKey[32]);
    if(iRet != 1)
      printf("\n Error cryptoVerifySignature!");
  }
  printf("\n\n %d*cryptoVerifySignature = %d ms", i, GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoMakeSignature(pEccBase,hash,privKey,r,s);
    if(iRet != 1)
      printf("\n Error cryptoMakeSignature!");
    iRet = cryptoVerifySignature(pEccBase,r,s,hash,&pubKey[0],&pubKey[32]);
    if(iRet != 1)
      printf("\n Error cryptoVerifySignature!");
  }
  printf("\n\n %d*(cryptoMakeSignature+cryptoVerifySignature) = %d ms", i, GetTickCountMy() - iTstart);
  
  getRNG(ukm,sizeof(ukm));
  ukm[0]&=0x7F;

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoVkoGenerate(pEccBase,privKey,ukm,&pubKey[0],&pubKey[32]);
    if(iRet != 1)
      printf("\n Error cryptoVkoGenerate!");
  }
  printf("\n\n %d*cryptoVkoGenerate = %d ms", i, GetTickCountMy() - iTstart);

  // cryptoVkoGenerate
  test_cryptoVkoGenerate(pEccBase);
  
  cryptoEccBaseDelete(&pEccBase);
}

void test_cryptoVkoGenerate(PTBaseEccGp pEccBase)
{
	int iRet;
	uint8_t pubKey1[2*64],privKey1[64],ukm[64];
	uint8_t pubKey2[2*64],privKey2[64];
	printf("\n\n--- Flow test cryptoVkoGenerate ---\n");
	if(!pEccBase){
		printf("\n Error pEccBase");
		return;
	}
	getRNG(ukm,pEccBase->iModSize);
	ukm[0]&=0x7F;

	iRet = cryptoGenerateKeyPair(pEccBase,&pubKey1[0],&pubKey1[0+pEccBase->iModSize],privKey1);
	if(1 == iRet)
	{
		printf("\n\nprivKey1 = ");
		print_hex(privKey1,pEccBase->iModSize);
		printf("\nPublKey1X = ");
		print_hex(&pubKey1[0],pEccBase->iModSize);
		printf("\nPublKey1Y = ");
		print_hex(&pubKey1[0+pEccBase->iModSize],pEccBase->iModSize);
	}
	else
		printf("\n Error cryptoGenerateKeyPair!");
	
	iRet = cryptoGenerateKeyPair(pEccBase,&pubKey2[0],&pubKey2[0+pEccBase->iModSize],privKey2);
	if(1 == iRet)
	{
		printf("\n\nprivKey2 = ");
		print_hex(privKey2,pEccBase->iModSize);
		printf("\nPublKey2X = ");
		print_hex(&pubKey2[0],pEccBase->iModSize);
		printf("\nPublKey2Y = ");
		print_hex(&pubKey2[0+pEccBase->iModSize],pEccBase->iModSize);
	}
	else
		printf("\n Error cryptoGenerateKeyPair!");

	printf("\n\nukm = ");
	print_hex(ukm,pEccBase->iModSize);

	printf("\n\ncryptoVkoGenerate privKey1 + pubKey2");
	iRet = cryptoVkoGenerate(pEccBase,privKey1,ukm,&pubKey2[0],&pubKey2[0+pEccBase->iModSize]);
	if(iRet != 1)
		printf("\n Error cryptoVkoGenerate!");
	//
	printf("\n\nVkoKeyX = ");
	print_hex(&pubKey2[0],pEccBase->iModSize);
	printf("\nVkoKeyY = ");
	print_hex(&pubKey2[0+pEccBase->iModSize],pEccBase->iModSize);
	//
	printf("\n\ncryptoVkoGenerate privKey2 + pubKey1");
	iRet = cryptoVkoGenerate(pEccBase,privKey2,ukm,&pubKey1[0],&pubKey1[0+pEccBase->iModSize]);
	if(iRet != 1)
		printf("\n Error cryptoVkoGenerate!");
	//
	printf("\n\nVkoKeyX = ");
	print_hex(&pubKey1[0],pEccBase->iModSize);
	printf("\nVkoKeyY = ");
	print_hex(&pubKey1[0+pEccBase->iModSize],pEccBase->iModSize);
}

void test_flow_gost3410_512(){
  int iRet,i,iTstart;
  uint8_t pubKey[2*64],privKey[64],hash[64],r[64],s[64],ukm[64];

  printf("\n\n--- Flow test id-tc26-gost-3410-12-512-paramSetA ---\n");
  
  PTBaseEccGp pEccBase = cryptoEccBase512Init();
  if(!pEccBase){
    printf("\n Error cryptoEccBase512Init!");
    return;
  }
  //*
  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoGenerateKeyPair(pEccBase,&pubKey[0],&pubKey[64],privKey);
    if(iRet != 1)
      printf("\n Error cryptoGenerateKeyPair!");
  }
  printf("\n\n 1000*cryptoGenerateKeyPair = %d ms", GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<100000; i++){
  iRet = cryptoCheckPubKeyValue(pEccBase,&pubKey[0],&pubKey[64]);
    if(iRet != 1)
      printf("\n Error cryptoCheckPubKeyValue!");
  }
  printf("\n\n %d*cryptoCheckPubKeyValue = %d ms", i, GetTickCountMy() - iTstart);
  //*/
  getRNG(hash,sizeof(hash));
  
  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoMakeSignature(pEccBase,hash,privKey,r,s);
    if(iRet != 1)
      printf("\n Error cryptoMakeSignature!");
  }
  printf("\n\n %d*cryptoMakeSignature = %d ms", i, GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoVerifySignature(pEccBase,r,s,hash,&pubKey[0],&pubKey[64]);
    if(iRet != 1)
      printf("\n Error cryptoVerifySignature!");
  }
  printf("\n\n %d*cryptoVerifySignature = %d ms", i, GetTickCountMy() - iTstart);

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoMakeSignature(pEccBase,hash,privKey,r,s);
    if(iRet != 1)
      printf("\n Error cryptoMakeSignature!");
    iRet = cryptoVerifySignature(pEccBase,r,s,hash,&pubKey[0],&pubKey[64]);
    if(iRet != 1)
      printf("\n Error cryptoVerifySignature!");
  }
  printf("\n\n %d*(cryptoMakeSignature+cryptoVerifySignature) = %d ms", i, GetTickCountMy() - iTstart);
  
  getRNG(ukm,sizeof(ukm));
  ukm[0]&=0x7F;

  iTstart = GetTickCountMy();
  for(i = 0; i<1000; i++){
    iRet = cryptoVkoGenerate(pEccBase,privKey,ukm,&pubKey[0],&pubKey[64]);
    if(iRet != 1)
      printf("\n Error cryptoVkoGenerate!");
  }
  printf("\n\n %d*cryptoVkoGenerate = %d ms", i, GetTickCountMy() - iTstart);

  // cryptoVkoGenerate
  test_cryptoVkoGenerate(pEccBase);
  
  cryptoEccBaseDelete(&pEccBase);
}

/*
     163 30  159:  SEQUENCE {
     166 06    7:   OBJECT IDENTIFIER
                :    id-GostR3410-2001-CryptoPro-A-ParamSet
     175 30  147:   SEQUENCE {
a     178 02   33:    INTEGER
                :     00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
                :     FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FD
                :     94
b     213 02    2:    INTEGER 166
p     217 02   33:    INTEGER
                :     00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
                :     FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF FD
                :     97
q     252 02   33:    INTEGER
                :     00 FF FF FF FF FF FF FF FF FF FF FF FF FF FF FF
                :     FF 6C 61 10 70 99 5A D1 00 45 84 1B 09 B7 61 B8
                :     93
Gx     287 02    1:    INTEGER 1
Gy     290 02   33:    INTEGER
                :     00 8D 91 E4 71 E0 98 9C DA 27 DF 50 5A 45 3F 2B
                :     76 35 29 4F 2D DF 23 E3 B1 22 AC C9 9C 9E 9F 1E
                :     14
                :    }
                :   }


//------

p = 8000000000000000000000000000000000000000000000000000000000000431
a = 7
b = 5FBFF498AA938CE739B8E022FBAFEF40563F6E6A3472FC2A514C0CE9DAE23B7E
q = 8000000000000000000000000000000150FE8A1892976154C59CFC193ACCF5B3


Gx= 2
Gy= 08E2A8A0E65147D4BD6316030E16D19C85C97F0A9CA267122B96ABBCEA7E8FC8
    08E2A8A0E65147D4BD6316030E16D19C85C97F0A9CA267122B96ABBCEA7E8FC8

d = 7A929ADE789BB9BE10ED359DD39A72C11B60961F49397EEE1D19CE9891EC3B28

Q = [d]*G

Qx= 7F2B49E270DB6D90D8595BEC458B50C58585BA1D4E9B788F6689DBD8E56FD80B
Qy= 26F1B489D6701DD185C8413A977B3CBBAF64D1C593D26627DFFB101A87FF77DA
*/

#pragma warning(disable: 4996)

//------ rnd emulate -------

static
int my_rand(int *pnext){
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

static int rnd_seed = 0;

void rngInit() 
{
  rnd_seed ^= (int)(GetTickCountMy()*GetTickCountMy());
  //rnd_seed = 0x39232B61;
  //rnd_seed = 0x09D4B790;
  //rnd_seed = 0x9624AE24;
  //rnd_seed = 0xB164BE40;
  //rnd_seed = 0x5BE888A4;
  printf("\nrngInit:seed = 0x%08X\n",(unsigned int)rnd_seed);
}

int getRNG(uint8_t *rnd, int iSize)
{
  int i = 0;
  if(rnd){
    for(i = 0; i<iSize; i++){
      rnd[i] = (uint8_t)my_rand(&rnd_seed);
    }
  }
  return i;
}




void print_hex_str(const unsigned char *strHex, unsigned int uLenHex, char *strRes)
{
  unsigned int i;
#ifdef __BORLANDC__

  typedef WINUSERAPI int (WINAPIV *TwsprintfA)(LPSTR, LPCSTR, ...);
  TwsprintfA pwsprintfA;
  pwsprintfA = (TwsprintfA)::GetProcAddress(::GetModuleHandle("user32"),"wsprintfA");

  if(pwsprintfA)
#endif
  {
    for(i=0;i<uLenHex;i++)
    {
      #ifdef __BORLANDC__
        pwsprintfA(&strRes[i*2],"%02X",(unsigned char)strHex[i]);
      #else
        //wsprintfA(&strRes[i*2],"%02X",(unsigned char)strHex[i]);
        sprintf(&strRes[i*2],"%02X",(unsigned char)strHex[i]);
      #endif
      strRes[i*2+2] = 0;
    }
  };
};

void print_hex(const unsigned char *strHex, unsigned int uLenHex)
{
  char str[512];
  str[0] = 0;
  print_hex_str(strHex,uLenHex,str);
  printf("%s",str);
};

void print_vlong(PTVlong pvl){
  char str[256],strx[128];
  unsigned int i;
  str[0] = 0x30;
  str[1] = 0;
  i = sizeof(strx);
  if(vl_get_uchar_BE(strx, &i, pvl) == 1){
    print_hex_str(strx,i,str);
    printf("%s",str);
  }

}

#ifdef __GNUC__
#include <sys/time.h>
int GetTickCountMy()//time in ms
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}
#endif

#ifdef WIN32
static LARGE_INTEGER frequency = { 0 };
int GetTickCountMy() {//time in ms
  LARGE_INTEGER start;
  //
  QueryPerformanceCounter(&start);
  //
  if (frequency.QuadPart == 0)
    QueryPerformanceFrequency(&frequency);
  //
  float fstart_ms = start.QuadPart * 1000.0f / frequency.QuadPart;
  //
  return (int)fstart_ms;
}
#endif//WIN32
