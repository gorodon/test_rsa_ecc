// ECC_GP
#ifndef _ECC_GP
#define _ECC_GP

#include "Vlong_Gp.h"

//Эллиптические кривые вида y^2=x^3+a*x+b mod p   G(Fp) Ep(a,b)
//
typedef 
struct TPointEcc{//(x,y) афинные координаты
  TVlong vx,vy;
} TPointEcc, *PTPointEcc;
//
typedef 
struct TJacPointEcc{//(X*R,Y*R,Z*R)->(x,y) (X/Z^2,Y/Z^3) проективное представление
  TVlong vX,vY,vZ;
} TJacPointEcc, *PTJacPointEcc;
//
typedef 
struct TBaseEccGp{
  //длина модуля в байтах (32(256 бит) или 64(512бит))
  int iModSize;
  //Модуль эллиптической кривой p (в объекте)
  TMontExp mMe_p;
  //Порядок(базовый) эллиптической кривой q = #Ep(a,b)/h (h=1 - cofactor)
  TMontExp mMe_q;
  //
  //Коэффициенты эллиптической кривой a,b
  TVlong va_r;// = a*R mod p
  TVlong vb;
  //
  //Базовая точка на эллиптической кривой G(Gx*R,Gy*R,1*R) (генератор подмножества мощности q)
  TJacPointEcc mG;
  //Вспомогательные точки эллиптической кривой 3G,7G
  //TJacPointEcc m3G, m7G;
  //
}TBaseEccGp,*PTBaseEccGp;

//init-delete
void fEccPoint_init(PTPointEcc pEcc);
void fEccPoint_delete(PTPointEcc pEcc);

void fEccJacPoint_init(PTJacPointEcc pEccJac);
void fEccJacPoint_delete(PTJacPointEcc pEccJac);

//Проверка принадлежности точки pEcc эллиптической кривой pEccBase
int fEcc_CheckPoint(const PTBaseEccGp pEccBase,const PTPointEcc pEcc);
// Преобразование (x,y)->(X,Y,Z) (x,y,1) проективное представление
int fEcc_P2Jac(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTPointEcc pAPoint,       //(IN)  Affine point
          PTJacPointEcc pJacPoint);       //(OUT) Jacobian point
// Преобразование (X,Y,Z)->(x,y) (X/Z^2,Y/Z^3)
int fEcc_Jac2P(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTJacPointEcc pJacPoint,  //(IN)  Jacobian point
          PTPointEcc pAPoint);            //(OUT) Affine point
// copy 
int fEcc_Jac2Jac(
          PTJacPointEcc pJacPointDst,         //(OUT) Jacobian point
          const PTJacPointEcc pJacPointSrc);  //(In) Jacobian point
//
// Сложение двух точек на элл-й кривой 
int fEcc_AddJacP(          
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          const PTJacPointEcc pJacPointQ, //(IN)  Jacobian point Q
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=P+Q
          TVlong pvl[9]);                 //(IN)  буфер объектов TVlong
// Удвоение двух точек на элл-й кривой 
int fEcc_DubJacP(          
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=2*P
          TVlong pvl[9]);                 //(IN)  буфер объектов TVlong

// d-кратная композиция точки на элл-й кривой [d]*P
int fEcc_MulJacP(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          const PTVlong pvd,              //(IN)  [d]-длинное целое
          PTJacPointEcc pJacPointP);      //(IN/OUT)  Jacobian point P P=[d]*P
//
// d,l-кратная композиция точек на элл-й кривой [d]*P+[l]*Q
int fEcc_MulJacPQ(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          const PTVlong pvl,              //(IN)  [l]-длинное целое
          PTJacPointEcc pJacPointQ,       //(IN)  Jacobian point Q(will be normalise)
          const PTVlong pvd,              //(IN)  [d]-длинное целое
          PTJacPointEcc pJacPointP);      //(IN/OUT)  Jacobian point P P=[d]*P+[l]*Q
//          
#endif//#ifndef _ECC_GP