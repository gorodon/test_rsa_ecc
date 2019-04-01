// ECC_GP
#ifndef _ECC_GP
#define _ECC_GP

#include "Vlong_Gp.h"

//������������� ������ ���� y^2=x^3+a*x+b mod p   G(Fp) Ep(a,b)
//
typedef 
struct TPointEcc{//(x,y) ������� ����������
  TVlong vx,vy;
} TPointEcc, *PTPointEcc;
//
typedef 
struct TJacPointEcc{//(X*R,Y*R,Z*R)->(x,y) (X/Z^2,Y/Z^3) ����������� �������������
  TVlong vX,vY,vZ;
} TJacPointEcc, *PTJacPointEcc;
//
typedef 
struct TBaseEccGp{
  //����� ������ � ������ (32(256 ���) ��� 64(512���))
  int iModSize;
  //������ ������������� ������ p (� �������)
  TMontExp mMe_p;
  //�������(�������) ������������� ������ q = #Ep(a,b)/h (h=1 - cofactor)
  TMontExp mMe_q;
  //
  //������������ ������������� ������ a,b
  TVlong va_r;// = a*R mod p
  TVlong vb;
  //
  //������� ����� �� ������������� ������ G(Gx*R,Gy*R,1*R) (��������� ������������ �������� q)
  TJacPointEcc mG;
  //��������������� ����� ������������� ������ 3G,7G
  //TJacPointEcc m3G, m7G;
  //
}TBaseEccGp,*PTBaseEccGp;

//init-delete
void fEccPoint_init(PTPointEcc pEcc);
void fEccPoint_delete(PTPointEcc pEcc);

void fEccJacPoint_init(PTJacPointEcc pEccJac);
void fEccJacPoint_delete(PTJacPointEcc pEccJac);

//�������� �������������� ����� pEcc ������������� ������ pEccBase
int fEcc_CheckPoint(const PTBaseEccGp pEccBase,const PTPointEcc pEcc);
// �������������� (x,y)->(X,Y,Z) (x,y,1) ����������� �������������
int fEcc_P2Jac(
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTPointEcc pAPoint,       //(IN)  Affine point
          PTJacPointEcc pJacPoint);       //(OUT) Jacobian point
// �������������� (X,Y,Z)->(x,y) (X/Z^2,Y/Z^3)
int fEcc_Jac2P(
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTJacPointEcc pJacPoint,  //(IN)  Jacobian point
          PTPointEcc pAPoint);            //(OUT) Affine point
// copy 
int fEcc_Jac2Jac(
          PTJacPointEcc pJacPointDst,         //(OUT) Jacobian point
          const PTJacPointEcc pJacPointSrc);  //(In) Jacobian point
//
// �������� ���� ����� �� ���-� ������ 
int fEcc_AddJacP(          
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  ����������� ������������� ������
          const PTJacPointEcc pJacPointQ, //(IN)  Jacobian point Q
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=P+Q
          TVlong pvl[9]);                 //(IN)  ����� �������� TVlong
// �������� ���� ����� �� ���-� ������ 
int fEcc_DubJacP(          
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  ����������� ������������� ������
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=2*P
          TVlong pvl[9]);                 //(IN)  ����� �������� TVlong

// d-������� ���������� ����� �� ���-� ������ [d]*P
int fEcc_MulJacP(
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  ����������� ������������� ������
          const PTVlong pvd,              //(IN)  [d]-������� �����
          PTJacPointEcc pJacPointP);      //(IN/OUT)  Jacobian point P P=[d]*P
//
// d,l-������� ���������� ����� �� ���-� ������ [d]*P+[l]*Q
int fEcc_MulJacPQ(
          const PTMontExp pMe_p, //(IN)  ������ ��� ���������� ������� �������� a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  ����������� ������������� ������
          const PTVlong pvl,              //(IN)  [l]-������� �����
          PTJacPointEcc pJacPointQ,       //(IN)  Jacobian point Q(will be normalise)
          const PTVlong pvd,              //(IN)  [d]-������� �����
          PTJacPointEcc pJacPointP);      //(IN/OUT)  Jacobian point P P=[d]*P+[l]*Q
//          
#endif//#ifndef _ECC_GP