#include <string.h>
#include <stdlib.h>
#include "ecc_Gp.h"

//Warning: single-thread version !!!

//Проверка принадлежности точки pEcc эллиптической кривой pEccBase
//y^2=x^3+a*x+b mod p
int fEcc_CheckPoint(const PTBaseEccGp pEccBase,const PTPointEcc pEcc){
	int iRet = 0;
	if(pEccBase && pEcc){
		TVlong r1,r2;
		const PTVlong pvp = &pEccBase->mMe_p.M;
		//
		vl_init(&r1);
		vl_init(&r2);
		//r1 = Ex^2 mod p
		me_modmul(&r1,&pEcc->vx,&pEcc->vx,&pEccBase->mMe_p);
		//r1 = Ex^2 + a mod p
		vl_copy(&r2,&pEccBase->va_r);
		me_mont(&r2,&pEccBase->mMe_p);
		vl_add(&r1,&r2);
		//vl_add(&r1,&pEccBase->va);
		if(vl_cf(&r1,pvp)>=0)
			vl_sub(&r1,pvp);
		//r2 = Ex*(Ex^2 + a) mod p
		me_modmul(&r2,&r1,&pEcc->vx,&pEccBase->mMe_p);
		//r2 = Ex^3 + a*Ex + b mod p
		vl_add(&r2,&pEccBase->vb);
		if(vl_cf(&r2,pvp)>=0)
			vl_sub(&r2,pvp);
		//
		//r1 = Ey^2 mod p
		me_modmul(&r1,&pEcc->vy,&pEcc->vy,&pEccBase->mMe_p);
		//
		if(vl_cf(&r1,&r2)==0)
			iRet = 1;
		//
		vl_delete(&r1);
		vl_delete(&r2);
	}
	return iRet;
}

//init-delete
void fEccPoint_init(PTPointEcc pEcc){
  if(pEcc){
    vl_init(&pEcc->vx);
    vl_init(&pEcc->vy);
  }
}
void fEccPoint_delete(PTPointEcc pEcc){
  if(pEcc){
    vl_delete(&pEcc->vx);
    vl_delete(&pEcc->vy);
  }
}

void fEccJacPoint_init(PTJacPointEcc pEccJac){
  if(pEccJac){
    vl_init(&pEccJac->vX);
    vl_init(&pEccJac->vY);
    vl_init(&pEccJac->vZ);
  }
}

void fEccJacPoint_delete(PTJacPointEcc pEccJac){
  if(pEccJac){
    vl_delete(&pEccJac->vX);
    vl_delete(&pEccJac->vY);
    vl_delete(&pEccJac->vZ);
  }
}

// copy 
int fEcc_Jac2Jac(
          PTJacPointEcc pJacPointDst,         //(OUT) Jacobian point
          const PTJacPointEcc pJacPointSrc){  //(In) Jacobian point
  int iRet = 0;
  //
  if(pJacPointDst && pJacPointSrc){
    vl_copy(&pJacPointDst->vX,&pJacPointSrc->vX);
    vl_copy(&pJacPointDst->vY,&pJacPointSrc->vY);
    vl_copy(&pJacPointDst->vZ,&pJacPointSrc->vZ);
    iRet = 1;
  }
  //
  return iRet;
}
// Преобразование (x,y)->(X,Y,Z) (x*R,y*R,1*R) проективное представление
int fEcc_P2Jac(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b mod p
          const PTPointEcc pAPoint,       //(IN)  Affine point
          PTJacPointEcc pJacPoint){       //(OUT) Jacobian point
	int iRet = 0;
	if(pMe_p && pAPoint && pJacPoint){
		vl_copy(&pJacPoint->vX,&pAPoint->vx);
		me_mul(&pJacPoint->vX,&pMe_p->r2,pMe_p);
		//
		vl_copy(&pJacPoint->vY,&pAPoint->vy);
		me_mul(&pJacPoint->vY,&pMe_p->r2,pMe_p);
		//
		//
		vl_copy(&pJacPoint->vZ,&pMe_p->r2);
		me_mont(&pJacPoint->vZ,pMe_p);
		//
		iRet = 1;
	}
	return iRet;
}
// Преобразование (X,Y,Z)->(x,y) (X/Z^2,Y/Z^3)
int fEcc_Jac2P(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b mod p
          const PTJacPointEcc pJacPoint,  //(IN)  Jacobian point
          PTPointEcc pAPoint){            //(OUT) Affine point 
	int iRet = 0;
	TVlong t;
	if(pMe_p && pAPoint && pJacPoint){
	//
	if(pJacPoint->vZ.value.n == 0){//Z==0 - точка в бесконечночти (0,0)
		vl_clear(&pAPoint->vx);//x=0
		vl_clear(&pAPoint->vy);//y=0
	}else{
		vl_init(&t);
		//
		vl_modinv(&t,&pJacPoint->vZ,&pMe_p->M); //t = (Z*R)^-1 mod M
		me_mul(&t,&pMe_p->r2,pMe_p);// t = Z^-1*R^-1*R^2*R^-1 = Z^-1 mod M
		//
		me_modmul(&pAPoint->vy,&t,&t,pMe_p);    //pAPoint->vy = Z^-2 mod M
		vl_copy(&pAPoint->vx,&pJacPoint->vX);
		me_mul(&pAPoint->vx,&pAPoint->vy,pMe_p);// (X*R)*(Z^-2)*R^-1 mod M = X*Z^-2 mod M
		//
		me_modmul(&t,&pAPoint->vy,&t,pMe_p);//t = Z^-3 mod M
		vl_copy(&pAPoint->vy,&pJacPoint->vY);
		me_mul(&pAPoint->vy,&t,pMe_p);// (Y*R)*(Z^-3)*R^-1 mod M = Y*Z^-3 mod M
		//
		vl_delete(&t);
	}
	//
	iRet = 1;
	}
	return iRet;
}

// Отрицание (P =-P) точки в проективных координатах(X,Y,Z)->(X,-Y,Z)
int fEcc_JacInvert(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b mod p
          PTJacPointEcc pJacPoint,        //(IN/OUT)  Jacobian point
          TVlong pvl[9]){                 //(IN)  буфер объектов TVlong
  int iRet = 0;
  if(pMe_p && pJacPoint){
    iRet = 1;
    vl_copy(&pvl[0],&pMe_p->M);     //buf = p
    vl_usub(&pvl[0],&pJacPoint->vY);//buf = p - Vy
    vl_copy(&pJacPoint->vY,&pvl[0]);//Vy = buf
  }
  return iRet;
}

// Сложение двух точек на элл-й кривой 
int fEcc_AddJacP(          
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой (a*R)
          const PTJacPointEcc pJacPointQ, //(IN)  Jacobian point Q
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=P+Q
          TVlong pvl[9]){                 //(IN)  буфер объектов TVlong
	int iRet = 0;
	PTVlong pvp;
	//
	if(pMe_p && pva && pJacPointQ && pJacPointP){
		//
		iRet = 1;
		pvp = &(pMe_p->M);//mod p
		//
		if(pJacPointQ->vZ.value.n == 0){//Qz=0 - точка в бесконечночти (0,0) P=P+0=P
			return iRet;
		}
		//
		if(pJacPointP->vZ.value.n == 0){//Pz=0 - точка в бесконечночти (0,0) P=0+Q=Q
			//
			vl_copy(&pJacPointP->vX,&pJacPointQ->vX);
			vl_copy(&pJacPointP->vY,&pJacPointQ->vY);
			vl_copy(&pJacPointP->vZ,&pJacPointQ->vZ);
			//
			return iRet;
		}
		//
		//сложение в проективном представлении
		//l0 = Pz^2
		vl_copy(&pvl[0],&pJacPointP->vZ);
		me_mul(&pvl[0],&pvl[0],pMe_p);
		//l3 = Qy*Pz^3 = Qy*Pz*l0
		vl_copy(&pvl[3],&pJacPointP->vZ);
		me_mul(&pvl[3],&pJacPointQ->vY,pMe_p);
		me_mul(&pvl[3],&pvl[0],pMe_p);
		//l0 = Qx*Pz^2
		me_mul(&pvl[0],&pJacPointQ->vX,pMe_p);
		//l1 = Qz^2
		vl_copy(&pvl[1],&pJacPointQ->vZ);
		me_mul(&pvl[1],&pvl[1],pMe_p);
		//l4 = Py*Qz^3 = Py*Qz*l1
		vl_copy(&pvl[4],&pJacPointQ->vZ);
		me_mul(&pvl[4],&pJacPointP->vY,pMe_p);
		me_mul(&pvl[4],&pvl[1],pMe_p);
		//l1 = Px*Qz^2
		me_mul(&pvl[1],&pJacPointP->vX,pMe_p);
		// корректное сравнение точек в проективном представлении (X,Y,Z) = (X*t^2,Y*t^3,Z*t)
		// Qx*Pz^2 == Px*Qz^2
		if (vl_cf(&pvl[1], &pvl[0]) == 0) //Px == Qx (affine)
		{
			//Qy*Pz^3 == Py*Qz^3
			if (vl_cf(&pvl[3], &pvl[4]) == 0)//Py == Qy (affine)
			{
				//Q==P
				return fEcc_DubJacP(pMe_p, pva, pJacPointP, pvl);//P=2*P
			}
			else 
			{	//P=P-P=0
				vl_clear(&pJacPointP->vZ);//Pz_out = 0
				return iRet;
			}
		}
		//
		//l2 = l1 - l0
		vl_copy(&pvl[2],&pvl[1]);
		if(vl_cf(&pvl[2],&pvl[0]) < 0)//l2 < l0
			vl_add(&pvl[2],pvp);
		vl_usub(&pvl[2],&pvl[0]);
		//l6 = l1 + l0
		vl_copy(&pvl[6],&pvl[1]);
		vl_add(&pvl[6],&pvl[0]);
		if(vl_cf(&pvl[6],pvp)>=0)
			vl_usub(&pvl[6],pvp);
		//l5 = l4 - l3
		vl_copy(&pvl[5],&pvl[4]);
		if(vl_cf(&pvl[5],&pvl[3]) < 0)//l5 < l3
			vl_add(&pvl[5],pvp);
		vl_usub(&pvl[5],&pvl[3]);
		//l7 = l4 + l3
		vl_copy(&pvl[7],&pvl[4]);
		vl_add(&pvl[7],&pvl[3]);
		if(vl_cf(&pvl[7],pvp)>=0)
			vl_usub(&pvl[7],pvp);
		//Pz_out=Pz*Qz*l2
		me_mul(&pJacPointP->vZ,&pJacPointQ->vZ,pMe_p);
		me_mul(&pJacPointP->vZ,&pvl[2],pMe_p);
		//Px_out=l5*l5-l6*l2*l2
		vl_copy(&pJacPointP->vX,&pvl[5]);
		me_mul(&pJacPointP->vX,&pvl[5],pMe_p);//Px = l5*l5
		vl_copy(&pvl[0],&pvl[2]);
		me_mul(&pvl[0],&pvl[2],pMe_p);//l0=l2*l2
		vl_copy(&pvl[8],&pvl[0]);
		me_mul(&pvl[8],&pvl[6],pMe_p);//l8=l6*l2*l2
		
		if(vl_cf(&pJacPointP->vX,&pvl[8]) < 0)//x < l8
			vl_add(&pJacPointP->vX,pvp);
		vl_usub(&pJacPointP->vX,&pvl[8]);
		//l8=l6*l2*l2-2*Xout
		if(vl_cf(&pvl[8],&pJacPointP->vX) < 0)//l8 < x
			vl_add(&pvl[8],pvp);
		vl_usub(&pvl[8],&pJacPointP->vX);//l8 -= x
		if(vl_cf(&pvl[8],&pJacPointP->vX) < 0)//l8 < x
			vl_add(&pvl[8],pvp);
		vl_usub(&pvl[8],&pJacPointP->vX);//l8 -= x
		//Py_out=(l8*l5-l7*l2^3)/2
		vl_copy(&pJacPointP->vY,&pvl[5]);
		me_mul(&pJacPointP->vY,&pvl[8],pMe_p);//pY = l5*l8
		me_mul(&pvl[7],&pvl[0],pMe_p);//l7=l7*l2^2
		me_mul(&pvl[7],&pvl[2],pMe_p);//l7=l7*l2^3
		
		if(vl_cf(&pJacPointP->vY,&pvl[7]) < 0)//y < l7
			vl_add(&pJacPointP->vY,pvp);
		vl_usub(&pJacPointP->vY,&pvl[7]);//y -= l7
		//Деление на 2 mod p
		if(vl_test(&pJacPointP->vY,0))//нечетное
			vl_add(&pJacPointP->vY,pvp);
		vl_shr(&pJacPointP->vY);
		//if(vl_cf(&pJacPointP->vY,pvp)>=0)
		//  vl_usub(&pJacPointP->vY,pvp);
		//
	}  
	//
	return iRet;
}

// Удвоение точки на элл-й кривой 
int fEcc_DubJacP(          
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой (a*R)
          PTJacPointEcc pJacPointP,       //(IN/OUT)  Jacobian point P P=2*P
          TVlong pvl[9]){                 //(IN)  буфер объектов TVlong
  int iRet = 0;
  PTVlong pvp;
  //
  if(pMe_p && pva && pJacPointP){
    //
    iRet = 1;
    pvp = &(pMe_p->M);//mod p
    //
    if(pJacPointP->vZ.value.n == 0){//Pz=0 - точка в бесконечночти (0,0) P=2*0=0
      return iRet;
    }
    //удвоение в проективном представлении
    //
    //vl_copy(&pvl[3],pva);
    // l0 = 3*Px^2+a*Pz^4
	vl_copy(&pvl[1],&pJacPointP->vX);
	me_mul(&pvl[1],&pJacPointP->vX,pMe_p);//l1 = Px^2
    vl_copy(&pvl[0],&pvl[1]);
    vl_shl(&pvl[0]);//l0 = 2*Px^2
    if(vl_cf(&pvl[0],pvp)>=0)
      vl_usub(&pvl[0],pvp);
    vl_add(&pvl[0],&pvl[1]);//l0 = 3*Px^2
    if(vl_cf(&pvl[0],pvp)>=0)
      vl_usub(&pvl[0],pvp);
    //
	vl_copy(&pvl[1],&pJacPointP->vZ);
	me_mul(&pvl[1],&pJacPointP->vZ,pMe_p);//l1 = Pz^2
	me_mul(&pvl[1],&pvl[1],pMe_p);//l1 = Pz^4
	me_mul(&pvl[1],pva,pMe_p);//l1 = a*Pz^4
    //
    vl_add(&pvl[0],&pvl[1]);
    if(vl_cf(&pvl[0],pvp)>=0)
      vl_usub(&pvl[0],pvp);    //l0 = 3*Px^2 + a*Pz^4
      
    // l1 = 4*Px*Py^2
    vl_copy(&pvl[2],&pJacPointP->vY);
    vl_shl(&pvl[2]);
    if(vl_cf(&pvl[2],pvp)>=0)   //l2 = 2*Py
      vl_usub(&pvl[2],pvp);
    //
	me_mul(&pJacPointP->vZ,&pvl[2],pMe_p);// Pz_out = Pz*2*Py
    //
	me_mul(&pvl[2],&pvl[2],pMe_p);//l2 = l2^2 = 4*Py^2
    vl_copy(&pvl[1],&pvl[2]);     //l1 = 4*Py^2
	me_mul(&pvl[1],&pJacPointP->vX,pMe_p);//l1 = 4*Px*Py^2

    //Px_out = l0^2 - 2*l1
	vl_copy(&pJacPointP->vX,&pvl[0]);
	me_mul(&pJacPointP->vX,&pvl[0],pMe_p);//Px = l0^2
    
    if(vl_cf(&pJacPointP->vX,&pvl[1]) < 0)//x < l1
      vl_add(&pJacPointP->vX,pvp);
    vl_usub(&pJacPointP->vX,&pvl[1]);

    if(vl_cf(&pJacPointP->vX,&pvl[1]) < 0)//x < l1
      vl_add(&pJacPointP->vX,pvp);
    vl_usub(&pJacPointP->vX,&pvl[1]);//Px_out = l0^2 - 2*l1
    
    //l2 = 8*Py^4
	me_mul(&pvl[2],&pvl[2],pMe_p);//l2 = l2^2 = 16*Py^4
    // Деление на 2 mod p
    if(vl_test(&pvl[2],0))//нечетное
      vl_add(&pvl[2],pvp);
    vl_shr(&pvl[2]);
    //if(vl_cf(&pvl[2],pvp)>=0)
    //  vl_usub(&pvl[2],pvp);

    //Py_out = l0*(l1-Px_out)-l2
    if(vl_cf(&pvl[1],&pJacPointP->vX) < 0)//l1 < x
      vl_add(&pvl[1],pvp);
    vl_usub(&pvl[1],&pJacPointP->vX);//l1 = l1-Px_out
    
	vl_copy(&pJacPointP->vY,&pvl[0]);
	me_mul(&pJacPointP->vY,&pvl[1],pMe_p); //Py = l0*(l1-Px_out)
    
    if(vl_cf(&pJacPointP->vY,&pvl[2]) < 0)//y < l2
      vl_add(&pJacPointP->vY,pvp);
    vl_usub(&pJacPointP->vY,&pvl[2]);//Py_out = l0*(l1-Px_out)-l2
    //
  }
  //
  return iRet;
}

// d-кратная композиция точки на элл-й кривой [d]*P (P is normalise)
int fEcc_MulJacP_useXP(
          const PTMontExp pMe_p,        //(IN)  объект для выполнения быстрых операций a*b*R^-1 mod p
          const PTVlong pva,            //(IN)  коэффициент эллиптической кривой
          const PTVlong pvd,            //(IN)  [d]-длинное целое
          PTJacPointEcc* pJacPointP,    //(IN/OUT)  Jacobian point P P[0]=[d]*P[0];(const P[1]=3P; P[2]=7P;....)
          TVlong pvl_buff[9])           //(IN)  буфер объектов TVlong
{
		  int iRet = 0, i, k, ksize;
	TJacPointEcc jacPointR;
	//
	if(pMe_p && pva && pvd && pJacPointP[0])
	{
		//array [1P,3P,7P,15P,...0]
		for(ksize = 0; pJacPointP[ksize+1]; ++ksize);
		//
		//
		fEccJacPoint_init(&jacPointR);
		//
		i = vl_bits(pvd) - 1;
		//
		//iRet = 1;
		while(i >= 0)
		{
			//
			iRet = fEcc_DubJacP(pMe_p,pva,&jacPointR,pvl_buff);//R=2*R
			//
			if(iRet != 1)
				break;
			//
			if (!vl_test(pvd,i)){//i-s bit = 0
				i-=1;
				continue;
			}
			//
			k = 0;i-=1;
			while((i > 0)&&(k < ksize)){
				if (vl_test(pvd,i)){
					k+=1;i-=1;
					iRet = fEcc_DubJacP(pMe_p,pva,&jacPointR,pvl_buff);//R=2*R
					if(iRet != 1)
						break;
				}else{
					break;
				}
			}
			if(iRet != 1)
				break;
			//
			iRet = fEcc_AddJacP(pMe_p,pva,pJacPointP[k],&jacPointR,pvl_buff);//R = R + [k]P
			if(iRet != 1)
				break;
		}
		//
		if(iRet == 1)
			fEcc_Jac2Jac(pJacPointP[0],&jacPointR);
		//
		fEccJacPoint_delete(&jacPointR);
	}
	//
	return iRet;
}


// d-кратная композиция точки на элл-й кривой [d]*P
int fEcc_MulJacP(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          const PTVlong pvd,              //(IN)  [d]-длинное целое
          PTJacPointEcc pJacPointP){      //(IN/OUT)  Jacobian point P P=[d]*P
	int iRet = 0, i;
	TJacPointEcc jacPoint3P,jacPoint7P,jacPoint15P;
	TVlong pvl_buff[9];
	//
	if(pMe_p && pva && pvd && pJacPointP){
	//
	fEccJacPoint_init(&jacPoint3P);
	fEccJacPoint_init(&jacPoint7P);
	fEccJacPoint_init(&jacPoint15P);
	//
	for(i=0; i<9; i++)
		vl_init(&pvl_buff[i]);
	//
	//compute 3P, 7P,... points...
	iRet = fEcc_Jac2Jac(&jacPoint3P,pJacPointP);
	//
	if(iRet == 1)
		iRet = fEcc_DubJacP(pMe_p,pva,&jacPoint3P,pvl_buff);//2P
	//
	if(iRet == 1)
		iRet = fEcc_AddJacP(pMe_p,pva,pJacPointP,&jacPoint3P,pvl_buff);//3P
	//
	if(iRet == 1)
		iRet = fEcc_Jac2Jac(&jacPoint7P,&jacPoint3P);//3P
	//
	if(iRet == 1)
		iRet = fEcc_DubJacP(pMe_p,pva,&jacPoint7P,pvl_buff);//6P
	//
	if(iRet == 1)
		iRet = fEcc_AddJacP(pMe_p,pva,pJacPointP,&jacPoint7P,pvl_buff);//7P
	//
	if(iRet == 1)
		iRet = fEcc_Jac2Jac(&jacPoint15P,&jacPoint7P);//7P
	//
	if(iRet == 1)
		iRet = fEcc_DubJacP(pMe_p,pva,&jacPoint15P,pvl_buff);//14P
	//
	if(iRet == 1)
		iRet = fEcc_AddJacP(pMe_p,pva,pJacPointP,&jacPoint15P,pvl_buff);//15P
	//
	if(iRet == 1)
	{
		PTJacPointEcc pArrJacPoint[] = {pJacPointP,&jacPoint3P,&jacPoint7P,&jacPoint15P,0};
		iRet = fEcc_MulJacP_useXP(pMe_p,pva,pvd,pArrJacPoint,pvl_buff);
	}
	//
	for(i=0; i<9; i++)
		vl_delete(&pvl_buff[i]);
	//
	fEccJacPoint_delete(&jacPoint3P);
	fEccJacPoint_delete(&jacPoint7P);
	fEccJacPoint_delete(&jacPoint15P);
	}
	//
	return iRet;
}

// d,l-кратная композиция точек на элл-й кривой [l]*Q+[d]*P
int fEcc_MulJacPQ(
          const PTMontExp pMe_p, //(IN)  объект для выполнения быстрых операций a*b mod p
          const PTVlong pva,              //(IN)  коэффициент эллиптической кривой
          const PTVlong pvl,              //(IN)  [l]-длинное целое
          PTJacPointEcc pJacPointQ,       //(IN)  Jacobian point Q(will be normalise)
          const PTVlong pvd,              //(IN)  [d]-длинное целое
          PTJacPointEcc pJacPointP)      //(IN/OUT)  Jacobian point P P=[l]*Q+[d]*P
{
	int iRet = 0, i, bits;
	unsigned int bit_d, bit_l;
	TJacPointEcc jacPointPQ,jacPointR;
	TVlong pvl_buff[9];
	//
	if(pMe_p && pva && pvl&& pJacPointQ && pvd && pJacPointP){
	//
	for(i=0; i<9; i++)
		vl_init(&pvl_buff[i]);
	//
	fEccJacPoint_init(&jacPointPQ);
	//
	iRet = fEcc_Jac2Jac(&jacPointPQ,pJacPointP);//PQ = P
	//
	if(iRet == 1)
		iRet = fEcc_AddJacP(pMe_p,pva,pJacPointQ,&jacPointPQ,pvl_buff);//PQ = P + Q
	//
	if(iRet != 1){//jacPointPQ = P+Q
		fEccJacPoint_delete(&jacPointPQ);
		//
		for(i=0; i<9; i++)
			vl_delete(&pvl_buff[i]);
		//
		return 0;
	}
	//
	fEccJacPoint_init(&jacPointR);
	//
	bits = vl_bits(pvd);
	i = vl_bits(pvl);
	if(i > bits)
	  bits = i;
	//
	if(bits >0)
	for(i = bits-1; i >= 0; --i){
		//
		iRet = fEcc_DubJacP(pMe_p,pva,&jacPointR,pvl_buff);//R=2*R
		if(iRet != 1)
			break;
		//
		bit_d = vl_test(pvd,i);
		bit_l = vl_test(pvl,i);
		//
		if(bit_d){
			if(bit_l)//R = R + (Q+P)
				iRet = fEcc_AddJacP(pMe_p,pva,&jacPointPQ,&jacPointR,pvl_buff);
			else//R = R + (0+P)
				iRet = fEcc_AddJacP(pMe_p,pva,pJacPointP,&jacPointR,pvl_buff);
		}else{
			if(bit_l)//R = R + (Q+0)
				iRet = fEcc_AddJacP(pMe_p,pva,pJacPointQ,&jacPointR,pvl_buff);
		}
		//
		if(iRet != 1)
			break;
	}
	//
	if(iRet == 1)
		fEcc_Jac2Jac(pJacPointP,&jacPointR);
	//
	fEccJacPoint_delete(&jacPointR);
	fEccJacPoint_delete(&jacPointPQ);
	//
	for(i=0; i<9; i++)
	  vl_delete(&pvl_buff[i]);
	//
	}
	return iRet;
}
