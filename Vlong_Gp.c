#include "Vlong_Gp.h"
//test
#include <stdio.h>
void print_vlong(PTVlong pvl);
//test

typedef union unionu_tvlong{
	unsigned long long udub;
	TVLongUnit pu[2];
} unionu_tvlong;


int vl_gcd(PTVlong pV, const PTVlong pVx, const PTVlong pVy){//pV = gcd(pVx,pVy)
	int iRes = 0;
	TVlong x,y,r;

	if(pV && pVx && pVy){
		iRes = 1;
		//
		vl_init(&x);vl_copy(&x,pVx);
		vl_init(&y);vl_copy(&y,pVy);
		vl_init(&r);
		//
		while(1){
			if(y.value.n == 0){
				vl_copy(pV,&x);
				break;
			}
			vl_dive(0,&r,&x,&y);//r =  x%y
			vl_copy(&x,&r);//x = x%y
			
			if(x.value.n == 0){
				vl_copy(pV,&y);
				break;
			}
			vl_dive(0,&r,&y,&x);//r = y%x
			vl_copy(&y,&r);//y = y%x
		}
		//
		vl_delete(&x);
		vl_delete(&y);
		vl_delete(&r);
	}
	return iRes;
}

void vl_modinv_all(PTVlong pV, const PTVlong pVa, const PTVlong pVm){// modular inverse
	// returns i in range 1..m-1 such that i*a = 1 mod m
	// a must be in range 1..m-1
	TVlong j,i,b,c,x,y,t;
	//
	if(pV && pVa && pVm){
		vl_init(&j);vl_set(&j,0,1); // j = 1
		vl_init(&i);vl_set(&i,0,0); // i = 0
		vl_init(&b);vl_copy(&b,pVm);// b = m
		vl_init(&c);vl_copy(&c,pVa);// c = a
		vl_init(&x);
		vl_init(&y);
		vl_init(&t);
		//
		while(c.value.n != 0){
			vl_dive(&x,&y,&b,&c);//x = b / c; y = b - x*c;
			vl_copy(&b,&c);// b=c;
			vl_copy(&c,&y);// c=y
			vl_copy(&y,&j);// y=j
			
			vl_mule(&t,&j,&x);
			vl_sub(&i,&t);
			vl_copy(&j,&i);// j = i - j*x
			
			vl_copy(&i,&y);// i=y
		}
		//
		if(vl_is_negative(&i))
			vl_add(&i,pVm);
		//
		vl_copy(pV,&i);//pV = i
		//
		vl_delete(&j);
		vl_delete(&i);
		vl_delete(&b);
		vl_delete(&c);
		vl_delete(&x);
		vl_delete(&y);
		vl_delete(&t);
	}
}

void vl_modinv(PTVlong pV, const PTVlong pVa, const PTVlong pVm){ // modular inverse
	// m must be a prime(!)
	// returns pV in range 1..m-1 such that pV*a = 1 mod m
	// a must be in range 1..m-1
	TVlong u,v1,x1,x2;
	//
	if(pV && pVa && pVm){
		//
		vl_init(&v1);vl_set(&v1,0,1);  //v1 = 1
		vl_init(&x1);vl_set(&x1,0,1);vl_reserve(&x1,pVm->value.n + 1);   //x1 = 1
		vl_init(&x2);vl_set(&x2,0,0);vl_reserve(&x2,pVm->value.n + 1);   //x2 = 0
		vl_init(&u);vl_copy(&u,pVa);    // u = a
		vl_copy(pV,pVm);    // pV = m
		//
		while((vl_fast_compare(&u,&v1) == 0)&&(vl_fast_compare(pV,&v1) == 0)){//while(u!=1 && pV!=1)
			//
			while(!vl_test(&u,0))//while u is even (четное)
			{
				vl_shr(&u);//u=u/2
				//x1 = x1/2 mod m
				if(vl_test(&x1,0))//(нечетное)
					vl_add(&x1,pVm);//x1 = x1+m
				//
				vl_shr(&x1);//x1=x1/2
			}
			//
			while(!vl_test(pV,0))//while pV is even (четное)
			{
				vl_shr(pV);//pV=pV/2
				//x2 = x2/2 mod m
				if(vl_test(&x2,0))//(нечетное)
					vl_add(&x2,pVm);//x2 = x2+m
				//
				vl_shr(&x2);//x2=x2/2
			}
			//
			if(vl_cf(&u,pV) >= 0){//(u >= pV)
				vl_usub(&u,pV);//u-=pV
				if(vl_cf(&x1,&x2) < 0)//x1 < x2
					vl_add(&x1,pVm);
				vl_usub(&x1,&x2);//x1-=x2 mod m
			}else{
				vl_usub(pV,&u);//pV-=u
				if(vl_cf(&x2,&x1) < 0)//x2 < x1
					vl_add(&x2,pVm);
				vl_usub(&x2,&x1);//x2-=x1 mod m
			}
		}
		//
		if(vl_fast_compare(&u,&v1)){//u==1
			vl_copy(pV,&x1);//pV=x1
		}else{
			vl_copy(pV,&x2);//pV=x2
		}
		//
		while(vl_cf(pV,pVm) >= 0)//pV >= m
			vl_usub(pV,pVm);//pV-=m
		//
		while(vl_is_negative(pV))
			vl_add(pV,pVm);//pV+=m
		//
		vl_delete(&x2);
		vl_delete(&x1);
		vl_delete(&v1);
		vl_delete(&u);
	}
}

//--- BarrettReduction ---
void br_init(PTBarrettReduction pBr, const PTVlong pM){//init pBr with module pM
	if(pBr && pM){
		//if((pBr->uiN == 0) || (!vl_fast_compare(pM,&pBr->M)))
		{//check - may be init already
			//
			vl_init(&pBr->Q);
			vl_init(&pBr->T);vl_set(&pBr->T,0,1);  //T = 1
			vl_init(&pBr->res);
			vl_inite(&pBr->M,pM->value.n);vl_copy(&pBr->M,pM);
			vl_init(&pBr->mu);
			//
			pBr->uiN = vl_bits(pM);// + 1 !!!
			//T = 2^2n
			vl_shlx(&pBr->T,2*pBr->uiN);
			//mu = 2^2n / M
			vl_dive(&pBr->mu,0,&pBr->T,pM);
			pBr->uiMu = vl_bits(&pBr->mu);
			//
			//vl_clear(&t);
			//vl_set(&t,0,1);
			//vl_shlx(&t,pBr->uiN);
			//if(vl_cf(&pBr->mu,&t) >= 0)
			//  vl_usub(&pBr->mu,&t);
			//
			//vl_reserve(&pBr->Q,pBr->M.value.n+pBr->mu.value.n);
			//vl_reserve(&pBr->T,pBr->M.value.n+pBr->mu.value.n);
		}
	}
}

void br_init_uchar_BE(PTBarrettReduction pBr, 
                      const unsigned char *psM, unsigned int lenM,
                      const unsigned char *psMu, unsigned int lenMu){
	//
	if(pBr && psM && psMu){
		vl_init(&pBr->M);
		vl_init(&pBr->mu);
		vl_init(&pBr->Q);
		vl_init(&pBr->T);
		vl_init(&pBr->res);

		//
		vl_set_uchar_BE(&pBr->M,psM,lenM);
		vl_set_uchar_BE(&pBr->mu,psMu,lenMu);
		pBr->uiN = vl_bits(&pBr->M);
		pBr->uiMu = vl_bits(&pBr->mu);
		//
		//vl_reserve(&pBr->Q,pBr->M.value.n+pBr->mu.value.n);
		//vl_reserve(&pBr->T,pBr->M.value.n+pBr->mu.value.n);
	}
}

void br_init_copy(PTBarrettReduction pBr, const PTBarrettReduction pBrIni){//init pBr with pBrIni (copy)
	if(pBr && pBrIni && (pBr != pBrIni)){
		vl_inite(&pBr->M,pBrIni->M.value.n);
		vl_inite(&pBr->mu,pBrIni->mu.value.n);
		vl_inite(&pBr->Q,pBrIni->M.value.n+pBrIni->mu.value.n);
		vl_inite(&pBr->T,pBrIni->M.value.n+pBrIni->mu.value.n);
		vl_inite(&pBr->res,2*pBrIni->M.value.n);
		//
		vl_copy(&pBr->M,&pBrIni->M);
		vl_copy(&pBr->mu,&pBrIni->mu);
		pBr->uiN = pBrIni->uiN;
		pBr->uiMu = pBrIni->uiMu;
		//
		//vl_reserve(&pBr->Q,pBr->M.value.n+pBr->mu.value.n);
		//vl_reserve(&pBr->T,pBr->M.value.n+pBr->mu.value.n);
	}
}

void br_delete(PTBarrettReduction pBr){
	if(pBr){
		vl_delete(&pBr->M);
		vl_delete(&pBr->mu);
		vl_delete(&pBr->Q);
		vl_delete(&pBr->T);
		vl_delete(&pBr->res);
		pBr->uiN = pBr->uiMu = 0;
	}
}

void br_module(PTVlong pV, const PTBarrettReduction pBr){ //pV = pV mod pBr->M
	//
	unsigned int uiNV;
	//
	if(pV && pBr && (vl_cf(pV,&pBr->M) >= 0)){ //pV >= M
		//1 T = pV / 2^n
		vl_copy(&pBr->T,pV);
		vl_shrx(&pBr->T,pBr->uiN);//+1
		//2 Q = T*mu
		vl_fast_mule(&pBr->Q,&pBr->T,&pBr->mu,pBr->uiN + pBr->uiMu);
		//3 Q = [[pV / 2^n] * mu / 2^n]
		vl_shrx(&pBr->Q,pBr->uiN);//-1
		//4.1 pV = pV mod 2^(n+bpu)
		uiNV = pBr->M.value.n + 1;
		if(pV->value.n > uiNV){
			pV->value.n = uiNV;//(pBr->uiN+BPU) bits
			vl_normal(pV);
		}
		//4.2 T = Q*M mod 2^(n+bpu)
		vl_fast_mule(&pBr->T,&pBr->Q,&pBr->M,pBr->uiN + BPU);
		//5 pV = pV - Q*M mod 2^(n+bpu)
		if(vl_cf(pV,&pBr->T) < 0)//(pV < T)
		{//pV == 0 - 99.9% or 100% ?
			//pV = pV + 2^(n+bpu)
			vl_clear(&pBr->Q);vl_set(&pBr->Q,0,1);vl_shlx(&pBr->Q,pBr->uiN+BPU);
			vl_add(pV,&pBr->Q);
			while(vl_cf(pV,&pBr->T) < 0)//pV < T
				vl_add(pV,&pBr->M);// pV += M
		}
		vl_usub(pV,&pBr->T);
		//6
		while(vl_cf(pV,&pBr->M) >= 0)//(pV >= M)
			vl_usub(pV,&pBr->M);// pV = pV - M
	}
}

/*
void br_module(PTVlong pV, const PTBarrettReduction pBr){ //pV = pV mod pBr->M
	//
	unsigned int uiNV;
	//
	if(pV && pBr && (vl_cf(pV,&pBr->M) >= 0)){ //pV >= M
		//1 T = pV / 2^n
		vl_copy(&pBr->T,pV);
		vl_shrx(&pBr->T,pBr->uiN);//+1
		//2 Q = T*mu
		vl_fast_mule(&pBr->Q,&pBr->T,&pBr->mu,pBr->uiN + pBr->uiMu);
		//3 Q = [[pV / 2^n] * mu / 2^n]
		vl_shrx(&pBr->Q,pBr->uiN);//-1
		//4.1 pV = pV mod 2^(n+bpu)
		uiNV = pBr->M.value.n + 1;
		if(pV->value.n > uiNV){
			pV->value.n = uiNV;//(pBr->uiN+BPU) bits
			vl_normal(pV);
		}
		//4.2 T = Q*M mod 2^(n+bpu)
		vl_fast_mule(&pBr->T,&pBr->Q,&pBr->M,pBr->uiN + BPU);
		//5 pV = pV - Q*M mod 2^(n+bpu)
		if(vl_cf(pV,&pBr->T) < 0)//(pV < T)
		{//pV == 0 - 99.9% or 100% ?
			//pV = pV + 2^(n+bpu)
			vl_clear(&pBr->Q);vl_set(&pBr->Q,0,1);vl_shlx(&pBr->Q,pBr->uiN+BPU);
			vl_add(pV,&pBr->Q);
			while(vl_cf(pV,&pBr->T) < 0)//pV < T
				vl_add(pV,&pBr->M);// pV += M
		}
		vl_usub(pV,&pBr->T);
		//6
		while(vl_cf(pV,&pBr->M) >= 0)//(pV >= M)
			vl_usub(pV,&pBr->M);// pV = pV - M
	}
}
//*/

void br_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTBarrettReduction pBr){
	if(pV && pVx && pVy && pBr){
		vl_fast_mule(&pBr->res,pVx,pVy,2*pBr->uiN);
		br_module(&pBr->res,pBr);
		vl_copy(pV,&pBr->res);
	}
}

//--- end BarrettReduction ---

void me_mont(PTVlong pV, const PTMontExp pMe);

static TVLongUnit me_mod_inv_Hensel(TVLongUnit m)//return m^-1 mod 2^BPU
{// 1)ret*m = 1 mod 2^0 ..... BPU)ret*m = 1 mod 2^BPU
	unsigned int i;
	TVLongUnit ret = 1, s, mask = 1;//m1
	for(i = 1; i < BPU; ++i){
		//m1 = m mod 2^(i+1)
		mask = (mask<<1) | 1;
		//m1 = m & mask;
		s = (m*ret) & mask;
		//for(j = 0,s = 0; j < i; ++j){
		//	s = (s<<1);
		//	if((ret>>(i-j-1)) & 0x01)
		//		s+= m1;
		//	//
		//	s &= mask;
		//}
		//
		if(s != 1)
			ret |= (TVLongUnit)1<<i;
		//
	}		
	return ret;
}
//--- Montgomery modular exponentiation and multiplication ---
void me_init(PTMontExp pMe, const PTVlong pM){   //init pMe with module pM
	if(pMe && pM){
		//if((pMe->uiN == 0) || (!vl_fast_compare(pM,&pMe->M))){//check - may be init already
		{
			vl_init(&pMe->M);vl_copy(&pMe->M,pM);
			//vl_init(&pMe->n1);
			vl_inite(&pMe->r2,pMe->M.value.n+1);
			//
			vl_inite(&pMe->T,2*pMe->M.value.n+1);
			//
			pMe->uiN = vl_bits(pM);
			//
			//n1short calculate
			pMe->n1short = 0 - me_mod_inv_Hensel(vl_get(pM,0));// 2^BPU - M^-1 mod 2^BPU
			//printf("\nme_init:n2s^= %08X",pMe->n1short);
			//n1short calculate

			/*
			//n1 = R - modinv_all(M,R);
			vl_clear(&R);vl_set(&R,0,1);vl_shlx(&R,pMe->uiN);//R = 2^bits(M)
			vl_copy(&pMe->r2,&R);
			vl_modinv_all(&pMe->n1,pM,&R);//n1 = modinv_all(M,R)
			//
			//vl_clear(&R);vl_set(&R,0,1);vl_shlx(&R,pMe->uiN);//R = 2^bits(M)
			
			vl_usub(&R,&pMe->n1);//R = R-n1
			vl_copy(&pMe->n1,&R);
			
			printf("\nme_init: n1 = ");
			print_vlong(&pMe->n1);
			//*/
			
			//r2 = R^2 mod M
			//vl_clear(&pMe->T);vl_set(&pMe->T,0,1);vl_shlx(&pMe->T,2*pMe->uiN);//T = 2^2n
			//vl_dive(0,&pMe->r2,&pMe->T,pM);//r2 = R^2 mod M

			vl_inite(&pMe->K,2*pMe->M.value.n + 1);
			
			//r2 = R^2 mod M (4 TVlong economy!)
			PTVlong pvT = &pMe->T;
			PTVlong pvK = &pMe->K;
			vl_set(pvK,0,1);vl_shlx(pvK,2*pMe->uiN);//pvK = 2^2n
			do{
				vl_copy(pvT,pM);
				vl_shlx(pvT,vl_bits(pvK) - pMe->uiN);
				//
				if(vl_cf(pvT,pvK) > 0)
					vl_shr(pvT);
				//
				vl_usub(pvK,pvT);
				//
			}while(vl_cf(pvK,pM) >= 0);//(pvK >= M)
			vl_copy(&pMe->r2,pvK);
			//
			//
			//vl_clear(&pMe->T);vl_set(&pMe->T,0,1);
			//me_mont(&pMe->T,pMe);
			//me_mont(&pMe->T,pMe);//T = R^-2 mod M
			//vl_modinv(&pMe->r2,&pMe->T,&pMe->M);
		}
	}
}

void me_delete(PTMontExp pMe){
	if(pMe){
		vl_delete(&pMe->M);
		//vl_delete(&pMe->n1);
		vl_delete(&pMe->r2);
		vl_delete(&pMe->T);
		vl_delete(&pMe->K);
	}
}

/*
static void me_mul(PTVlong pV, const PTVlong pVy, const PTMontExp pMe){
	if(pV && pVy && pMe){
		// T = pV*pVy;
		vl_fast_mule(&pMe->T,pV,pVy,2*pMe->uiN);
		// K = ( T * n1 ) % R;
		vl_fast_mule(&pMe->K,&pMe->n1,&pMe->T,pMe->uiN);
		// pV = ( T + K*M ) / R;
		//vl_fast_mule(pV,&pMe->K,&pMe->M,2*pMe->uiN);
		//vl_add(pV,&pMe->T);
		//vl_shrx(pV,pMe->uiN);
		vl_fast_add_mule(&pMe->T,&pMe->K,&pMe->M,2*pMe->uiN +1);
		vl_shrx(&pMe->T,pMe->uiN);
		vl_copy(pV,&pMe->T);
		//
		while(vl_cf(pV,&pMe->M) >= 0)//(pV >= M)
			vl_usub(pV,&pMe->M);// pV = pV - M
	}
}
/*/
//pV = pV*pVy*R^-1 mod M
void me_mul(PTVlong pV, const PTVlong pVy, const PTMontExp pMe){
	if(pV && pVy && pMe){
		// T = pV*pVy;
		vl_fast_mule(&pMe->T,pV,pVy,2*pMe->uiN);
		//
		me_mont(&pMe->T,pMe);
		//
		vl_copy(pV,&pMe->T);
	}
}
//*/

/*
static void me_mont(PTVlong pV, const PTMontExp pMe){
	if(pV && pMe){
		// T = pV;
		//vl_copy(&pMe->T,pV);
		// K = ( T * n1 ) % R;
		vl_fast_mule(&pMe->K,pV,&pMe->n1,pMe->uiN);
		// pV = ( T + K*M ) / R;
		vl_fast_add_mule(pV,&pMe->K,&pMe->M,2*pMe->uiN +1);
		//show
		printf("\nme_mont: ( T + K*M ) = ");
		print_vlong(pV);
		//show
		vl_shrx(pV,pMe->uiN);
		//
		while(vl_cf(pV,&pMe->M) >= 0)//(pV >= M)
			vl_usub(pV,&pMe->M);// pV = pV - M
	}
}
/*/
//pV = pV*R^-1 mod M
void me_mont(PTVlong pV, const PTMontExp pMe){
	PTVlong pVres;
	TVLongUnit *paV, *paU, *paUl, *paM;
	TVLongUnit m,c,n1short;
	unsigned int i,j,Nm;
	unionu_tvlong ww;
	
	if(pV && pMe){
		pVres = &pMe->K;
		vl_clear(pVres);
		vl_copy(pVres,pV);//pVres = pV
		paU = pVres->value.pa;
		paM = pMe->M.value.pa;
		paV = pV->value.pa;
		//
		Nm = pMe->M.value.n;
		n1short = pMe->n1short;
		//
		for(i = 0; i < Nm; ++i, ++paU){
			//
			c = 0;
			m = (*paU)*n1short;//mod 2^BPU
			//
			for(j = 0, paUl = paU; j < Nm; ++j, ++paUl){
				ww.udub = (unsigned long long)m;
				ww.udub *= (unsigned long long)paM[j];
				//
				ww.udub += (unsigned long long)c;
				ww.udub += (unsigned long long)*paUl;
				//
				*paUl = ww.pu[0];
				c = ww.pu[1];
			}
			//continue carrying
			while(c){// && (j < Nm+Nm-i)){
				ww.udub = (unsigned long long)c;
				ww.udub += (unsigned long long)*paUl;
				//
				*paUl = ww.pu[0];
				c = ww.pu[1];
				//
				++paUl;
			}
		}
		//show
		//pVres->value.n = 2*Nm + 1;
		//vl_normal(pVres);
		//printf("\nme_mont: ( T + K*M ) = ");
		//print_vlong(pVres);
		//show
		paUl = pVres->value.pa;
		for(i = 0; i < Nm+1; ++i){//paU = &paU[Nm]
			paUl[i] = paU[i];
		}
		while (i && paUl[i-1]==0) i-=1; // normalise
		pVres->value.n = i;
		//
		while(vl_cf(pVres,&pMe->M) >= 0)//(pVres >= M)
			vl_usub(pVres,&pMe->M);// pVres = pVres - M
		vl_copy(pV,pVres);//pV = pVres
	}
}
//*/

//pV = pVx*pVy mod pMe->M       WARNINNG pV != pVy
void me_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTMontExp pMe){
	//
	if(pV && pVx && pVy && pMe){
		if(pV != pVy){
			vl_copy(pV,pVx);
			me_mul(pV,pVy,pMe);//pV = pVx*pVy*R^-1 mod M
		}else{
			//pV == pVy
			me_mul(pV,pVx,pMe);//pV = pVx*pVy*R^-1 mod M
		}
		//
		me_mul(pV,&pMe->r2,pMe);//pV = (pVx*pVy*R^-1)*(R^2)*R^-1  mod M = pVx*pVy mod M
	}
}

void me_module(PTVlong pV, const PTMontExp pMe){
	//
	if(pV && pMe){
		if(pV->value.n >= pMe->M.value.n){
			me_mont(pV,pMe);
			//
			me_mul(pV,&pMe->r2,pMe);//pV = (pV*R^-1)*(R^2)*R^-1  mod M = pVx*pVy mod M
		}
	}
}

//*
//pV = pVx^pVe mod pMe->M
void me_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTMontExp pMe){
	int i;
	TVlong v,v3;
	//
	if(pV && pVx && pVe && pMe){
		vl_init(&v);
		vl_init(&v3);
		//
		vl_copy(&v,pVx);
		me_mul(&v,&pMe->r2,pMe);//v = pVx*(R^2)*R^-1  mod M = pVx*R mod M
		//
		vl_copy(&v3,&v);
		me_mul(&v3,&v,pMe);//v3 = pVx^2*R mod M
		me_mul(&v3,&v,pMe);//v3 = pVx^3*R mod M
		//pV = R mod M
		vl_copy(pV,&pMe->r2);
		me_mont(pV,pMe);//pV = (R^2)*R^-1  mod M = R mod M
		//
		i = vl_bits(pVe) - 1;
		while (1)
		{
			if(i < 0) break;
			//
			me_mul(pV,pV,pMe);
			//
			if (!vl_test(pVe,i)){//i-s bit = 0
				i-=1;
				continue;
			}
			//i-s bit = 1
			if(i > 0){
				if (vl_test(pVe,i-1)){//(i-1)-s bit = 1
					me_mul(pV,pV,pMe);
					me_mul(pV,&v3,pMe);
					i-=2;
					continue;
				}
			}
			//i-s bit = 1
			me_mul(pV,&v,pMe);
			i-=1;
		}
		//
		me_mont(pV,pMe);
		//
		vl_delete(&v);
		vl_delete(&v3);
	}
}
//*/
//pV = 2^pVe mod pMe->M  //M length must be >= 256 bits
void me_2modexp(PTVlong pV, const PTVlong pVe, const PTMontExp pMe){
	int i;
	unsigned int ex;
	TVlong v;
	const TVLongUnit *paE;
	//
	if(pV && pVe && pMe){
		vl_init(&v);
		//
		//pV = R mod M
		vl_copy(pV,&pMe->r2);
		me_mont(pV,pMe);//pV = (R^2)*R^-1  mod M = R mod M
		//
		for(paE = pVe->value.pa + pVe->value.n - 1; paE >= pVe->value.pa; --paE){
			for(i = sizeof(TVLongUnit)-1; i >= 0; --i){
				ex = ((*paE) >> (i*8))&0xFF;
				//
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				me_mul(pV,pV,pMe);
				//
				if(ex){
					vl_copy(&v,&pMe->r2);
					vl_shlx(&v,ex);//v = R^2 mod M *(2^ex)
					me_mont(&v,pMe);//v = R*2^ex mod M
					me_mul(pV,&v,pMe);//
				}
			}
		}
		//
		me_mont(pV,pMe);
		//
		vl_delete(&v);
	}
}


/*
//pV = pVx^pVe mod pMe->M //simple release
void me_modexp(PTVlong pV, const PTVlong pVx, const PTVlong pVe, const PTMontExp pMe){
	unsigned int bits,i;
	TVlong t;
	//
	if(pV && pVx && pVe && pMe){
		vl_init(&t);
		//
		vl_copy(&t,pVx);
		me_mul(&t,&pMe->r2,pMe);//t = pVx*(R^2)*R^-1  mod M = pVx*R mod M
		//pV = R mod M
		vl_copy(pV,&pMe->r2);
		me_mont(pV,pMe);//pV = (R^2)*R^-1  mod M = R mod M
		//
		i = 0;
		bits = vl_bits(pVe);
		while (1)
		{
			if (vl_test(pVe,i))
				me_mul(pV,&t,pMe);
			i += 1;
			if ( i == bits ) break;
			me_mul(&t,&t,pMe);
		}
		//
		me_mont(pV,pMe);
		//
		vl_delete(&t);
	}
}
//*/
//--- end Montgomery modular exponentiation ---

//------------------------------------------------------------------------------

//--- Quisquater reduction ---
//init pQr with module pM
void qr_init(PTQuisquaterReduction pQr, const PTVlong pM)
{
	if(pQr && pM){
		pQr->uiN = vl_bits(pM);
		vl_init(&pQr->M);vl_copy(&pQr->M,pM);
		//vl_init(&pQr->del);
		vl_inite(&pQr->N,pM->value.n + 2);
		vl_inite(&pQr->Nsh,pM->value.n + 2);
		vl_inite(&pQr->res,pM->value.n + 3);
		//
		//vl_init(&pQr->q);
		//
		vl_set(&pQr->N,0,1);vl_shlx(&pQr->N,pQr->uiN + BPU);
		vl_copy(&pQr->Nsh,&pQr->N);//Nsh = 2^(n+bpu)
		vl_dive(&pQr->res,0,&pQr->N,pM);//res = 2^(n+bpu) / M
		//
		//vl_copy(&pQr->del,&pQr->res);
		//printf("\n 2^(n+bpu) / M = ");
		//print_vlong(&pQr->res);
		//
		vl_clear(&pQr->N);
		vl_fast_mule(&pQr->N,&pQr->res,pM,pQr->uiN + BPU);
		vl_usub(&pQr->Nsh,&pQr->N);//Nsh = 2^(n+bpu) - N
	}
}
//delete pQr
void qr_delete(PTQuisquaterReduction pQr)
{
	if(pQr){
		vl_delete(&pQr->M);
		vl_delete(&pQr->N);
		//vl_delete(&pQr->del);
		vl_delete(&pQr->Nsh);
		//vl_delete(&pQr->q);
		vl_delete(&pQr->res);
		pQr->uiN = 0;
	}
}

void print_vlong(PTVlong pvl);
//pV = (pV mod pQr->N) mod pQr->M
void qr_module(PTVlong pV, const PTQuisquaterReduction pQr)
{
	PTVlong pVres;
	if(pV && pQr && (vl_cf(pV,&pQr->M) >= 0)){ //pV >= M
		pVres = &pQr->res;
		//
		while(vl_cf(pV,&pQr->M) >= 0)//(pV >= M)
		{
			vl_copy(pVres,&pQr->M);
			vl_shlx(pVres,vl_bits(pV) - pQr->uiN);
			//
			if(vl_cf(pVres,pV) > 0)
				vl_shr(pVres);
			//
			vl_usub(pV,pVres);
			//printf("\nqr_module: V mod M = ");
			//print_vlong(pV);
		}
	}
}
//
/*
//pV = pVx*pVy mod pQr->N(!) pVx,pVy may len = len(M) + BPU
void qr_modmul_(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTQuisquaterReduction pQr)
{
	TVLongUnit *paX;
	PTVlong pVres;
	if(pV && pVx && pVy && pQr){
		pVres = &pQr->res;
		vl_clear(pVres);
		for(paX = pVx->value.pa + pVx->value.n - 1; paX >= pVx->value.pa; --paX){
			//
			vl_clear(&pQr->q);
			vl_set(&pQr->q,0,*paX);
			//1 ----------------
			vl_shlx(pVres,BPU);//
			vl_fast_add_mule(pVres,&pQr->q,pVy,pQr->uiN+BPU+BPU+BPU);
			//show
			//printf("\nqr_modmul:1)pVres= ");
			//print_vlong(pVres);
			//show
			//2 ----------------
			//vl_copy(&pQr->q,pVres);
			//vl_shrx(&pQr->q,pQr->uiN+BPU);
			vl_set(&pQr->q,0,vl_get(pVres,pQr->M.value.n + 1));
			vl_set(&pQr->q,1,vl_get(pVres,pQr->M.value.n + 2));
			//show
			//printf("\nqr_modmul:2) q   = ");
			//print_vlong(&pQr->q);
			//show
			//3 ----------------
			if(pVres->value.n > pQr->M.value.n + 1){
				pVres->value.n = pQr->M.value.n + 1;
				vl_normal(pVres);
			}
			//
			vl_fast_add_mule(pVres,&pQr->q,&pQr->Nsh,pQr->uiN+BPU+BPU);//+2 - very impotant!
			//
			//printf("\nqr_modmul:3)pVres= ");
			//print_vlong(pVres);
			//
			while(vl_cf(pVres,&pQr->N) >= 0)//(pV >= N)
				vl_usub(pVres,&pQr->N);// pV = pV - N
			//
		}
		//
		//qr_module(pVres,pQr);
		vl_copy(pV,pVres);
	}
}
//*/

//pV = pVx*pVy mod pQr->N(!) pVx,pVy may len = len(M) + BPU
void qr_modmul(PTVlong pV, const PTVlong pVx, const PTVlong pVy, const PTQuisquaterReduction pQr)
{
	PTVlong pVres;
	TVLongUnit *paX, *paY, *paU, *paNsh;// *paUl,
	TVLongUnit Qh,Ql,Run,c,c1;
	unsigned int j,N,Nm;
	unionu_tvlong ww;
	//paY[i] = B[i]; paX[i] = A[i];
	if(pV && pVx && pVy && pQr){
		pVres = &pQr->res;
		vl_clear(pVres);
		paU = pVres->value.pa;
		paY = pVy->value.pa;
		paNsh = pQr->Nsh.value.pa;
		Nm = pQr->M.value.n;
		//
		for(paX = pVx->value.pa + pVx->value.n - 1; paX >= pVx->value.pa; --paX)
		{
			// 1 ---
			// pVres = (pVres << BPU) + (*paX)*pVy
			Run = c = 0;
			N = pVy->value.n;
			for(j = 0; j < N; ++j){
				ww.udub = (unsigned long long)paY[j];
				ww.udub *= (unsigned long long)*paX;
				//
				ww.udub += (unsigned long long)c;
				ww.udub += (unsigned long long)Run;
				//
				Run = paU[j];
				paU[j] = ww.pu[0];
				c = ww.pu[1];
			}
			//
			N = Nm + 3;
			//while((c || Run)&&(j < N)){//max size = n+3 (c || Ru)
			while(j < N){//max size = n+3 // здесь по-любому надо сдвигать до упора!
				ww.udub = (unsigned long long)c;
				ww.udub += (unsigned long long)Run;
				Run = paU[j];
				paU[j] = ww.pu[0];
				c = ww.pu[1];
				++j;
			}
			//
			//show
			//printf("\nqr_modmul:1)pVres= ");
			//pVres->value.n = pQr->M.value.n + 3;
			//vl_normal(pVres);
			//print_vlong(pVres);
			//show
			// 2 --- 
			Ql = paU[Nm + 1];
			Qh = paU[Nm + 2];
			//show
			//printf("\nqr_modmul:2)q   = %08x %08X",Qh,Ql);
			//if(Qh >= 2)
			//	printf("\n!!!qr_modmul:Qh = %8X",Qh);
			//show
			// 3 --- 
			//pVres = pVres mod 2^(n+bpu) + (Qh,Ql)*Nsh
			//vl_fast_add_mule(pVres,&pQr->q,&pQr->Nsh,pQr->uiN+BPU+BPU);//+2 - very impotant!
			if(Ql || Qh)
			{
				paU[Nm + 1] = paU[Nm + 2] = 0;//pVres = pVres mod 2^(n+bpu)
				//
				c = c1 = 0;
				N = pQr->Nsh.value.n;
				if(Qh){
					for(j = 0; j < N; ++j){
						Run = paNsh[j];
						ww.udub = (unsigned long long)Run;
						ww.udub *= (unsigned long long)Ql;
						//
						//ww.udub += (unsigned long long)c;
						ww.udub += (unsigned long long)paU[j];
						//
						paU[j] = ww.pu[0];
						c = ww.pu[1];
						//Qh = 1 or 2
						ww.udub = (unsigned long long)Run;
						ww.udub += (unsigned long long)paU[j+1];
						ww.udub += (unsigned long long)c;
						ww.udub += (unsigned long long)c1;
						if(2 == Qh)
							ww.udub += (unsigned long long)Run;
						//
						paU[j+1] = ww.pu[0];
						c1 = ww.pu[1];
						//c = 0;
					}
					//continue carrying
					while((c1)&&(j < Nm + 2)){// c1!=0   +1
						++j;
						ww.udub = (unsigned long long)c1;
						ww.udub += (unsigned long long)paU[j];
						paU[j] = ww.pu[0];
						c1 = ww.pu[1];
					}
				}else{
					for(j = 0; j < N; ++j){
						Run = paNsh[j];
						ww.udub = (unsigned long long)Run;
						ww.udub *= (unsigned long long)Ql;
						//
						ww.udub += (unsigned long long)c;
						ww.udub += (unsigned long long)paU[j];
						//
						paU[j] = ww.pu[0];
						c = ww.pu[1];
					}
					//continue carrying
					while((c)&&(j < Nm + 2)){// c!=0
						ww.udub = (unsigned long long)c;
						ww.udub += (unsigned long long)paU[j];
						paU[j] = ww.pu[0];
						c = ww.pu[1];
						++j;
					}
					
				}
				//continue carrying
				// if c!=0 then c1=0
				// if c1!=0 then c=0
				//if(c1)
				//	printf("\n!!!qr_modmul:c1 = %8X, c = %8X",c1,c);
				
				//show
				//printf("\nqr_modmul:3)pVres= ");
				//pVres->value.n = pQr->M.value.n + 2;
				//print_vlong(pVres);
				//show
			}
			// 5 ---
			//vl_normal(pVres);
			j = Nm + 2;
			while (j && paU[j-1]==0) j-=1; // normalise
			pVres->value.n = j;
			// 6 ---
			while(vl_cf(pVres,&pQr->N) >= 0){//while(pV >= N)
				vl_usub(pVres,&pQr->N);// pV = pV - N
				//printf("\nqr_modmul:6)pV = pV - N; - %d",(int)(paX - pVx->value.pa));
			}
		}
		//
		vl_copy(pV,pVres);
	}
}
//--- end Quisquater reduction ---