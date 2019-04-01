#include <stdlib.h>
#include "Vlong.h"

//extern int iMemVlongAlloc;

//--- TVlongValue implementation ---

static
void vlv_clear(PTVlongValue pV){ // set n to zero
	TVLongUnit *pa,*pa_end;
	pV->n = 0;
	if(pV->pa){
		pa = pV->pa;
		for(pa_end = pa + pV->z; pa < pa_end; ++pa)
			*pa = 0;
	}
}

static
void vlv_normal(PTVlongValue pV){ //
	TVLongUnit *pa,*pa_end;
	if(pa = pV->pa){
		pa_end = pa + pV->z;
		pa += pV->n;
		for( ; pa < pa_end; ++pa)
			*pa = 0;
		//
		pa = pV->pa;
		while (pV->n && pa[pV->n - 1]==0) pV->n-=1; // normalise
	}
}

static inline
TVLongUnit vlv_get(const PTVlongValue pV, unsigned int i ){     // vlv_get i-th unit
	if((i < pV->n) && pV->pa){
			return pV->pa[i];
	}
	return 0;
}

static
void vlv_reserve(PTVlongValue pV, unsigned int x)
{
	unsigned int i,N;
	if (x > pV->z)
	{
		TVLongUnit *na = (TVLongUnit*)malloc(sizeof(TVLongUnit)*x);//iMemVlongAlloc+=sizeof(TVLongUnit)*x;
		if(pV->n > x) pV->n = x;
		N = pV->n;
		if(pV->pa){
			TVLongUnit *pa = pV->pa;
			for(i = 0; i < N; ++i, ++pa){
				na[i] = *pa;
				*pa = 0;
			}
			free(pV->pa);//iMemVlongAlloc-=sizeof(TVLongUnit)*pV->z;
		}
		else{
			for(i = 0; i < N; ++i){
				na[i] = 0;
			}
		}
		for(i = N; i < x; ++i) na[i] = 0;
		pV->pa = na;
		pV->z = x;
	}
}

static
void vlv_set(PTVlongValue pV, unsigned int i, TVLongUnit x ){   // set i-th unit to x
	TVLongUnit *pa;
	if ( i < pV->n )
	{
		pa = pV->pa;
		if(pa){
			pa[i] = x;
			if (x==0) 
				while (pV->n && pa[pV->n - 1]==0) pV->n-=1; // normalise
		}
	}
	else if ( x )
	{
		if(i <= pV->z) vlv_reserve(pV, i+1);
		//pa = pV->pa;
		//for (j=pV->n;j<pV->z;j+=1) pa[j] = 0;
		pV->pa[i] = x;
		pV->n = i+1;
	}
}

static inline
void vlv_init(PTVlongValue pV, TVLongUnit x ){
	pV->pa = 0;
	pV->n = 0;
	pV->z = 0;
	vlv_set(pV,0,x);
}

static inline
void vlv_inite(PTVlongValue pV, unsigned int N){
	pV->pa = 0;
	pV->n = 0;
	pV->z = 0;
	vlv_reserve(pV, N);
}

static
void vlv_delete(PTVlongValue pV){
	if(pV->pa){
		vlv_clear(pV);
		free(pV->pa);//iMemVlongAlloc-=sizeof(TVLongUnit)*pV->z;
		pV->pa = 0;
	}
	pV->z = 0;
}

static
void vlv_copy(PTVlongValue pVdst, const PTVlongValue pVsrc)
{
	unsigned int i,N;
	TVLongUnit *src_pa,*dst_pa;
	pVdst->n = 0;
	if(pVdst->z < pVsrc->n)
		vlv_reserve(pVdst, pVsrc->n);
	//
	src_pa = pVsrc->pa;
	dst_pa = pVdst->pa;
	N = pVsrc->n;
	if(src_pa){
		for(i = 0; i < N; i++, dst_pa++, src_pa++) *dst_pa = *src_pa;
		//for(     ; i < pVdst->z; i++) pVdst->pa[i] = 0;
	}
	pVdst->n = N;
}

//
static inline
int vlv_is_zero(PTVlongValue pV){
	return (pV->n == 0);
}

static inline
int vlv_test(PTVlongValue pV, unsigned int i){  //test i-th bit
	return ( vlv_get(pV,i/BPU) & ((TVLongUnit)1<<(i%BPU)) ) != 0;
}

static
unsigned int vlv_bits(PTVlongValue pV){
	unsigned int x = 0;
	TVLongUnit a;
	if(pV->n && pV->pa){
		x = pV->n;
		do{
			//a = vlv_get(pV,x-1);
			a = pV->pa[x - 1];
			if(a){
				x *= BPU;
				while(!(a & ((TVLongUnit)1<<(BPU-1)))){
					x-=1;
					a<<=1;
				};
				break;
			}else{
				x-=1;
			}
		}while(x);
		//while (x && vlv_test(pV,x-1)==0) x -= 1;//may be faster
	}
	return x;
}

static
int vlv_cf(PTVlongValue pVl, PTVlongValue pVr){
	unsigned int i;
	TVLongUnit *lpa,*rpa, l = 0, r = 0;
	if ( pVl->n != pVr->n ){
		if ( pVl->n > pVr->n ) return +1;
		if ( pVl->n < pVr->n ) return -1;
	}
	i = pVl->n;//pVl->n == pVr->n
	lpa = pVl->pa;
	rpa = pVr->pa;
	while (i)
	{
		i -= 1;
		l = lpa[i];//vlv_get(pVl,i);
		r = rpa[i];//vlv_get(pVr,i);
		if ( l != r )
			break;
	}
	//
	if ( l > r ) return +1;
	if ( l < r ) return -1;
	//
	return 0;
}

static
void vlv_shl(PTVlongValue pV)
{
	unsigned int i, N;
	TVLongUnit *pa, u, carry = 0;
	N = pV->n;
	pa = pV->pa;
	for (i=0;i < N;i+=1)
	{
		u = pa[i];
		pa[i] = (u<<1) | carry;
		carry = u>>(BPU-1);
	}
	//
	if(carry)
		vlv_set(pV,i,carry);
}

static
void vlv_shlx(PTVlongValue pV, unsigned int x){
	unsigned int i, N, BPUmx, delta = 0;//carry = 0,
	TVLongUnit *pa, u1, u2;
	//
	if(pV->n == 0 || x == 0)
		return;
	//
	if(x >= BPU)
	{
		delta = x/BPU;
		x = x%BPU;
	}
	BPUmx = BPU - x;
	//
	N = pV->n + delta;
	if(x) N+=1;
	if(N > pV->z) vlv_reserve(pV, N);
	vlv_normal(pV);
	//set locale vars
	pa = pV->pa;
	N = pV->n;
	//
	i = N;
	if(x){
		do
		{
			i--;
			//
			u2 = u1 = pa[i];
			u1 >>= BPUmx;
			u2 <<= x;
			pa[i] = 0;
			pa[i + delta + 1] |= u1;
			pa[i + delta] = u2;
			//
		}while (i > 0);
		N += delta + 1;
	}else{
		do
		{
			i--;
			//
			pa[i + delta] = pa[i];
			pa[i] = 0;
			//
		}while (i > 0);
		N += delta;
	}
	//
	//N += delta;if(x) N+=1;
	while (N && pa[N - 1]==0) N-=1; // normalise
	pV->n = N;
}

static
void vlv_shr(PTVlongValue pV)
{
	unsigned int i;
	TVLongUnit *pa, u, carry = 0;
	i = pV->n;
	pa = pV->pa;
	while (i)
	{
		i -= 1;
		u = pa[i];//vlv_get(pV,i);
		pa[i] = (u>>1) | carry;//vlv_set(pV,i,(u>>1)+carry);
		carry = u<<(BPU-1);
	}
	while (pV->n && pa[pV->n - 1]==0) pV->n-=1; // normalise
}

static inline
void vlv_shrx(PTVlongValue pV, unsigned int x)
{
	unsigned int delta = x/BPU, i, j, N;
	TVLongUnit u, u_delta1,*pa;
	x %= BPU;
	pa = pV->pa;
	N = pV->n;
	//
	if(x){
		if(delta < N) u_delta1 = pa[delta];
		else u_delta1 = 0;
		//
		for (i=0,j=delta+1;i<N;++i,++j){
			u = u_delta1;
			if(j < N) u_delta1 = pa[j];
			else u_delta1 = 0;
			u >>= x;
			u |= (u_delta1 << (BPU-x));
			pa[i] = u;
		}
	}else{
		for (i=0,j=delta;j<N;++i,++j){
			pa[i] = pa[j];
		}
		for (;i<N;++i){
			pa[i] = 0;
		}
	}
	//
	if(delta <= N) N-=delta;
	while (N && !pa[N - 1]) N-=1; // normalise
	pV->n = N;
}

static
void vlv_add(PTVlongValue pV, const PTVlongValue pVx)//pV = pV + pVx
{
	unsigned int i;
	unsigned int max;
	TVLongUnit *pa,u,ux,carry = 0;
	max = pV->n;
	if (max < pVx->n) max = pVx->n;
	max++;
	if(pV->z < max)
		vlv_reserve(pV,max);
	pa = pV->pa;
	for (i=0;i < max; i+=1)
	{
		u = vlv_get(pV,i);
		u = u + carry; carry = ( u < carry );
		ux = vlv_get(pVx,i);
		u = u + ux; carry += ( u < ux );
		pa[i] = u;
	}
	while (max && !pa[max - 1]) max-=1; // normalise
	pV->n = max;
}

static
void vlv_subtract(PTVlongValue pV, const PTVlongValue pVx)//pV = pV - pVx
{
	unsigned int N, Nx, i;
	TVLongUnit *pa, *pax, u, ux, nu, carry = 0;
	N = pV->n;
	pa = pV->pa;
	Nx = pVx->n;
	pax = pVx->pa;
	for (i=0;i<N;++i)
	{
		//ux = vlv_get(pVx,i);
		if(i < Nx) ux = pax[i] + carry;
		else ux = carry;
		//
		//ux += carry;
		if(ux)// ( ux >= carry )
		{
			u = pa[i];//vlv_get(pV,i);
			nu = u - ux;
			carry = nu > u;
			pa[i] = nu;
		}
	}
	while (N && !pa[N - 1]) N-=1; // normalise
	pV->n = N;
}


typedef union unionu_tvlong{
	unsigned long long udub;
	TVLongUnit pu[2];
} unionu_tvlong;

static inline
void vlv_fast_mul(PTVlongValue pVres, const PTVlongValue pVx, const PTVlongValue pVy, unsigned int keep)
{
	// pVres = (pVx*pVy) % (2^keep)
	unsigned int i, limit = (keep + BPU - 1) / BPU; // size of result in ints
	unsigned int min, j, min_;
	TVLongUnit *res_pa, *lres_pa, *x_pa, *y_pa, m, c;
	unionu_tvlong ww;

	if (pVres && pVx && pVy) {
		//
		if (pVres->z < limit) {
			pVres->n = 0;
			vlv_reserve(pVres, limit);
		}
		vlv_clear(pVres);
		//
		res_pa = pVres->pa;
		x_pa = pVx->pa;
		//y_pa = pVy->pa;
		//
		min = pVx->n; if (min > limit) min = limit;

		for (i = 0; i < min; ++i, ++x_pa)
		{
			if(m = *x_pa)//m!=0
			{
				c = 0;  //carry
				min_ = i + pVy->n; if (min_ > limit) min_ = limit;
				y_pa = pVy->pa;
				lres_pa = res_pa + i;
				//
				for (j = i; j < min_; j += 1, ++y_pa, ++lres_pa)
				{
					// This is the critical loop
					// Machine dependent code could help here
					// c:a[j] = a[j] + c + m*y.a[j-i];
					ww.udub = (unsigned long long)m;
					ww.udub *= (unsigned long long)*y_pa;//p;//pVy->pa[j-i];
					//
					ww.udub += (unsigned long long)c;
					ww.udub += (unsigned long long)*lres_pa;
					//
					*lres_pa = ww.pu[0];//(unsigned int)ww;
					c = ww.pu[1];//(unsigned int)(ww >> BPU);
				}

				ww.pu[0] = c;
				ww.pu[1] = 0;
				
				while (ww.pu[0] && (j < limit))
				{
					ww.udub += (unsigned long long)*lres_pa;
					*lres_pa = ww.pu[0];
					ww.pu[0] = ww.pu[1];
					ww.pu[1] = 0;
					++j;
					++lres_pa;
				}
			}
		}

		// eliminate unwanted bits
		keep %= BPU; if (keep) res_pa[limit - 1] &= (1 << keep) - 1;

		// calculate n
		while (limit && !res_pa[limit - 1]) limit -= 1;
		pVres->n = limit;
	}
}

static inline
void vlv_fast_add_mul(PTVlongValue pVres, const PTVlongValue pVx, const PTVlongValue pVy, unsigned int keep)
{
	// pVres = (pVx*pVy) % (2^keep)
	unsigned int i, limit = (keep + BPU - 1) / BPU; // size of result in ints
	unsigned int min, j, min_;
	TVLongUnit *res_pa, *lres_pa, *x_pa, *y_pa, m, c;
	unionu_tvlong ww;

	if (pVres && pVx && pVy) {
		//
		if (pVres->z < limit) {
			vlv_reserve(pVres, limit);
		}
		//
		res_pa = pVres->pa;
		x_pa = pVx->pa;
		//y_pa = pVy->pa;
		//
		min = pVx->n; if (min > limit) min = limit;

		for (i = 0; i < min; ++i, ++x_pa)
		{
			if(m = *x_pa)//m!=0
			{
				c = 0;  //carry
				min_ = i + pVy->n; if (min_ > limit) min_ = limit;
				y_pa = pVy->pa;
				lres_pa = res_pa + i;
				//
				for (j = i; j < min_; j += 1, ++y_pa, ++lres_pa)
				{
					// This is the critical loop
					// Machine dependent code could help here
					// c:a[j] = a[j] + c + m*y.a[j-i];
					ww.udub = (unsigned long long)m;
					ww.udub *= (unsigned long long)*y_pa;//p;//pVy->pa[j-i];
					//
					ww.udub += (unsigned long long)c;
					ww.udub += (unsigned long long)*lres_pa;
					//
					*lres_pa = ww.pu[0];//(unsigned int)ww;
					c = ww.pu[1];//(unsigned int)(ww >> BPU);
				}

				while (c && j < limit)
				{
					*lres_pa += c;
					c = *lres_pa < c;
					++j;
					++lres_pa;
				}
			}
		}

		// eliminate unwanted bits
		keep %= BPU; if (keep) res_pa[limit - 1] &= (1 << keep) - 1;

		// calculate n
		while (limit && !res_pa[limit - 1]) limit -= 1;
		pVres->n = limit;
	}
}

static
void vlv_mul(PTVlongValue pVres, const PTVlongValue pVx, const PTVlongValue pVy)
{
	if(pVres && pVx && pVy)
		//vlv_fast_mul(pVres, pVx, pVy, vlv_bits(pVx)+vlv_bits(pVy) );
		vlv_fast_mul(pVres, pVx, pVy, (pVx->n+pVy->n)*BPU );
}

static inline
TVLongUnit vlv_mod_word(const PTVlongValue pVres, TVLongUnit w){
	int i;
	unionu_tvlong ww;
	const TVLongUnit *res_pa;
	ww.udub = 0;
	if(pVres && (w > 1)){
		res_pa = pVres->pa;
		for(i = pVres->n - 1; i>=0; i-=1){
			ww.pu[1] = ww.pu[0];//<<
			ww.pu[0] = res_pa[i];
			ww.udub %= (unsigned long long)w;
		}
	}
	return ww.pu[0];
}

static
void vlv_divide(PTVlongValue pVres, PTVlongValue pVrem, const PTVlongValue pVx, const PTVlongValue pVy)//x = res*y + rem
{
	int i;
	TVlongValue m,s;

	vlv_clear(pVres);
	vlv_copy(pVrem,pVx);

	vlv_inite(&m,pVy->n);
	vlv_copy(&m,pVy);

	vlv_init(&s,1);
	//
	i = vlv_bits(pVrem) - vlv_bits(&m);
	if(i > 0){
		vlv_shlx(&m,i);
		vlv_shlx(&s,i);
	}
	while ( vlv_cf(pVrem,&m) > 0 ){
		vlv_shl(&m);
		vlv_shl(&s);
	}
	//
	while ( vlv_cf(pVrem,pVy) >= 0 )
	{
		while ( vlv_cf(pVrem,&m) < 0 )
		{
			vlv_shr(&m);
			vlv_shr(&s);
		}
		vlv_subtract(pVrem,&m);
		vlv_add(pVres,&s);
	}
	//
	vlv_delete(&m);
	vlv_delete(&s);
}

//--------------------------------------
static
void vlv_set_uchar_BE(PTVlongValue pV, const unsigned char *s, unsigned int len)
{
	unsigned int i;

	i = (len + sizeof(TVLongUnit) - 1)/sizeof(TVLongUnit);//size in uint
	vlv_clear(pV);
	vlv_reserve(pV,i+1);

	i = 0;
	while (i<len)
	{
		vlv_shlx(pV,8);
		vlv_set(pV,0,vlv_get(pV,0)|(TVLongUnit)s[i]);
		i++;
	}
}

static
int vlv_get_uchar_BE(unsigned char *s, unsigned int *plen, const PTVlongValue pV){
	unsigned int i = (vlv_bits(pV)+7)/8;//trunc from octet
	TVlongValue x;
	//
	if(!s)
	{
		*plen = i;
		return 1;
	}
	//
	if(*plen < i)
	{
		*plen = i;
		return 0;
	}
	//
	if(i<1){
		*plen = 1;
		s[0] = 0;
		return 1;
	}
	//
	vlv_init(&x,0);
	*plen = i;
	//
	vlv_copy(&x,pV);
	//
	while (i)
	{
		i--;
		s[i] = (unsigned char)(vlv_get(&x,0) & 0xFF);
		vlv_shrx(&x,8);
	}
	//
	vlv_delete(&x);
	//
	return 1;
}

static
int vlv_get_uchar_LE(unsigned char *s, unsigned int *plen, const PTVlongValue pV){
	unsigned int i = (vlv_bits(pV)+7)/8;//trunc from octet
	TVlongValue x;
	//
	if(!s)
	{
		*plen = i;
		return 1;
	}
	//
	if(*plen < i)
	{
		*plen = i;
		return 0;
	}
	//
	if(i<1){
		*plen = 1;
		s[0] = 0;
		return 1;
	}
	//
	vlv_init(&x,0);
	*plen = i;
	//
	vlv_copy(&x,pV);
	//
	for(i = 0; i< *plen; i++)
	{
		s[i] = (unsigned char)(vlv_get(&x,0) & 0xFF);
		vlv_shrx(&x,8);
	}
	//
	vlv_delete(&x);
	//
	return 1;
}

//--- end TVlongValue implementation ---

//--- Implementation of TVlong ---

void vl_init(PTVlong pV){
	if(pV){
		vlv_init(&pV->value,0);
		pV->negative = 0;
	}
}

void vl_inite(PTVlong pV, unsigned int N){
	if(pV){
		vlv_inite(&pV->value,N);
		pV->negative = 0;
	}
}

void vl_delete(PTVlong pV){
	if(pV){
		vlv_delete(&pV->value);
		pV->negative = 0;
	}
}

void vl_clear(PTVlong pV){
	if(pV){
		vlv_clear(&pV->value);
		pV->negative = 0;
	}
}

void vl_normal(PTVlong pV){
	if(pV)
		vlv_normal(&pV->value);
}

void vl_reserve(PTVlong pV, unsigned int x){
	if(pV){
		vlv_reserve(&pV->value,x);
	}
}

void vl_set(PTVlong pV, unsigned int i, TVLongUnit x){
	if(pV){
		vlv_set(&pV->value, i, x);
	}
}

TVLongUnit vl_get(PTVlong pV, unsigned int i){
	if(pV){
		return vlv_get(&pV->value, i);
	}
	return 0;
}

unsigned int vl_bits(const PTVlong pV){
	if(pV){
		return vlv_bits(&pV->value);
	}
	return 0;
}

int vl_test(const PTVlong pV, unsigned int i){ //test i-th bit
	if(pV){
		return vlv_test(&pV->value,i);
	}
	return 0;
}

int vl_is_negative(const PTVlong pV){
	if(pV){
		return pV->negative;
	}
	return 0;
}

void vl_shr(PTVlong pV){
	if(pV){
		vlv_shr(&pV->value);
	}
}

void vl_shrx(PTVlong pV, unsigned int x){
	if(pV){
		vlv_shrx(&pV->value,x);
	}
}

void vl_shl(PTVlong pV){
	if(pV){
		vlv_shl(&pV->value);
	}
}

void vl_shlx(PTVlong pV, unsigned int x){
	if(pV){
		vlv_shlx(&pV->value,x);
	}
}

int vl_cf(const PTVlong pVl, const PTVlong pVr)
{
	int neg;
	//
	if(pVl && pVr){
		neg = pVl->negative && pVl->value.n;
		
		if ( neg == (pVr->negative && pVr->value.n) ){
			if(neg)
				return -vlv_cf(&pVl->value,&pVr->value);//(* -1 )(if neg == 1)!!!
			else
				return vlv_cf(&pVl->value,&pVr->value);
		}
		else if ( neg ) return -1;
		else return +1;
	}
	return 0;//err
}

void vl_copy(PTVlong pVdst, const PTVlong pVsrc){
	if((pVdst != pVsrc) && pVdst && pVsrc){
		vlv_copy(&pVdst->value,&pVsrc->value);
		pVdst->negative = pVsrc->negative;
	}
}

void vl_add(PTVlong pV, const PTVlong pVx){//pV = pV + pVx
	TVlong Vtmp;
	if(pV && pVx){
		if ( pV->negative == pVx->negative )
		{
			vlv_add(&pV->value,&pVx->value);
		}
		else if ( vlv_cf(&pV->value,&pVx->value) >= 0 )
		{
			vlv_subtract(&pV->value,&pVx->value);
		}
		else
		{
			vl_init(&Vtmp);
			vl_copy(&Vtmp,pVx);
			vl_add(&Vtmp,pV);
			vl_copy(pV,&Vtmp);
			vl_delete(&Vtmp);
		}
	}
}

void vl_word_uadd(PTVlong pV, TVLongUnit a){
	if(pV && a){
		unsigned int i;
		TVLongUnit m = 0;
		for(i = 0; i < pV->value.z; ++i){
			m += vl_get(pV, i) + a;
			vl_set(pV,i,m);
			if(a < m){
				break;
			}
			m = 1;//carry
		}
	}
}

void vl_inc(PTVlong pV){//pV+=1
	TVlong t;
	if(pV){
		vl_init(&t);
		vl_set(&t,0,1);
		//vl_add(pV,&t);
		vlv_add(&pV->value,&t.value);
		vl_delete(&t);
	}
}

void vl_dec(PTVlong pV){//pV-=1
	TVlong t;
	if(pV){
		vl_init(&t);
		vl_set(&t,0,1);
		vlv_subtract(&pV->value,&t.value);
		vl_delete(&t);
	}
}

void vl_usub(PTVlong pV, const PTVlong pVx){//unsigned pV = pV - pVx (pV must >= pVx)
	if(pV && pVx){
		vlv_subtract(&pV->value,&pVx->value);
	}
}

void vl_sub(PTVlong pV, const PTVlong pVx){//pV = pV - pVx
	TVlong Vtmp;
	if(pV && pVx){
		if ( pV->negative != pVx->negative )
		{
			vlv_add(&pV->value,&pVx->value);
		}
		else if ( vlv_cf(&pV->value,&pVx->value) >= 0 )
		{
			vlv_subtract(&pV->value,&pVx->value);
		}
		else
		{
			vl_init(&Vtmp);
			vl_copy(&Vtmp,pVx);
			vl_sub(&Vtmp,pV);
			vl_copy(pV,&Vtmp);
			vl_delete(&Vtmp);
			pV->negative = 1 - pV->negative;
		}
	}
}

void vl_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy){
	if(pVres && pVx && pVy){
		vlv_mul(&pVres->value,&pVx->value,&pVy->value);
		pVres->negative = pVx->negative ^ pVy->negative;
	}
}

void vl_mul(PTVlong pVx, const PTVlong pVy){//Vx*=Vy
	TVlong res;
	if(pVx && pVy){
		vl_init(&res);
		vlv_mul(&res.value,&pVx->value,&pVy->value);
		res.negative = pVx->negative ^ pVy->negative;
		vl_copy(pVx,&res);
		vl_delete(&res);
	}
}

void vl_fast_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy, unsigned int keep ){
	if(pVres && pVx && pVy){
		vlv_fast_mul(&pVres->value,&pVx->value,&pVy->value,keep);
		pVres->negative = pVx->negative ^ pVy->negative;
	}
}

void vl_fast_add_mule(PTVlong pVres, const PTVlong pVx, const PTVlong pVy, unsigned int keep ){
	if(pVres && pVx && pVy){
		vlv_fast_add_mul(&pVres->value,&pVx->value,&pVy->value,keep);
		pVres->negative = pVx->negative ^ pVy->negative;//by modulus
	}
}

TVLongUnit vl_mod_word(const PTVlong pVres, TVLongUnit w){//return pVres%w
	//
	return vlv_mod_word(&pVres->value,w);
}

void vl_dive(PTVlong pVres, PTVlong pVrem, const PTVlong pVx, const PTVlong pVy){//x = res*y + rem; res = x/y rem = x%y
	TVlong Vres,Vrem;
	if(pVx && pVy){
		vl_init(&Vres);
		vl_init(&Vrem);
		//
		vlv_divide(&Vres.value,&Vrem.value,&pVx->value,&pVy->value);
		//
		if(pVres){
			vlv_copy(&pVres->value,&Vres.value);
			pVres->negative = pVx->negative ^ pVy->negative;
		}
		if(pVrem){
			vlv_copy(&pVrem->value,&Vrem.value);
			pVrem->negative = pVx->negative;//?
		}
		//
		vl_delete(&Vres);
		vl_delete(&Vrem);
	}
}

int vl_fast_compare(const PTVlong pVx, const PTVlong pVy){
	unsigned int i;
	int iRes = 0;
	//
	if(pVx && pVy){
		iRes = 0;
		if(pVx->negative != pVy->negative)
			return iRes;
		//
		if(pVx->value.n != pVy->value.n)
			return iRes;
		//
		iRes = 1;
		//
		for(i = 0; i < pVx->value.n; i++)
			if(vlv_get(&pVx->value,i) != vlv_get(&pVy->value,i))
			{
				iRes = 0;
				break;
			}
	}
	//
	return iRes;
}

//--- in-out ---

void vl_set_uchar_BE(PTVlong pV, const unsigned char *s, unsigned int len){
	if(pV && s){
		vl_clear(pV);
		vlv_set_uchar_BE(&pV->value,s,len);
	}
}

int vl_get_uchar_BE(unsigned char *s, unsigned int *plen, const PTVlong pV){
	if(pV){
		return vlv_get_uchar_BE(s,plen,&pV->value);
	}
	return -1;
}

int vl_get_uchar_LE(unsigned char *s, unsigned int *plen, const PTVlong pV){
	if(pV){
		return vlv_get_uchar_LE(s,plen,&pV->value);
	}
	return -1;
}

int vl_set_vlong(PTVlong pV, const TVLongUnit *puiVL, unsigned int uiLength){
	if(pV && puiVL && uiLength>0){
		vlv_clear(&pV->value);
		vlv_reserve(&pV->value,uiLength);
		pV->negative = 0;
		while (uiLength) { uiLength -= 1; vlv_set(&pV->value, uiLength, puiVL[uiLength]); }
		return 1;
	}

	return 0;
}

int vl_get_vlong(TVLongUnit *puiVL, unsigned int *puiLength, const PTVlong pV){//if puiVL==0, puiLength - sizeof puiVL
	unsigned int i;
	if(pV){
		i = pV->value.n;
		if(puiVL == 0)
		{
			*puiLength = i;
			return 1;
		}
		//
		if(*puiLength < i) return 0;
		*puiLength = i;
		//
		while (i) { i -= 1; puiVL[i] = vlv_get(&pV->value,i); }
		return 1;
	}
	return 0;
}

//--- end Implementation of TVlong ---

