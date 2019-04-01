#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WIN32
#include <windows.h>
#endif

#ifdef _MSC_VER
#   pragma warning( disable : 4996 )
#endif


#ifdef __GNUC__
  //
  #include <fcntl.h>
  #include <sys/mman.h>
  #include <sys/stat.h>
  #include <sys/time.h>
  #include <sys/types.h>
  #include <time.h>
  #include <unistd.h>
  #include <sys/unistd.h>
  //
  #define __stdcall
  //
#endif


#include "Vlong.h"
#include "rng.h"

void print_hex_str(const unsigned char *strHex, unsigned int uLenHex, char *strRes);
void print_hex(const unsigned char *strHex, unsigned int uLenHex);
void print_vlong(PTVlong pvl);
void print_vlong_LE(PTVlong pvl);

//
int GetTickCountMy();//time in ms
//
int main(int argc, char *argv[])
{
	
  #ifdef _MSC_VER
  SetProcessAffinityMask(GetCurrentProcess(),0x01);
  #endif

  int iTstart;
  int i;
  TVlong p,a,b,c;
  vl_init(&p);
  vl_init(&a);
  vl_init(&b);
  vl_init(&c);

  unsigned char str_p[] = {
    0x80,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
    0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x04,0x31};
  unsigned char str_a[] = {
	0x7A,0x92,0x9A,0xDE,0x78,0x9B,0xB9,0xBE,0x10,0xED,0x35,0x9D,0xD3,0x9A,0x72,0xC1,
	0x1B,0x60,0x96,0x1F,0x49,0x39,0x7E,0xEE,0x1D,0x19,0xCE,0x98,0x91,0xEC,0x3B,0x28,
	0x7A,0x92,0x9A,0xDE,0x78,0x9B,0xB9,0xBE,0x10,0xED,0x35,0x9D,0xD3,0x9A,0x72,0xC1,
	0x1B,0x60,0x96,0x1F,0x49,0x39,0x7E,0xEE,0x1D,0x19,0xCE,0x98,0x91,0xEC,0x3B,0x28
  };
  unsigned char str_b[]= {
	0x08,0xE2,0xA8,0xA0,0xE6,0x51,0x47,0xD4,0xBD,0x63,0x16,0x03,0x0E,0x16,0xD1,0x9C,
	0x85,0xC9,0x7F,0x0A,0x9C,0xA2,0x67,0x12,0x2B,0x96,0xAB,0xBC,0xEA,0x7E,0x8F,0xC8,
	0x08,0xE2,0xA8,0xA0,0xE6,0x51,0x47,0xD4,0xBD,0x63,0x16,0x03,0x0E,0x16,0xD1,0x9C,
	0x85,0xC9,0x7F,0x0A,0x9C,0xA2,0x67,0x12,0x2B,0x96,0xAB,0xBC,0xEA,0x7E,0x8F,0xC8
  };
  //
  vl_set_uchar_BE(&p,str_p,sizeof(str_p));
  vl_set_uchar_BE(&a,str_a,sizeof(str_a));
  vl_set_uchar_BE(&b,str_b,sizeof(str_b));


  printf("\na        = ");
  print_vlong(&a);
  printf("\nb        = ");
  print_vlong(&b);
  
  //1
  vl_mule(&c,&a,&b);
  printf("\na*b(raw) = ");
  print_vlong(&c);

  //2
  vl_mule_k3(&c,&a,&b);
  printf("\na*b( k3) = ");
  print_vlong(&c);

  //3
  iTstart = GetTickCountMy();
  for(i = 0; i < 1000000; ++i){
   vl_mule(&c,&a,&b);
  }
  printf("\n1000000*vl_mule(&c,&a,&b) = %d ms", GetTickCountMy() - iTstart);
  printf("\nvl_mule(&c,&a,&b) = ");
  print_vlong(&c);

  //4
  iTstart = GetTickCountMy();
  for(i = 0; i < 1000000; ++i){
   vl_mule_k3(&c,&a,&b);
  }
  printf("\n1000000*vl_mule_k3(&c,&a,&b) = %d ms", GetTickCountMy() - iTstart);
  printf("\nvl_mule_k3(&c,&a,&b) = ");
  print_vlong(&c);

  //5
  iTstart = GetTickCountMy();
  for(i = 0; i < 1000000; ++i){
   vl_add(&a,&b);
  }
  printf("\n1000000*vl_add(&a,&b); = %d ms", GetTickCountMy() - iTstart);
  printf("\nvl_add(&a,&b); = ");
  print_vlong(&a);


  //
  vl_delete(&p);
  vl_delete(&a);
  vl_delete(&b);
  vl_delete(&c);

  return 0;
}

// ------------------------------------------------------------------------


// ------------------------------------------------------------------------

void print_hex(const unsigned char *strHex, unsigned int uLenHex)
{
  char str[1024];
  str[0] = 0;
  print_hex_str(strHex,uLenHex,str);
  printf("%s",str);
};

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

void print_vlong(PTVlong pvl){
  unsigned char strx[1024];
  char str[1024];
  unsigned int i;
  str[0] = 0x30;
  str[1] = 0;
  i = sizeof(strx);
  if(vl_get_uchar_BE(strx, &i, pvl) == 1){
    print_hex_str(strx,i,str);
    printf("%s",(char*)str);
  }

}

void print_vlong_LE(PTVlong pvl){
  unsigned char strx[1024];
  char str[1024];
  unsigned int i;
  str[0] = 0x30;
  str[1] = 0;
  i = sizeof(strx);
  if(vl_get_uchar_LE(strx, &i, pvl) == 1){
    print_hex_str(strx,i,str);
    printf("%s",(char*)str);
  }

}


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

RngError rngInit(RngMode mode, uint8_t seed[40]) {
  rnd_seed ^= (int)GetTickCountMy();
  return NO_RNG_ERROR;
}

RngError rngGet(uint8_t rng[32]) {
  int i;
  RngError ret = NO_RNG_ERROR;
  if(rng)
    for(i = 0; i<32; i++){
      rng[i] = (uint8_t)my_rand(&rnd_seed);
    }
  else
	  ret = RNG_FAILED;
  //
  return ret;
}


#ifdef __GNUC__
int GetTickCountMy(){//time in ms
  struct timespec tp;
  
  clock_gettime(1, &tp);
  
  return (tp.tv_sec*1000) + (tp.tv_nsec/1000000);
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
