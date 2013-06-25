#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//int foo(const* const char* const a, const int num);

#define CONST_CHAR_PTR_ARRAY (const char* const* )
typedef char* const* char_const_ptr_t;
//typedef const_char_array_t (const char)* const*;

//int foo(const char* const* a, const int num)
int foo(const char_const_ptr_t a, const int num)
{
   int cnt=0;
   int i,len;

   for(i=0;i<num;i++) {
     len=strlen(a[i]);
     printf("(%d) String:%s has %d characters\n",i,a[i],len);
     cnt+=len;
   }

   return cnt;
}

int bar(const int x[], const int num)
{
  int i;

  for(i=0;i<num-1;i++)
    printf("%d,",x[i]);

  printf("%d\n",x[num-1]);

  return 0;
}

int main(int argc, char ** argv)
{
  int char_cnt;
  int b[8] = { 0, 1, 2, 3, 4, 5, 6, 7};

  char_cnt = foo((const char * const*)argv,argc);
  bar(b,8);

  printf("Number of chars in argv %d\n", char_cnt);

  return 0;
}


