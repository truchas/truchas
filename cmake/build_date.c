#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main ( int argc, char **argv)
{
  char date[128];
  time_t now = time(0);
  strftime(date,sizeof(date),"%e %b %Y %H:%M:%S",gmtime(&now));
  printf("%s",date);
  return 0;
}
