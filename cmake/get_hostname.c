#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv)
{
  char hostname[128];
  int exit;
  exit=gethostname(hostname,128);
  printf("%s",hostname);
  return exit;
}
