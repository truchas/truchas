#include <stdlib.h>
#include <stdio.h>


#include <danu_error.h>

int main(int argc, char **argv)
{
    char dum_mess[] = "This is a test message";
    int dummy = 1;
    char dum_mess2[] = "Another dummy message";

    DANU_ERROR_MESS(dum_mess);
    DANU_WARN_MESS(dum_mess);

    danu_error_printf("This is a test message %d",dummy);
    danu_warn_printf("This is a test message %d",dummy);
    danu_debug_printf("This is a test message %d %s",dummy,dum_mess2);

    return 0;
}

