/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/* Exception Handling Code for SWIG wrappers */

#include <string.h>

#define DANU_PYTHON_OK    0
#define DANU_PYTHON_FAIL -1
#define MAX_ERROR_MESSAGE 256

static char error_message[MAX_ERROR_MESSAGE];
static int error_status = 0;

void throw_exception(char *msg) {
  strncpy(error_message,msg,MAX_ERROR_MESSAGE);
  error_status = DANU_PYTHON_FAIL;
}

void clear_exception() {
  error_status = DANU_PYTHON_OK;
}

char *check_exception() {
  if (error_status == DANU_PYTHON_FAIL) {
    return error_message;
  } 
  else {
    return NULL;
  }

}


