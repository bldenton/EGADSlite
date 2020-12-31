/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Memory Handling Functions
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "egadsTypes.h"
#include "egadsStack.h"


/*@null@*/ /*@out@*/ /*@only@*/ void *
EG_alloc(size_t nbytes)
{
  return malloc(nbytes);
}


/*@null@*/ /*@only@*/ void *
EG_calloc(size_t nele, size_t size)
{
  return calloc(nele, size);
}


/*@null@*/ /*@only@*/ void *
EG_reall(/*@null@*/ /*@only@*/ /*@returned@*/ void *ptr, size_t nbytes)
{
  return realloc(ptr, nbytes);
}


void
EG_free(/*@null@*/ /*@only@*/ void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}


/*@null@*/ /*@only@*/ char *
EG_strdup(/*@null@*/ const char *str)
{
  int  i, len;
  char *dup;

  if (str == NULL) return NULL;

  len = strlen(str) + 1;
  dup = (char *) EG_alloc(len*sizeof(char));
  if (dup != NULL)
    for (i = 0; i < len; i++) dup[i] = str[i];

  return dup;
}


int
EG_stackInit(objStack *stack)
{
  stack->stkPtr = -1;
  stack->stkSiz = 256;
  stack->stack  = (egObject **) EG_alloc(256*sizeof(egObject *));
  if (stack->stack == NULL) return EGADS_MALLOC;
  
  return EGADS_SUCCESS;
}


void
EG_stackFree(objStack *stack)
{
  if (stack->stack != NULL) EG_free(stack->stack);
  stack->stack  = NULL;
  stack->stkPtr = -1;
}


int
EG_stackPush(objStack *stack, egObject *entity)
{
  egObject **tmp;
  
  if (stack->stack == NULL) return EGADS_MALLOC;
  if (entity == NULL) {
    printf(" EGADS Internal: EG_stackPush with NULL Entity!\n");
    return EGADS_NULLOBJ;
  }
  if (stack->stkSiz == stack->stkPtr+1) {
    tmp = (egObject **) EG_reall(stack->stack,
                                 (stack->stkSiz+256)*sizeof(egObject *));
    if (tmp == NULL) return EGADS_MALLOC;
    stack->stack   = tmp;
    stack->stkSiz += 256;
  }
  stack->stkPtr++;
  stack->stack[stack->stkPtr] = entity;
  
  return EGADS_SUCCESS;
}


void
EG_stackPop(objStack *stack, egObject **entity)
{
  *entity = NULL;
  if (stack->stkPtr < 0) return;
  *entity = stack->stack[stack->stkPtr];
  stack->stkPtr--;
}
