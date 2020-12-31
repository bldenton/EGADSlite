#ifndef EGADSSTACK_H
#define EGADSSTACK_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Stack function Header
 *
 *      Copyright 2011-2020, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */


typedef struct {
  int      stkPtr;
  int      stkSiz;
  egObject **stack;
} objStack;


#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__
#else
#define __ProtoExt__ extern
#endif
  
__ProtoExt__ int  EG_stackInit( objStack *stack );
__ProtoExt__ void EG_stackFree( objStack *stack );
__ProtoExt__ int  EG_stackPush( objStack *stack, egObject *entity );
__ProtoExt__ void EG_stackPop( objStack *stack, egObject **entity );

#ifdef __cplusplus
}
#endif

#endif
