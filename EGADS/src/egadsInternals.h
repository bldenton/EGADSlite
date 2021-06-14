#ifndef EGADSINTERNALS_H
#define EGADSINTERNALS_H
/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             Internal function Header
 *
 *      Copyright 2011-2021, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#define PI       3.1415926535897931159979635
#define KNDIFF   1.0e-5         // knot difference


#ifdef __HOST_AND_DEVICE__
#undef __HOST_AND_DEVICE__
#endif

#ifdef __CUDACC__
#define __HOST_AND_DEVICE__ __host__ __device__
#else
#define __HOST_AND_DEVICE__
#endif

#ifdef __ProtoExt__
#undef __ProtoExt__
#endif
#ifdef __cplusplus
extern "C" {
#define __ProtoExt__ __HOST_AND_DEVICE__
#else
#define __ProtoExt__ __HOST_AND_DEVICE__ extern
#endif


__ProtoExt__ /*@null@*/ /*@out@*/ /*@only@*/ 
             void *EG_alloc( size_t nbytes );
__ProtoExt__ /*@null@*/ /*@only@*/ 
             void *EG_calloc( size_t nele, size_t size );
__ProtoExt__ /*@null@*/ /*@only@*/ 
             void *EG_reall( /*@null@*/ /*@only@*/ /*@returned@*/ void *ptr,
                             size_t nbytes );
__ProtoExt__ void EG_free( /*@null@*/ /*@only@*/ void *pointer );
__ProtoExt__ /*@null@*/ /*@only@*/ 
             char *EG_strdup( /*@null@*/ const char *str );

__ProtoExt__ /*@kept@*/ /*@null@*/ egObject *
                  EG_context( const egObject *object );
__ProtoExt__ int  EG_sameThread( const egObject *object );
__ProtoExt__ int  EG_outLevel( const egObject *object );
__ProtoExt__ int  EG_makeObject( /*@null@*/ egObject *context, egObject **obj );
__ProtoExt__ int  EG_deleteObject( egObject *object );
__ProtoExt__ int  EG_dereferenceObject( egObject *object,
                                        /*@null@*/ const egObject *ref );
__ProtoExt__ int  EG_dereferenceTopObj( egObject *object,
                                        /*@null@*/ const egObject *ref );
__ProtoExt__ int  EG_referenceObject( egObject *object, 
                                      /*@null@*/ const egObject *ref );
__ProtoExt__ int  EG_referenceTopObj( egObject *object, 
                                      /*@null@*/ const egObject *ref );
__ProtoExt__ int  EG_removeCntxtRef( egObject *object );

__ProtoExt__ int  EG_attributeDel( egObject *obj, /*@null@*/ const char *name );
__ProtoExt__ int  EG_attributeDup( const egObject *src, egObject *dst );
__ProtoExt__ int  EG_attributeXDup( const egObject *src,
                                    /*@null@*/ const double *xform,
                                          egObject *dst );
__ProtoExt__ int  EG_attributePrint( const egObject *src );

#ifdef __cplusplus
}
#endif

#undef __HOST_AND_DEVICE__

#endif
