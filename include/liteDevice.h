#ifndef LITE_DEVICE_H
#define LITE_DEVICE_H

#include <stdlib.h>
#include <string.h>

#include "cudaUtil.h"

/**
 * The following macros are EGADS specific to perform Device operations
 * on EGADS data structures.
 **/

#if defined(__NVCC__) && !defined(__CUDA_ARCH__)

# define EG_NEW(ptraddr, type, nmemb) DEVICE_MALLOC((void**)(ptraddr), nmemb*sizeof(type))

# define EG_REALLOC(dstaddr, src, type, old_nmemb, nmemb) DEVICE_REALLOC(dstaddr, src, type, old_nmemb, nmemb)

# define EG_FREE(ptr) DEVICE_FREE((ptr))

# define EG_COPY(ptr_d, ptr, type, nmemb) MEMCPY_HOST_TO_DEVICE(ptr_d, ptr, nmemb*sizeof(type))

# define EG_GET(ptr, ptr_d, type, nmemb) MEMCPY_DEVICE_TO_HOST(ptr, ptr_d, nmemb*sizeof(type))

# define EG_SET_STR(string_daddr, string) STRDUP_HOST_TO_DEVICE(string_daddr, string)

# define EG_GET_STRN(string, string_daddr, N) STRNCPY_DEVICE_TO_HOST(string, string_daddr, N)

# define EG_NEW_OBJECT(objectaddr) DEVICE_MALLOC((void**)(objectaddr), sizeof(egObject))

# define EG_FREE_OBJECT(object) DEVICE_FREE((object))

# define EG_GET_OBJECT_PTR(objectaddr, object_daddr) MEMCPY_DEVICE_TO_HOST(objectaddr, object_daddr, sizeof(egObject*))

# define EG_GET_OBJECT(object, object_d) MEMCPY_DEVICE_TO_HOST(object, object_d, sizeof(egObject))

# define EG_SET_OBJECT_PTR(object_daddr, objectaddr) MEMCPY_HOST_TO_DEVICE(object_daddr, objectaddr, sizeof(egObject*))

# define EG_SET_OBJECT(object_daddr, object) MEMCPY_HOST_TO_DEVICE(*(object_daddr), object, sizeof(egObject))

# define EG_NEW_CNTXT(cntxaddr)                                                \
do {                                                                           \
printf("allocating new context\n");\
  DEVICE_MALLOC((void**)(cntxaddr), sizeof(egCntxt));                          \
  if (NULL != *(cntxaddr)) {                                                   \
printf("allocating new strings\n");\
    char **strings_d;                                                          \
    DEVICE_MALLOC(&strings_d, 2*sizeof(char*));                                \
    if (NULL == strings_d) {                                                   \
      DEVICE_FREE(*(cntxaddr));                                                \
      break;                                                                   \
    }                                                                          \
printf("copying strings\n");\
    MEMCPY_HOST_TO_DEVICE(&((*(cntxaddr))->signature), &strings_d, sizeof(char**));\
printf("have new context\n");\
  }                                                                            \
} while(0)

# define EG_FREE_CNTXT(cntx)                                                   \
do {                                                                           \
  if (NULL != cntx) {                                                          \
    char **strings_d;                                                          \
    MEMCPY_DEVICE_TO_HOST(&strings_d, &(cntx->signature), sizeof(char**));     \
    if (EGADS_SUCCESS != status) break;                                        \
    if (NULL != strings_d) DEVICE_FREE(strings_d);                             \
    if (EGADS_SUCCESS != status) break;                                        \
  }                                                                            \
  DEVICE_FREE((cntx));                                                         \
} while(0)

# define EG_GET_CNTXT(cntx, cntx_d) MEMCPY_DEVICE_TO_HOST(cntx, cntx_d, sizeof(egCntxt))

# define EG_SET_CNTXT(cntx_d, cntx)                                               \
do {                                                                           \
  char **strings_d;                                                            \
  char *string_d;                                                              \
  long zero = 0L;                                                              \
  void *nil = NULL;                                                            \
  MEMCPY_DEVICE_TO_HOST(&strings_d, &(cntx_d->signature), sizeof(char**));     \
  if (EGADS_SUCCESS != status) break;                                          \
  MEMCPY_HOST_TO_DEVICE(cntx_d, cntx, sizeof(egCntxt));                        \
  if (EGADS_SUCCESS == status) {                                               \
    MEMCPY_HOST_TO_DEVICE(&(cntx_d->signature), &strings_d, sizeof(char**));   \
    if (EGADS_SUCCESS == status && NULL != strings_d) {                        \
      MEMCPY_DEVICE_TO_HOST(&string_d, &(strings_d[1]), sizeof(char*));        \
      if (EGADS_SUCCESS == status && NULL != string_d) {                       \
        DEVICE_FREE(string_d);                                                 \
        if (EGADS_SUCCESS != status) break;                                    \
      }                                                                        \
      STRDUP_HOST_TO_DEVICE(&(string_d), cntx->signature[1]);                  \
      MEMCPY_DEVICE_TO_HOST(&string_d, &(strings_d[0]), sizeof(char*));        \
      if (EGADS_SUCCESS == status && NULL != string_d) {                       \
        DEVICE_FREE(string_d);                                                 \
        if (EGADS_SUCCESS != status) break;                                    \
      }                                                                        \
      STRDUP_HOST_TO_DEVICE(&(string_d), cntx->signature[0]);                  \
    }                                                                          \
  }                                                                            \
  if (EGADS_SUCCESS != status) break;                                          \
  MEMCPY_HOST_TO_DEVICE(&(cntx_d->threadID), &zero, sizeof(long));             \
  if (EGADS_SUCCESS != status) break;                                          \
  MEMCPY_HOST_TO_DEVICE(&(cntx_d->mutex), &nil, sizeof(void*));                \
} while(0)

# define EG_GET_ATTRS_PTR(attrs, attrs_d) MEMCPY_DEVICE_TO_HOST(&(attrs), &(attrs_d), sizeof(egAttrs*))

# define EG_SET_ATTRS_PTR(attrs_daddr, attrs) MEMCPY_HOST_TO_DEVICE(attrs_daddr, &(attrs), sizeof(egAttrs*))

# define EG_GET_ATTRS(attrs, attrs_d) MEMCPY_DEVICE_TO_HOST(attrs, attrs_d, sizeof(egAttrs))

# define EG_SET_ATTRS(attrs_daddr, attrs) MEMCPY_HOST_TO_DEVICE(attrs_daddr, attrs, sizeof(egAttrs))

# define EG_GET_ATTR_PTR(attr, attr_d) MEMCPY_DEVICE_TO_HOST(&(attr), &(attr_d), sizeof(egAttr*))

# define EG_GET_ATTR(attr, attr_daddr) MEMCPY_DEVICE_TO_HOST(attr, attr_daddr, sizeof(egAttr))

# define EG_SET_ATTR(attr_daddr, attr) MEMCPY_HOST_TO_DEVICE(attr_daddr, attr, sizeof(egAttr))

# define EG_GET_TESSEL(tess, tess_d) MEMCPY_DEVICE_TO_HOST(tess, tess_d, sizeof(egTessel))

# define EG_GET_TESS1D(tess1d, tess1d_daddr)                                   \
do {                                                                           \
  MEMCPY_DEVICE_TO_HOST(tess1d, tess1d_daddr, sizeof(egTess1D));               \
} while(0)

# define EG_GET_TESS2D(tess2daddr, tess2d_daddr)                               \
do {                                                                           \
  MEMCPY_DEVICE_TO_HOST(tess2daddr, tess2d_daddr, sizeof(egTess2D));           \
} while(0)

# define EG_GET_PATCH(patch, patch_daddr)                                      \
do {                                                                           \
  MEMCPY_DEVICE_TO_HOST(patch, patch_daddr, sizeof(egPatch));                  \
} while(0)

# define EG_GET_GEOM(geom, geom_d) MEMCPY_DEVICE_TO_HOST(geom, geom_d, sizeof(liteGeometry))

# define EG_SET_GEOM(geom_d, geom) MEMCPY_HOST_TO_DEVICE(geom_d, geom, sizeof(liteGeometry))

# define EG_GET_NODE(node, node_d) MEMCPY_DEVICE_TO_HOST(node, node_d, sizeof(liteNode))

# define EG_SET_NODE(node_d, node) MEMCPY_HOST_TO_DEVICE(node_d, node, sizeof(liteNode))

# define EG_GET_EDGE(edge, edge_d) MEMCPY_DEVICE_TO_HOST(edge, edge_d, sizeof(liteEdge))

# define EG_SET_EDGE(edge_d, edge) MEMCPY_HOST_TO_DEVICE(edge_d, edge, sizeof(liteEdge))

# define EG_GET_LOOP(loop, loop_d) MEMCPY_DEVICE_TO_HOST(loop, loop_d, sizeof(liteLoop))

# define EG_SET_LOOP(loop_d, loop) MEMCPY_HOST_TO_DEVICE(loop_d, loop, sizeof(liteLoop))

# define EG_GET_FACE(face, face_d) MEMCPY_DEVICE_TO_HOST(face, face_d, sizeof(liteFace))

# define EG_SET_FACE(face_d, face) MEMCPY_HOST_TO_DEVICE(face_d, face, sizeof(liteFace))

# define EG_GET_SHELL(shell, shell_d) MEMCPY_DEVICE_TO_HOST(shell, shell_d, sizeof(liteShell))

# define EG_SET_SHELL(shell_d, shell) MEMCPY_HOST_TO_DEVICE(shell_d, shell, sizeof(liteShell))

# define EG_GET_BODY(body, body_d) MEMCPY_DEVICE_TO_HOST(body, body_d, sizeof(liteBody))

# define EG_SET_BODY(body_d, body) MEMCPY_HOST_TO_DEVICE(body_d, body, sizeof(liteBody))

# define EG_GET_MODEL(model, model_d) MEMCPY_DEVICE_TO_HOST(model, model_d, sizeof(liteModel))

# define EG_SET_MODEL(model_d, model) MEMCPY_HOST_TO_DEVICE(model_d, model, sizeof(liteModel))

# define EG_GET_HEADERS_AT(header, headers_d_start, nmemb) MEMCPY_DEVICE_TO_HOST(header, headers_d_start, nmemb*sizeof(int)); 

# define EG_GET_DATA_AT(data, datas_d_start, nmemb) MEMCPY_DEVICE_TO_HOST(data, datas_d_start, nmemb*sizeof(double)); 


#else	/* Witout NVCC */


# define EG_NEW(ptraddr, type, nmemb) *((type**)ptraddr) = (type*) EG_alloc(nmemb*sizeof(type))

# define EG_REALLOC(dstaddr, src, type, old_nmemb, nmemb) *((type**)dstaddr) = (type*) EG_reall(src, nmemb*sizeof(type))

# define EG_FREE(ptr) EG_free(ptr)

# define EG_COPY(ptr_d, ptr, type, nmemb) memcpy((void*)ptr_d, (const void*)ptr, nmemb*sizeof(type))

# define EG_GET(ptr, ptr_d, type, nmemb) memcpy((void*)ptr, (const void*)ptr_d, nmemb*sizeof(type))

# define EG_SET_STR(string_daddr, string) *((char**)string_daddr) = EG_strdup(string)

# define EG_GET_STRN(string, string_daddr, N) strncpy(string,*((char**)string_daddr),N)

# define EG_NEW_OBJECT(objectaddr) *((egObject**)objectaddr) = (egObject *) EG_alloc(sizeof(egObject))

# define EG_FREE_OBJECT(object) EG_free(object)

# define EG_GET_OBJECT_PTR(objectaddr, object_daddr) *((egObject**)objectaddr) = *((egObject**)object_daddr)

# define EG_GET_OBJECT(object, object_d) memcpy((void*)object, (const void*)(object_d), sizeof(egObject))

# define EG_SET_OBJECT_PTR(object_daddr, objectaddr) *((egObject**)object_daddr) = *((egObject**)objectaddr)

# define EG_SET_OBJECT(object_daddr, object) *(*((egObject**)object_daddr)) = *(object)

# define EG_NEW_CNTXT(cntxaddr) *((egCntxt**)cntxaddr) = (egCntxt *) EG_alloc(sizeof(egCntxt))

# define EG_FREE_CNTXT(cntx) EG_free(cntx)

# define EG_GET_CNTXT(cntx, cntx_d) memcpy((void*)cntx, (const void*)(cntx_d), sizeof(egCntxt))

# define EG_SET_CNTXT(cntx_d, cntx) *((egCntxt*)cntx_d) = *(cntx);

# define EG_GET_ATTRS_PTR(attrs, attrs_d) attrs = attrs_d

# define EG_SET_ATTRS_PTR(attrs_daddr, attrs) *((egAttrs**)attrs_daddr) = attrs

# define EG_GET_ATTRS(attrs, attrs_d) *(attrs) = *(attrs_d)

# define EG_SET_ATTRS(attrs_daddr, attrs) *((egAttrs*)attrs_daddr) = *(attrs)

# define EG_GET_ATTR_PTR(attr, attr_d) attr = attr_d

# define EG_GET_ATTR(attr, attr_daddr) *(attr) = *((egAttr*)attr_daddr)

# define EG_SET_ATTR(attr_daddr, attr) *((egAttr*)attr_daddr) = *(attr)

# define EG_GET_TESSEL(tess, tess_d) memcpy((void*)tess, (const void*)(tess_d), sizeof(egTessel))

# define EG_GET_TESS1D(tess1d, tess1d_daddr) memcpy((void*)tess1d, (const void*)(tess1d_daddr), sizeof(egTess1D))

# define EG_GET_TESS2D(tess2d, tess2d_daddr) memcpy((void*)tess2d, (const void*)(tess2d_daddr), sizeof(egTess2D))

# define EG_GET_PATCH(patch, patch_daddr) memcpy((void*)patch, (const void*)(patch_daddr), sizeof(egPatch))

# define EG_GET_GEOM(geom, geom_d) *(geom) = *(geom_d)

# define EG_SET_GEOM(geom_d, geom) *(geom_d) = *(geom)

# define EG_GET_NODE(node, node_d) *(node) = *(node_d)

# define EG_SET_NODE(node_d, node) *(node_d) = *(node)

# define EG_GET_EDGE(edge, edge_d) *(edge) = *(edge_d)

# define EG_SET_EDGE(edge_d, edge) *(edge_d) = *(edge)

# define EG_GET_LOOP(loop, loop_d) *(loop) = *(loop_d)

# define EG_SET_LOOP(loop_d, loop) *(loop_d) = *(loop)

# define EG_GET_FACE(face, face_d) *(face) = *(face_d)

# define EG_SET_FACE(face_d, face) *(face_d) = *(face)

# define EG_GET_SHELL(shell, shell_d) *(shell) = *(shell_d)

# define EG_SET_SHELL(shell_d, shell) *(shell_d) = *(shell)

# define EG_GET_BODY(body, body_d) *(body) = *(body_d)

# define EG_SET_BODY(body_d, body) *(body_d) = *(body)

# define EG_GET_MODEL(model, model_d) *(model) = *(model_d)

# define EG_SET_MODEL(model_d, model) *(model_d) = *(model)

# define EG_GET_HEADERS_AT(header, headers_d_start, nmemb) memcpy((void*)header, (const void*)(headers_d_start), nmemb*sizeof(int))

# define EG_GET_DATA_AT(data, datas_d_start, nmemb) memcpy((void*)data, (const void*)(datas_d_start), nmemb*sizeof(double))

#endif

#endif
