#ifndef __UVMAP_STRUCT_H__

#define __UVMAP_STRUCT_H__

typedef struct _uvmap_struct uvmap_struct;

struct _uvmap_struct
{
  INT_ ndef;
  INT_ mdef;
  INT_ idef;
  INT_ isrch;
  INT_ ibface;
  INT_ nbface;
  INT_ nnode;
  INT_ *idibf;
  INT_ *msrch;
  INT_3D *ibfibf;
  INT_3D *inibf;
  DOUBLE_2D *u;
};

#endif
