INT_ uvmap_struct_add_entry (
  INT_ idef,
  INT_ nbface,
  INT_ *idibf,
  INT_3D *inibf,
  INT_3D *ibfibf,
  DOUBLE_2D *u,
  uvmap_struct **uvmap_struct_ptr);

INT_ uvmap_struct_get_entry (
  INT_ idef,
  INT_ *index,
  INT_ *isrch,
  INT_ *ibface,
  INT_ *nbface,
  INT_ **idibf,
  INT_ **msrch,
  INT_3D **inibf,
  INT_3D **ibfibf,
  DOUBLE_2D **u,
  uvmap_struct *uvmap_struct_ptr);

void uvmap_struct_set_srch_data (
  INT_ index,
  INT_ isrch,
  INT_ ibface,
  uvmap_struct *uvmap_struct_ptr);

INT_ uvmap_struct_find_entry (
  INT_ idef,
  INT_ *index,
  uvmap_struct *uvmap_struct_ptr);

void uvmap_struct_free (void *ptr);

void uvmap_struct_free_idef (INT_ idef, uvmap_struct *uvmap_struct_ptr);

void uvmap_struct_free_index (INT_ index, uvmap_struct *uvmap_struct_ptr);
