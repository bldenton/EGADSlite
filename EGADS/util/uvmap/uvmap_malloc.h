void uvmap_register_ext_free (void (*ext_free_routine) (void *ptr));

void uvmap_register_ext_malloc (void * (*ext_malloc_routine) (INT_ *err_flag, size_t size));

void uvmap_register_ext_realloc (void * (*ext_realloc_routine) (INT_ *err_flag, void *ptr, size_t size));

void uvmap_free (void *ptr);

void * uvmap_malloc (INT_ *err_flag, size_t size);

void * uvmap_realloc (INT_ *err_flag, void *ptr, size_t size);
