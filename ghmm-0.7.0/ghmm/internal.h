#ifndef INTERNAL_H
#define INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

#ifndef DECLARE_DEPRECATED
#  if ( __GNUC__ == 3 && __GNUC_MINOR__ > 0 ) || __GNUC__ > 3 
#    define DECLARE_DEPRECATED __attribute__((deprecated))
#  else
#    define DECLARE_DEPRECATED /* empty, no compile-time warnings */
#  endif
#endif

#ifndef ARRAY_MALLOC
#define ARRAY_MALLOC(ptr, entries) { if (!((ptr) = mes_malloc(sizeof(*(ptr)) * (entries)))) \
                                         {mes_proc(); goto STOP;}                           \
				   }
#endif

#ifndef ARRAY_CALLOC
#define ARRAY_CALLOC(ptr, entries) { if (!((ptr) = mes_calloc(sizeof(*(ptr)) * (entries))))   \
                                         {mes_proc(); goto STOP;}                             \
				   }
#endif

#ifndef ARRAY_REALLOC
#define ARRAY_REALLOC(ptr, entries) { if (mes_realloc ((void**)&(ptr), sizeof(*(ptr))*(entries))) \
                                          {mes_proc(); goto STOP;}                                \
				    }
#endif


#ifdef __cplusplus
}
#endif /*__cplusplus*/

#endif
