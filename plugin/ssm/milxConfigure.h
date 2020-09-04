/*
 * here is where system computed values get stored these values should only
 * change when the target compile platform changes
 */

/* what threading system are we using */
/* #undef CMAKE_USE_PTHREADS */
#ifdef CMAKE_USE_PTHREADS
#define MILX_USE_PTHREADS
#endif

/* #undef CMAKE_USE_SPROC */
#ifdef CMAKE_USE_SPROC
#define MILX_USE_SPROC
#endif

/* #undef CMAKE_HP_PTHREADS */
#ifdef CMAKE_HP_PTHREADS
#define MILX_HP_PTHREADS
#endif

/* #undef CMAKE_USE_WIN32_THREADS */
#ifdef CMAKE_USE_WIN32_THREADS
#define MILX_USE_WIN32_THREADS
#endif

#define BUILD_SHARED_LIBS
#ifdef BUILD_SHARED_LIBS
#define MILXDLL
#else
#define MILXSTATIC
#endif

/* #undef CMAKE_NO_STD_NAMESPACE */
/* #undef CMAKE_NO_ANSI_STREAM_HEADERS */
/* #undef CMAKE_NO_ANSI_STRING_STREAM */
/* #undef CMAKE_NO_ANSI_FOR_SCOPE */

#define MILX_VERSION_MAJOR 
#define MILX_VERSION_MINOR 
#define MILX_VERSION_PATCH 
#define MILX_VERSION_STRING ""

/* #undef USE_ITK181 */
/* #undef USE_ITK20 */
/* #undef USE_ITK22 */
/* #undef USE_ITK26 */
/* #undef USE_ITK28 */

