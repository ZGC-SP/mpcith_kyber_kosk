
#ifndef NTL_config__H
#define NTL_config__H

/*************************************************************************

                          NTL Configuration File
                          ----------------------

This file is automatically generated by the configure script.

You can also edit this file by hand, but that is not generally recommended.

To set a flag, just replace the pre-processor directive 
'if 0' by 'if 1' for that flag, which causes the appropriate macro 
to be defined.  Of course,  to unset a flag, just replace the 
'if 1' by an 'if 0'.

 *************************************************************************/



/*************************************************************************
 *
 * Basic Configuration Options
 *
 *************************************************************************/


 /* None of these flags are set by the configuration wizard;
  * they must be set by hand, before installation begins.
  */


#if 0
#define NTL_LEGACY_NO_NAMESPACE

/* 
 * By default, NTL components are declared inside the namespace NTL.
 * Set this flag if you want to instead have these components
 * declared in the global namespace.  This is for backward
 * compatibility only -- not recommended.
 *
 */

#endif


#if 0
#define NTL_LEGACY_INPUT_ERROR

/*
 * Also for backward compatibility.  Set if you want input 
 * operations to abort on error, instead of just setting the
 * "fail bit" of the input stream.
 *
 */


#endif

#if 1
#define NTL_TLS_HACK

/* Set if you want to compile NTL with "TLS hack"
 *
 */

#endif

#if 1
#define NTL_THREADS

/* Set if you want to compile NTL as a thread-safe library.
 *
 */

#endif


#if 0
#define NTL_EXCEPTIONS

/* Set if you want to compile NTL with exceptions enabled
 *
 */

#endif

#if 1
#define NTL_THREAD_BOOST

/* Set if you want to compile NTL to exploit threads internally.
 *
 */

#endif


#if 0
#define NTL_GMP_LIP

/* 
 * Use this flag if you want to use GMP as the long integer package.
 * This can result in significantly faster code on some platforms.
 * It requires that the GMP package (version >= 3.1) has already been
 * installed.  You will also have to set the variables GMP_OPT_INCDIR,
 * GMP_OPT_LIBDIR, GMP_OPT_LIB in the makefile (these are set automatically
 * by the confiuration script when you pass the flag NTL_GMP_LIP=on
 * to that script.
 *
 * Beware that setting this flag can break some very old NTL codes.
 *
 * You may also have to edit the makefile to modify the variables
 * GMP_OPT_INCDIR, GMP_OPT_LIBDIR, and GMP_OPT_LIB.
 */

#endif

#if 0
#define NTL_GF2X_LIB

/* 
 * Use this flag if you want to use the gf2x library for
 * faster GF2X arithmetic.
 * This can result in significantly faster code, especially
 * when working with polynomials of huge degree.
 * You will also have to set the variables GF2X_OPT_INCDIR,
 * GF2X_OPT_LIBDIR, GF2X_OPT_LIB in the makefile (these are set automatically
 * by the confiuration script when you pass the flag NTL_GF2X_LIB=on
 * to that script.
 *
 * You may also have to edit the makefile to modify the variables
 * GF2X_OPT_INCDIR, GF2X_OPT_LIBDIR, and GF2X_OPT_LIB.
 */

#endif


#if 1
#define NTL_STD_CXX11

/*
 * Set this flag if you want to enable C++11 features within NTL.
 */

#endif

#if 0
#define NTL_STD_CXX14

/*
 * Set this flag if you want to enable C++14 features within NTL.
 */

#endif

#if 1
#define NTL_DISABLE_MOVE_ASSIGN

/*
 * Set this flag if you want to disable move assignment
 * operators for vectors (and, by extension, polynomials)
 * and matrices.
 */

#endif

#if 0
#define NTL_DISABLE_MOVE

/*
 * This flag disables all move constructors and assignments. 
 */

#endif


#if 0
#define NTL_UNSIGNED_LONG_LONG_TYPE unsigned long long

/*
 *   NTL_UNSIGNED_LONG_LONG_TYPE will be used
 *   to declare 'double word' unsigned integer types.
 *   If left undefined, some "ifdef magic" will attempt
 *   to find the best choice for your platform, depending
 *   on the compiler and wordsize.  On 32-bit machines,
 *   this is usually 'unsigned long long'.
 *
 */

#endif


#if 0
#define NTL_CLEAN_INT

/*
 *   This will disallow the use of some non-standard integer arithmetic
 *   that may improve performance somewhat.
 *
 */

#endif

#if 1
#define NTL_CLEAN_PTR

/*
 *   This will disallow the use of some non-standard pointer arithmetic
 *   that may improve performance somewhat.
 *
 */

#endif

#if 1
#define NTL_SAFE_VECTORS

/*
 * This will compile NTL in "safe vector" mode, only assuming
 * the relocatability property for trivial types and types
 * explicitly declared relocatable.  See vector.txt for more details.
 */

#endif

#if 0
#define NTL_ENABLE_AVX_FFT

/*
 * This will compile NTL in a way that enables an AVX implemention
 * of the small-prime FFT.
 */

#endif


#if 0
#define NTL_AVOID_AVX512

/*
 * This will compile NTL in a way that avoids 512-bit operations,
 * even if AVX512 is available.
 */

#endif
 
#if 0
#define NTL_RANGE_CHECK

/*
 *   This will generate vector subscript range-check code.
 *   Useful for debugging, but it slows things down of course.
 *
 */

#endif





#if 1
#define NTL_NO_INIT_TRANS

/*
 *   Without this flag, NTL uses a special code sequence to avoid
 *   copying large objects in return statements.  However, if your
 *   compiler optimizes away the return of a *named* local object,
 *   this is not necessary, and setting this flag will result
 *   in *slightly* more compact and efficient code.  Although
 *   the emeriging C++ standard allows compilers to perform
 *   this optimization, I know of none that currently do.
 *   Most will avoid copying *temporary* objects in return statements,
 *   and NTL's default code sequence exploits this fact.
 *
 */

#endif


#if 0
#define NTL_X86_FIX

/*
 *  Forces the "x86 floating point fix", overriding the default behavior.
 *  By default, NTL will apply the "fix" if it looks like it is
 *  necessary, and if knows how to fix it.
 *  The problem addressed here is that x86 processors sometimes
 *  run in a mode where FP registers have more precision than doubles.
 *  This will cause code in quad_float.cpp some trouble.
 *  NTL can normally correctly detect the problem, and fix it,
 *  so you shouldn't need to worry about this or the next flag.
 *  
 */

#elif 0
#define NTL_NO_X86_FIX
/*
 *  Forces no "x86 floating point fix", overriding the default behavior.
 */

#endif



#if 0
#define NTL_LEGACY_SP_MULMOD

/* Forces legacy single-precision MulMod implementation.
 */

#endif


#if 0
#define NTL_DISABLE_LONGDOUBLE

/* Explicitly disables us of long double arithmetic
 */

#endif


#if 0
#define NTL_DISABLE_LONGLONG

/* Explicitly disables us of long long arithmetic 
 */

#endif

#if 0
#define NTL_DISABLE_LL_ASM

/* Explicitly disables us of inline assembly as a replacement
 * for long lobg arithmetic.
 */

#endif


#if 0
#define NTL_MAXIMIZE_SP_NBITS

/* Allows for 62-bit single-precision moduli on 64-bit platforms.
 * By default, such moduli are restricted to 60 bits, which
 * usually gives slightly better performance across a range of
 * of parameters.
 */

#endif

/*************************************************************************
 *
 *  Performance Options
 *
 *************************************************************************/



/* There are three strategies to implmement single-precision
 * modular multiplication with preconditioning (see the MulModPrecon
 * function in the ZZ module): the default and NTL_SPMM_ULL.
 * This plays a crucial role in the  "small prime FFT" used to 
 * implement polynomial arithmetic, and in other CRT-based methods 
 * (such as linear  algebra over ZZ), as well as polynomial and matrix 
 * arithmetic over zz_p.  
 */



#if 1
#define NTL_SPMM_ULL

/*    This also causes an "all integer"
 *    implementation of MulModPrecon to be used.
 *    It us usually a faster implementation,
 *    but it is not enturely portable.
 *    It relies on double-word unsigned multiplication
 *    (see NTL_UNSIGNED_LONG_LONG_TYPE above). 
 *
 */


#endif



/*
 * The following two flags provide additional control for how the 
 * FFT modulo single-precision primes is implemented.
 */

#if 1
#define NTL_FFT_BIGTAB

/*
 * Precomputed tables are used to store all the roots of unity
 * used in FFT computations. 
 *
 */


#endif


#if 1
#define  NTL_FFT_LAZYMUL

/*
 * When set, a "lazy multiplication" strategy due to David Harvey:
 * see his paper "FASTER ARITHMETIC FOR NUMBER-THEORETIC TRANSFORMS".
 *
 */


#endif



#if 1
#define NTL_AVOID_BRANCHING

/*
 *   With this option, branches are replaced at several 
 *   key points with equivalent code using shifts and masks.
 *   It may speed things up on machines with 
 *   deep pipelines and high branch penalities.
 *   This flag mainly affects the implementation of the
 *   single-precision modular arithmetic routines.
 *
 */

#endif



#if 1
#define NTL_TBL_REM

/*
 *
 *   With this flag, some divisions are avoided in the
 *   ZZ_pX multiplication routines.  
 *
 */

#endif



#if 1
#define NTL_CRT_ALTCODE

/*
 * Employs an alternative CRT strategy.
 * Only relevant with GMP.
 * Seems to be marginally faster on some x86_64 platforms.
 *
 */

#endif

#if 0
#define NTL_CRT_ALTCODE_SMALL

/*
 * Employs an alternative CRT strategy for small moduli.
 * Only relevant with GMP.
 * Seems to be marginally faster on some x86_64 platforms.
 *
 */

#endif


#if 0
#define NTL_GF2X_ALTCODE

/*
 * With this option, the default strategy for implmenting low-level
 * GF2X multiplication is replaced with an alternative strategy.
 * This alternative strategy seems to work better on RISC machines
 * with deep pipelines and high branch penalties (like a powerpc),
 * but does no better (or even worse) on x86s.
 *
 */

#elif 1
#define NTL_GF2X_ALTCODE1


/*
 * Yest another alternative strategy for implementing GF2X
 * multiplication.
 *
 */


#endif

#if 0
#define NTL_GF2X_NOINLINE

/*
 * By default, the low-level GF2X multiplication routine in inlined.
 * This can potentially lead to some trouble on some platforms,
 * and you can override the default by setting this flag.
 *
 */

#endif

#if 0
#define NTL_RANDOM_AES256CTR

/*
 * By default, the random-number generator is based on ChaCha20.
 * From a performance perspective, this choice may not be optimal
 * for platforms featuring AES hardware support.
 * By setting this flag you can override the default and use an
 * AES-256-CTR based random-number generator.
 *
 */

#endif


/* sanity checks */

#if (defined(NTL_THREAD_BOOST) && !defined(NTL_THREADS))
#error "NTL_THREAD_BOOST defined but not NTL_THREADS"
#endif


#if (defined(NTL_THREADS) && !(defined(NTL_STD_CXX11) || defined(NTL_STD_CXX14)))
#error "NTL_THREADS defined but not NTL_STD_CXX11 or NTL_STD_CXX14"
#endif


#if (defined(NTL_EXCEPTIONS) && !(defined(NTL_STD_CXX11) || defined(NTL_STD_CXX14)))
#error "NTL_EXCEPTIONS defined but not NTL_STD_CXX11 or NTL_STD_CXX14"
#endif


#if (defined(NTL_SAFE_VECTORS) && !(defined(NTL_STD_CXX11) || defined(NTL_STD_CXX14)))
#error "NTL_SAFE_VECTORS defined but not NTL_STD_CXX11 or NTL_STD_CXX14"
#endif





#endif
