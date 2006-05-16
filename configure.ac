AC_INIT(src/ica.c)
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
FLIBS=`${R_HOME}/bin/R CMD config FLIBS`
BLAS_LIBS=`${R_HOME}/bin/R CMD config BLAS_LIBS`
LAPACK_LIBS=`${R_HOME}/bin/R CMD config LAPACK_LIBS`

# dummy
AC_CHECK_HEADERS(math.h)

if test -n "`grep HAVE_F77_UNDERSCORE ${R_HOME}/include/Rconfig.h | \
             grep '#define'`"; then
  have_f77_underscore=true
elif test -n "`grep F77 ${R_HOME}/include/Rconfig.h | grep '_$'`"; then
  have_f77_underscore=true
fi
if ${have_f77_underscore}; then
  scopy_func=scopy_
  lsame_func=lsame_
  sgesdd_func=sgesdd_
  xerbla=xerbla_
else
  scopy_func=scopy
  lsame_func=lsame
  sgesdd_func=sgesdd
  xerbla=xerbla
fi

if test -z "${BLAS_LIBS}"; then
  echo "it seems you are not using a BLAS library: using our BLAS and LAPACK routines"
  (cd src; ln -s sblas.unused sblas.f; ln -s lsame.unused lsame.f; ln -s slapack.unused slapack.f; ln -s slamc.unused slamc.f)
else
  LIBS="${BLAS_LIBS} ${FLIBS}"
  AC_MSG_CHECKING([for scopy in ${BLAS_LIBS}])
  AC_TRY_LINK([void ${xerbla}(char *srname, int *info){}], ${scopy_func}(),
      [blas_ok=yes], [blas_ok=no])
  AC_MSG_RESULT([${blas_ok}])
  if test "${blas_ok}" = no; then
    AC_MSG_WARN([could not find scopy in your BLAS library, using our BLAS and LAPACK routines])
    (cd src; ln -s sblas.unused sblas.f; ln -s lsame.unused lsame.f; ln -s slapack.unused slapack.f; ln -s slamc.unused slamc.f)
  else
    LIBS="${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS}"
    AC_CHECK_FUNC(${sgesdd_func}, BLAS_LIBS="${LAPACK_LIBS} ${BLAS_LIBS}", 
                  NO_LIBLAPACK=1)
    if test ${NO_LIBLAPACK}; then
      AC_MSG_WARN([could not find sgesdd: using our LAPACK routines])
      (cd src; ln -s slapack.unused slapack.f; ln -s slamc.unused slamc.f)
    fi
    AC_CHECK_FUNC(${lsame_func}, , NO_LSAME=1)
    if test ${NO_LSAME}; then
      AC_MSG_WARN([could not find lsame in your BLAS or LAPACK library, using our copy])
      (cd src; ln -s lsame.unused lsame.f)
    fi
  fi
fi


AC_SUBST(BLAS_LIBS)
AC_OUTPUT(src/Makevars)
