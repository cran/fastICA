AC_INIT(src/ica.c)
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi

CC=`${R_HOME}/bin/R CMD config CC`
CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
FLIBS=`${R_HOME}/bin/R CMD config FLIBS`
LIBS=`grep "^LIBS" ${R_HOME}/etc/Makeconf | sed "s/^LIBS *= *//"`
LIBS0=${LIBS}
BLAS_LIBS=`${R_HOME}/bin/R CMD config BLAS_LIBS`
LAPACK_LIBS=`${R_HOME}/bin/R CMD config LAPACK_LIBS`

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
else
  scopy_func=scopy
  lsame_func=lsame
  sgesdd_func=sgesdd
fi

if test -z "${BLAS_LIBS}"; then
  echo "it seems you are not using a BLAS library: using our BLAS routines"
  (cd src; ln -s sblas.unused sblas.f; ln -s lsame.unused lsame.f; ln -s slapack.unused slapack.f)
else
  LIBS="${BLAS_LIBS} ${FLIBS} ${LIBS}"
  AC_CHECK_FUNC(${scopy_func}, , NO_SBLAS=1)
  if test ${NO_SBLAS}; then
    echo "could not find SCOPY in your BLAS library, using our BLAS routines"
    (cd src; ln -s sblas.unused sblas.f; ln -s lsame.unused lsame.f; ln -s slapack.unused slapack.f)
  else
    LIBS="${LAPACK_LIBS} ${BLAS_LIBS} ${FLIBS} ${LIBS0}"
    AC_CHECK_FUNC(${sgesdd_func}, BLAS_LIBS="${LAPACK_LIBS} ${BLAS_LIBS}", 
                  NO_LIBLAPACK=1)
    if test ${NO_LIBLAPACK}; then
      echo "could not find SGESDD: using our LAPACK routines"
      (cd src; ln -s slapack.unused slapack.f)
    fi
    AC_CHECK_FUNC(${lsame_func}, , NO_LSAME=1)
    if test ${NO_LSAME}; then
      echo "could not find LSAME in your BLAS or LAPACK library, using our copy"
      (cd src; ln -s lsame.unused lsame.f)
    fi
  fi
fi


AC_SUBST(BLAS_LIBS)
AC_OUTPUT(src/Makevars)
