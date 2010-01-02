#include "stdio.h"
#include "R.h"
#include "Rinternals.h"

SEXP pgm( SEXP filename ) {
  SEXP pgm, class, names, dim;
  int r, c, h, w, m, *d;
  char * b = R_alloc( 64, sizeof( char ) );
  FILE * f = fopen( CHAR(STRING_ELT(filename,0)), "rb" );
  if( f == NULL ) { error( "file open failed" ); }
  PROTECT(pgm = allocVector(VECSXP, 4));
  PROTECT(class = allocVector(STRSXP, 1));
  PROTECT(names = allocVector(STRSXP, 4));
  SET_STRING_ELT(class, 0, mkChar("pgm"));
  SET_STRING_ELT(names, 0, mkChar("width"));
  SET_STRING_ELT(names, 1, mkChar("height"));
  SET_STRING_ELT(names, 2, mkChar("maxval"));
  SET_STRING_ELT(names, 3, mkChar("data"));
  setAttrib(pgm, R_NamesSymbol, names);
  setAttrib(pgm, R_ClassSymbol, class);

  SET_VECTOR_ELT(pgm, 0, allocVector(INTSXP, 1));
  SET_VECTOR_ELT(pgm, 1, allocVector(INTSXP, 1));
  SET_VECTOR_ELT(pgm, 2, allocVector(INTSXP, 1));

  //0. check for magic number (ascii 'P5')
  fread( b, 1, 2, f );
  if( b[0] != 'P' || b[1] != '5' ) { error( "file not PGM" ); }
  //1. read width, height, maxval, set dim, alloc data
  fscanf( f, "%d", INTEGER(VECTOR_ELT(pgm, 0)) );
  w = INTEGER(VECTOR_ELT(pgm, 0))[0];
  fscanf( f, "%d", INTEGER(VECTOR_ELT(pgm, 1)) );
  h = INTEGER(VECTOR_ELT(pgm, 1))[0];
  fscanf( f, "%d", INTEGER(VECTOR_ELT(pgm, 2)) );
  m = INTEGER(VECTOR_ELT(pgm, 2))[0];
  SET_VECTOR_ELT(pgm, 3, allocVector(INTSXP, h*w));
  d = INTEGER(VECTOR_ELT(pgm, 3));
  PROTECT(dim = allocVector(INTSXP, 2));
  INTEGER(dim)[0] = h;
  INTEGER(dim)[1] = w;
  setAttrib(VECTOR_ELT(pgm, 3), R_DimSymbol, dim);
  //2. get the whitespace character, read data
  fgetc( f );
  for(r = 0; r < h; r++) {
    for(c = 0; c < w; c++ ) {
      d[c*h+r] = 0;
      if(m>=256) { fread( ((void *)(d+(c*h+r)))+1, 1, 1, f ); }
      fread( ((void *)(d+(c*h+r))), 1, 1, f );
      if( ferror( f ) ) { error("file read error"); }
    }
  }
  fclose(f);
  UNPROTECT(4);
  return pgm;
}
  
  
