//------------------------------------------------------------------------------
// SPEX_Update/Test/test.c: support functions for the test programs
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "test.h"

//------------------------------------------------------------------------------
// SPEX_mmread:  read a matrix from a Matrix Market file, modified from LAGraph
//------------------------------------------------------------------------------

// The file format used here is compatible with all variations of the Matrix
// Market "coordinate" and "array" format (http://www.nist.gov/MatrixMarket).
// The format is fully described in SPEX/Doc/MatrixMarket.pdf, and
// summarized here (with extensions for SPEX).

// The first line of the file starts with %%MatrixMarket, with the following
// format:

//      %%MatrixMarket matrix <fmt> <type> <storage>

//      <fmt> is one of: coordinate or array.  The former is a sparse matrix in
//      triplet form.  The latter is a dense matrix in column-major form.
//      Both formats are returned as a SPEX_Matrix.

//      <type> is one of: real, complex, pattern, or integer.  The real,
//      integer, and pattern formats are returned as SPEX_FP64, SPEX_INT64, and
//      SPEX_BOOL, respectively, but these types are modified the %GraphBLAS
//      structured comment described below.  Complex matrices are returned
//      using the SPEX_ComplexFP64 type (which is a GraphBLAS type corresponding
//      to the ANSI C11 double complex type).

//      <storage> is one of: general, Hermitian, symmetric, or skew-symmetric.
//      The Matrix Market format is case-insensitive, so "hermitian" and
//      "Hermitian" are treated the same).

//      Not all combinations are permitted.  Only the following are meaningful:

//      (1) (coordinate or array) x (real, integer, or complex)
//          x (general, symmetric, or skew-symmetric)

//      (2) (coordinate or array) x (complex) x (Hermitian)

//      (3) (coodinate) x (pattern) x (general or symmetric)

// Any other lines starting with "%" are treated as comments, and are ignored.
// Comments may be interspersed throughout the file.  Blank lines are ignored.
// The Matrix Market header is optional in this routine (it is not optional in
// the Matrix Market format).  If not present, the <fmt> defaults to
// coordinate, <type> defaults to real, and <storage> defaults to general.  The
// remaining lines are space delimited, and free format (one or more spaces can
// appear, and each field has arbitrary width).

// The Matrix Market file <fmt> can be coordinate or array:

//      coordinate:  for this format, the first non-comment line must appear,
//          and it must contain three integers:

//              nrows ncols nvals

//          For example, a 5-by-12 matrix with 42 entries would have:

//              5 12 42

//          Each of the remaining lines defines one entry.  The order is
//          arbitrary.  If the Matrix Market <type> is real or integer, each
//          line contains three numbers: row index, column index, and value.
//          For example, if A(3,4) is equal to 5.77, a line:

//              3 4 5.77

//          would appear in the file.  The indices in the Matrix Market are
//          1-based, so this entry becomes A(2,3) in the SPEX_Matrix returned to
//          the caller.  If the <type> is pattern, then only the row and column
//          index appears.  If <type> is complex, four values appear.  If
//          A(8,4) has a real part of 6.2 and an imaginary part of -9.3, then
//          the line is:

//              8 4 6.2 -9.3

//          and since the file is 1-based but a GraphBLAS matrix is always
//          0-based, one is subtracted from the row and column indices in the
//          file, so this entry becomes A(7,3).

//      array: for this format, the first non-comment line must appear, and
//          it must contain just two integers:

//              nrows ncols

//          A 5-by-12 matrix would thus have the line

//              5 12

//          Each of the remaining lines defines one entry, in column major
//          order.  If the <type> is real or integer, this is the value of the
//          entry.  An entry if <type> of complex consists of two values, the
//          real and imaginary part.  The <type> cannot be pattern in this
//          case.

//      For both coordinate and array formats, real and complex values may use
//      the terms INF, +INF, -INF, and NAN to represent floating-point infinity
//      and NaN values.

// The <storage> token is general, symmetric, skew-symmetric, or Hermitian:

//      general: the matrix has no symmetry properties (or at least none
//          that were exploited when the file was created).

//      symmetric:  A(i,j) == A(j,i).  Only entries on or below the diagonal
//          appear in the file.  Each off-diagonal entry in the file creates
//          two entries in the SPEX_Matrix that is returned.

//      skew-symmetric:  A(i,j) == -A(i,j).  There are no entries on the
//          diagonal.  Only entries below the diagonal appear in the file.
//          Each off-diagonal entry in the file creates two entries in the
//          SPEX_Matrix that is returned.

//      Hermitian: square complex matrix with A(i,j) = conj (A(j,i)).
//          All entries on the diagonal are real.  Each off-diagonal entry in
//          the file creates two entries in the SPEX_Matrix that is returned.
//          NOTE: THIS IS NOT HANDLED BY THIS FUNCTION, SINCE IT IS MEANINGLESS
//                FOR LINEAR PROGRAMS

// According to the Matrix Market format, entries are always listed in
// column-major order.  This rule is follwed by SPEX_mmwrite.  However,
// SPEX_mmread can read the entries in any order.

// Parts of this code are from SuiteSparse/CHOLMOD/Check/cholmod_read.c, and
// are used here by permission of the author of CHOLMOD/Check (T. A. Davis).

//------------------------------------------------------------------------------
// get_line
//------------------------------------------------------------------------------

// Read one line of the file, return true if successful, false if EOF.
// The string is returned in buf, converted to lower case.

static inline bool get_line
(
    FILE *f,        // file open for reading
    char *buf       // size MAXLINE+1
)
{

    // check inputs
    ASSERT (f != NULL) ;
    ASSERT (buf != NULL) ;

    // read the line from the file
    buf [0] = '\0' ;
    buf [1] = '\0' ;
    if (fgets (buf, MAXLINE, f) == NULL)
    {
        // EOF or other I/O error
        return (false) ;
    }
    buf [MAXLINE] = '\0' ;

    // convert the string to lower case
    for (int k = 0 ; k < MAXLINE && buf [k] != '\0' ; k++)
    {
        buf [k] = tolower (buf [k]) ;
    }
    return (true) ;
}

//------------------------------------------------------------------------------
// is_blank_line
//------------------------------------------------------------------------------

// returns true if buf is a blank line or comment, false otherwise.

static inline bool is_blank_line
(
    char *buf       // size MAXLINE+1, never NULL
)
{

    // check inputs
    ASSERT (buf != NULL) ;

    // check if comment line
    if (buf [0] == '%')
    {
        // line is a comment
        return (true) ;
    }

    // check if blank line
    for (int k = 0 ; k <= MAXLINE ; k++)
    {
        int c = buf [k] ;
        if (c == '\0')
        {
            // end of line
            break ;
        }
        if (!isspace (c))
        {
            // non-space character; this is not an error
            return (false) ;
        }
    }

    // line is blank
    return (true) ;
}

//------------------------------------------------------------------------------
// read_double
//------------------------------------------------------------------------------

// Read a single double value from a string.  The string may be any string
// recognized by sscanf, or inf, -inf, +inf, or nan.  The token infinity is
// also OK instead of inf (only the first 3 letters of inf* or nan* are
// significant, and the rest are ignored).

static inline bool read_double      // true if successful, false if failure
(
    char *p,        // string containing the value
    double *rval    // value to read in
)
{
    while (*p && isspace (*p)) p++ ;   // skip any spaces

    if ((strncmp (p, "inf", 3) == 0) || (strncmp (p, "+inf", 4) == 0))
    {
        (*rval) = INFINITY ;
    }
    else if (strncmp (p, "-inf", 4) == 0)
    {
        (*rval) = -INFINITY ;
    }
    else if (strncmp (p, "nan", 3) == 0)
    {
        (*rval) = NAN ;
    }
    else
    {
        if (sscanf (p, "%lg", rval) != 1)
        {
            // bad file format, EOF, or other I/O error
            return (false) ;
        }
    }
    return (true) ;
}

//------------------------------------------------------------------------------
// read_entry
//------------------------------------------------------------------------------

static inline bool read_entry   // true if successful, false if failure
(
    char *p,        // string containing the value
    SPEX_type type, // type of value to read, either SPEX_INT64 or SPEX_FP64
    bool pattern,   // if true, then the value is 1
    char *x         // value read in, a pointer to space of size of the type
)
{

    int64_t ival = 1 ;
    double rval = 1;

    while (*p && isspace (*p)) p++ ;   // skip any spaces

    // printf ("read entry [%s]: ", p) ;

    if (type == SPEX_INT64)
    {
        if (!pattern && sscanf (p, "%" SCNd64, &ival) != 1) return (false) ;
        // printf ("%" PRId64 "\n", ival) ;
        int64_t *result = (int64_t *) x ;
        result [0] = (int64_t) ival ;
    }
    else if (type == SPEX_FP64)
    {
        if (!pattern && !read_double (p, &rval)) return (false) ;
        // printf ("%g\n", rval) ;
        double *result = (double *) x ;
        result [0] = rval ;
    }
    else
    {
        // type not supported
        printf ("SPEX_mmread: read_entry: type not supported\n") ;
        return (false) ;
    }

    return (true) ;
}

//------------------------------------------------------------------------------
// negate_scalar: negate a scalar value
//------------------------------------------------------------------------------

// negate the scalar x.  Do nothing for bool or uint*.

static inline void negate_scalar
(
    SPEX_type type,
    void *x
)
{
    if (type == SPEX_INT64)
    {
        int64_t *value = (int64_t *) x ;
        (*value) = - (*value) ;
    }
    else if (type == SPEX_FP64)
    {
        double *value = (double *) x ;
        (*value) = - (*value) ;
    }
}

//------------------------------------------------------------------------------
// set_value
//------------------------------------------------------------------------------

// A(i,j) = x using SPEX_Matrix_setElement_<type>.  No typecasting is done.

static inline SPEX_info set_value
(
    SPEX_matrix A,
    SPEX_type type,
    int64_t i,
    int64_t j,
    char *x
)
{
    int64_t Anz = A->nz++;
    A->i[Anz] = i;
    A->j[Anz] = j;

    if (type == SPEX_INT64)
    {
        int64_t *value = (int64_t *) x ;
        A->x.int64[Anz] = *value;
        return SPEX_OK;
    }
    else if (type == SPEX_FP64)
    {
        double *value = (double *) x ;
        A->x.fp64[Anz] = *value;
        return SPEX_OK;
    }
    else
    {
        // type not supported
        return (SPEX_PANIC) ;
    }
}

//------------------------------------------------------------------------------
// SPEX_mmread
//------------------------------------------------------------------------------

SPEX_info SPEX_mmread
(
    SPEX_matrix *A_handle,// handle of matrix to create
    FILE *f,             // file to read from, already open
    SPEX_options option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (A_handle == NULL || f == NULL)
    {
        // input arguments invalid
        printf ("SPEX_mmread: bad args\n") ;
        return (SPEX_INCORRECT_INPUT) ;
    }
    *A_handle = NULL;
    SPEX_matrix A = NULL;

    //--------------------------------------------------------------------------
    // set the default properties
    //--------------------------------------------------------------------------

    MM_fmt_enum     MM_fmt     = MM_coordinate ;
    MM_type_enum    MM_type    = MM_real ;
    MM_storage_enum MM_storage = MM_general ;
    SPEX_type type = SPEX_FP64 ;
    uint64_t nrows = 0,   ncols = 0, nvals = 0 ;

    //--------------------------------------------------------------------------
    // read the Matrix Market header
    //--------------------------------------------------------------------------

    // Read the header.  This consists of zero or more comment lines (blank, or
    // starting with a "%" in the first column), followed by a single data line
    // containing two or three numerical values.  The first line is normally:
    //
    //          %%MatrixMarket matrix <fmt> <type> <storage>
    //
    // but this is optional.
    // If the %%MatrixMarket line is not present, then the <fmt> <type> and
    // <storage> are implicit.  If the first data line contains 3 items,
    // then the implicit header is:
    //
    //          %%MatrixMarket matrix coordinate real general
    //
    // If the first data line contains 2 items (nrows ncols), then the implicit
    // header is:
    //
    //          %%MatrixMarket matrix array real general
    //
    // The implicit header is an extension of the Matrix Market format.

    char buf [MAXLINE+1] ;

    bool got_mm_header = false ;

    for (int64_t line = 1 ; get_line (f, buf) ; line++)
    {

        //----------------------------------------------------------------------
        // parse the line
        //----------------------------------------------------------------------

        if ((line == 1) && (strncmp (buf, "%%matrixmarket", 14) == 0))
        {

            //------------------------------------------------------------------
            // read a Matrix Market header
            //------------------------------------------------------------------

            //  %%MatrixMarket matrix <fmt> <type> <storage>
            //  if present, it must be the first line in the file.

            got_mm_header = true ;
            char *p = buf + 14 ;

            //------------------------------------------------------------------
            // get "matrix" token and discard it
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            // printf ("header now [%s]\n", p) ;
            // printf ("compare %d\n", (strncmp (p, "matrix", 6))) ;

            if (strncmp (p, "matrix", 6) != 0)
            {
                // invalid Matrix Market object
                printf ("SPEX_mmread: bad object\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }
            p += 6 ;                                // skip past token "matrix"

            //------------------------------------------------------------------
            // get the fmt token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "coordinate", 10) == 0)
            {
                MM_fmt = MM_coordinate ;
                p += 10 ;
            }
            else if (strncmp (p, "array", 5) == 0)
            {
                MM_fmt = MM_array ;
                p += 5 ;
            }
            else
            {
                // invalid Matrix Market format
                printf ("SPEX_mmread: bad format\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }

            //------------------------------------------------------------------
            // get the Matrix Market type token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "real", 4) == 0)
            {
                MM_type = MM_real ;
                type = SPEX_FP64 ;
                p += 4 ;
            }
            else if (strncmp (p, "integer", 7) == 0)
            {
                MM_type = MM_integer ;
                type = SPEX_INT64 ;
                p += 7 ;
            }
            else if (strncmp (p, "pattern", 7) == 0)
            {
                MM_type = MM_pattern ;
                type = SPEX_INT64 ;
                p += 7 ;
            }
            else
            {
                // invalid Matrix Market type
                printf ("SPEX_mmread: bad type\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }

            //------------------------------------------------------------------
            // get the storage token
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any leading spaces

            if (strncmp (p, "general", 7) == 0)
            {
                MM_storage = MM_general ;
            }
            else if (strncmp (p, "symmetric", 9) == 0)
            {
                MM_storage = MM_symmetric ;
            }
            else if (strncmp (p, "skew-symmetric", 14) == 0)
            {
                MM_storage = MM_skew_symmetric ;
            }
            else
            {
                // invalid Matrix Market storage
                printf ("SPEX_mmread: bad type\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }

            //------------------------------------------------------------------
            // ensure the combinations are valid
            //------------------------------------------------------------------

            if (MM_type == MM_pattern)
            {
                // (coodinate) x (pattern) x (general or symmetric)
                if (! (MM_fmt == MM_coordinate &&
                        (MM_storage == MM_general ||
                         MM_storage == MM_symmetric)))
                {
                    // invalid combination
                    printf ("SPEX_mmread: bad pattern combo\n") ;
                    return (SPEX_INCORRECT_INPUT) ;
                }
            }

        }
        else if (is_blank_line (buf))
        {

            // -----------------------------------------------------------------
            // blank line or comment line
            // -----------------------------------------------------------------

            continue ;

        }
        else
        {

            // -----------------------------------------------------------------
            // read the first data line and return
            // -----------------------------------------------------------------

            // format: [nrows ncols nvals] or just [nrows ncols]

            int nitems = sscanf (buf, "%" SCNu64 " %" SCNu64 " %" SCNu64,
                &nrows, &ncols, &nvals) ;

            if (nitems == 2)
            {
                // a dense matrix
                if (!got_mm_header)
                {
                    // if no header, treat it as if it were
                    // %%MatrixMarket matrix array real general
                    MM_fmt = MM_array ;
                    MM_type = MM_real ;
                    MM_storage = MM_general ;
                    type = SPEX_FP64 ;
                }
                nvals = nrows * ncols ;
            }
            else if (nitems == 3)
            {
                // a sparse matrix
                if (!got_mm_header)
                {
                    // if no header, treat it as if it were
                    // %%MatrixMarket matrix coordinate real general
                    MM_fmt = MM_coordinate ;
                    MM_type = MM_real ;
                    MM_storage = MM_general ;
                    type = SPEX_FP64 ;
                }
            }
            else
            {
                // wrong number of items in first data line
                printf ("SPEX_mmread: bad 1st line\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }

            if (nrows != ncols)
            {
                if (! (MM_storage == MM_general))
                {
                    // a rectangular matrix must be in the general storage
                    printf ("SPEX_mmread: bad rectangular\n") ;
                    return (SPEX_INCORRECT_INPUT) ;
                }
            }

            //------------------------------------------------------------------
            // header has been read in
            //------------------------------------------------------------------

            break ;
        }
    }
    if (MM_storage == MM_symmetric || MM_storage == MM_skew_symmetric)
    {
        nvals *= 2;
    }

    //--------------------------------------------------------------------------
    // create the matrix
    //--------------------------------------------------------------------------

    SPEX_info info = SPEX_matrix_allocate (&A, SPEX_TRIPLET, type, nrows,
        ncols, nvals, false, true, option) ;
    if (info != SPEX_OK)
    {
        // failed to construct matrix
        // printf ("mmread: failed to construct A\n") ;
        return (info) ;
    }

    //--------------------------------------------------------------------------
    // quick return for empty matrix
    //--------------------------------------------------------------------------

    if (nrows == 0 || ncols == 0 || nvals == 0)
    {
        // success: return an empty matrix.  This is not an error.
        return (SPEX_OK) ;
    }

    //--------------------------------------------------------------------------
    // read the entries
    //--------------------------------------------------------------------------

    for (uint64_t k = 0 ; k < nvals ; k++)
    {

        //----------------------------------------------------------------------
        // get the next triplet, skipping blank lines and comment lines
        //----------------------------------------------------------------------

        uint64_t i, j ;
        char x [MAXLINE] ;

        for ( ; ; )
        {

            //------------------------------------------------------------------
            // read the file until finding the next triplet
            //------------------------------------------------------------------

            if (!get_line (f, buf))
            {
                // premature end of file - not enough triplets read in
                SPEX_matrix_free (&A, option) ;
                printf ("SPEX:mmread: premature EOF\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }
            if (is_blank_line (buf))
            {
                // blank line or comment
                continue ;
            }

            //------------------------------------------------------------------
            // get the row and column index
            //------------------------------------------------------------------

            char *p ;
            if (MM_fmt == MM_array)
            {
                // array format, column major order
                i = k % nrows ;
                j = k / nrows ;
                p = buf ;
                // printf ("array now [%s]\n", p) ;
            }
            else
            {
                // coordinate format; read the row and column index
                p = buf ;
                if (sscanf (p, "%" SCNu64 " %" SCNu64, &i, &j) != 2)
                {
                    // EOF or other I/O error
                    SPEX_matrix_free (&A, option) ;
                    printf ("SPEX_mmread: I/O error on indices\n") ;
                    return (SPEX_INCORRECT_INPUT) ;
                }
                // convert from 1-based to 0-based.
                i-- ;
                j-- ;
                // printf ("got (%g,%g)\n", (double) i, (double) j) ;
                // advance p to the 3rd token to get the value of the entry
                while (*p &&  isspace (*p)) p++ ;   // skip any leading spaces
                while (*p && !isspace (*p)) p++ ;   // skip nrows
                while (*p &&  isspace (*p)) p++ ;   // skip any spaces
                while (*p && !isspace (*p)) p++ ;   // skip nrows
                // printf ("now [%s]\n", p) ;
            }

            //------------------------------------------------------------------
            // read the value of the entry
            //------------------------------------------------------------------

            while (*p && isspace (*p)) p++ ;        // skip any spaces

            if (!read_entry (p, type, MM_type == MM_pattern, x))
            {
                // EOF or other I/O error, or value of entry out of range
                SPEX_matrix_free (&A, option) ;
                printf ("SPEX_mmread: I/O error on value\n") ;
                return (SPEX_INCORRECT_INPUT) ;
            }

            //------------------------------------------------------------------
            // set the value in the matrix
            //------------------------------------------------------------------

            info = set_value (A, type, i, j, x) ;
            if (info != SPEX_OK)
            {
                // unable to set element: invalid indices, or out of memory
                printf ("mmread: unable to set element\n") ;
                SPEX_matrix_free (&A, option) ;
                return (info) ;
            }

            // GxB_fprint (*A, GxB_COMPLETE, stdout) ;

            //------------------------------------------------------------------
            // also set the A(j,i) entry, if symmetric
            //------------------------------------------------------------------

            if (i != j && MM_storage != MM_general)
            {
                if (MM_storage == MM_symmetric)
                {
                    info = set_value (A, type, j, i, x) ;
                }
                else if (MM_storage == MM_skew_symmetric)
                {
                    negate_scalar (type, x) ;
                    info = set_value (A, type, j, i, x) ;
                }
                if (info != SPEX_OK)
                {
                    // unable to set element: invalid indices, or out of memory
                    SPEX_matrix_free (&A, option) ;
                    // printf ("mmread: unable to set symmetric element\n") ;
                    return (info) ;
                }
            }

            // one more entry has been read in
            break ;
        }
    }

    *A_handle = A;
    return (SPEX_OK) ;
}


// -------------------------------------------------------------------------
// read the matrices and construct LP obj
// -------------------------------------------------------------------------
#define MY_MAT_X(M, i) \
    (((M)->type == SPEX_FP64) ? (M)->x.fp64[i] : (M)->x.int64[i])

#define MY_SHIFT_X(M, i, j)       \
{                                 \
    if ((M)->type == SPEX_FP64)     \
    {                             \
        (M)->x.fp64[i] = (M)->x.fp64[j];\
    }                             \
    else                          \
    {                             \
        (M)->x.int64[i] = (M)->x.int64[j];\
    }                             \
}

#define MY_SHIFT_NEG_X(M, i, j)       \
{                                 \
    if ((M)->type == SPEX_FP64)     \
    {                             \
        (M)->x.fp64[i] = -1*((M)->x.fp64[j]);\
    }                             \
    else                          \
    {                             \
        (M)->x.int64[i] = -1*((M)->x.int64[j]);\
    }                             \
}

#define MY_FREE_WORK              \
{                                 \
    SPEX_matrix_free(&MA, option);\
    SPEX_matrix_free(&Mlb, option);\
    SPEX_matrix_free(&Mub, option);\
    SPEX_FREE(I);                 \
    SPEX_FREE(J);                 \
    SPEX_FREE(X);                 \
    SPEX_matrix_free(&MA_CSC, option);\
    SPEX_matrix_free(&Mb, option);\
    SPEX_matrix_free(&Mc, option);\
    if (File != NULL) fclose(File);\
}

#define MY_FREE_ALL               \
{                                 \
    MY_FREE_WORK;                 \
}

#define OK1(method)               \
{                                 \
    info = method;                \
    if (info != SPEX_OK)          \
    {                             \
        MY_FREE_ALL;              \
        return info;              \
    }                             \
}

SPEX_info SPEX_construct_LP
(
    glp_prob *LP,
    SPEX_matrix *A_handle,
    SPEX_matrix *b_handle,
    SPEX_matrix *c_handle,
    double *z0_handle,
    char *file_name,
    SPEX_options option
)
{
    SPEX_info info;
    double lb, ub;
    int file_name_len = strlen(file_name);
    SPEX_matrix MA = NULL, Mb = NULL, Mc = NULL, Mlb = NULL, Mub = NULL;
    SPEX_matrix MA_CSC = NULL, tmp = NULL;
    int *I = NULL, *J = NULL;
    double *X = NULL;
    FILE *File = NULL;
    char *suffix = "";

    // -------------------------------------------------------------------------
    // read A matrix
    // -------------------------------------------------------------------------
    suffix = ".mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }

    // Read matrix from given file as a triplet matrix
    OK1(SPEX_mmread(&MA, File, option));
    // convert to a CSC matrix
    OK1(SPEX_matrix_copy(&MA_CSC, SPEX_CSC, MA->type, MA, option));
    int64_t nz = MA->nz;
    SPEX_matrix_free(&MA, option);
    fclose(File); File = NULL;

    // -------------------------------------------------------------------------
    // read b matrix
    // -------------------------------------------------------------------------
    suffix = "_b.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }
    OK1(SPEX_mmread(&Mb, File, option));
    fclose(File); File = NULL;
    // convert to a dense matrix
    OK1(SPEX_matrix_copy(&tmp, SPEX_DENSE, Mb->type, Mb, option));
    SPEX_matrix_free(&Mb, option);
    Mb = tmp;

    // -------------------------------------------------------------------------
    // read c matrix
    // -------------------------------------------------------------------------
    suffix = "_c.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }
    OK1(SPEX_mmread(&Mc, File, option));
    fclose(File);File = NULL;
    // convert to a dense matrix
    OK1(SPEX_matrix_copy(&tmp, SPEX_DENSE, Mc->type, Mc, option));
    SPEX_matrix_free(&Mc, option);
    Mc = tmp;

    // -------------------------------------------------------------------------
    // read lb and ub matrix
    // -------------------------------------------------------------------------
    //read lower bound file
    suffix = "_lo.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }
    OK1(SPEX_mmread(&Mlb, File, option));
    fclose(File); File = NULL;
    // convert to a dense matrix
    OK1(SPEX_matrix_copy(&tmp, SPEX_DENSE, Mlb->type, Mlb, option));
    SPEX_matrix_free(&Mlb, option);
    Mlb = tmp;

    // read upper bound file
    suffix = "_hi.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }
    OK1(SPEX_mmread(&Mub, File, option));
    fclose(File);File = NULL;
    // convert to a dense matrix
    OK1(SPEX_matrix_copy(&tmp, SPEX_DENSE, Mub->type, Mub, option));
    SPEX_matrix_free(&Mub, option);
    Mub = tmp;

    // -------------------------------------------------------------------------
    // check dimension
    // -------------------------------------------------------------------------
    int64_t nvars = MA_CSC->n;
    int64_t neqs  = MA_CSC->m;
    if (Mlb->m != nvars || Mub->m != nvars || Mc->m != nvars || Mb->m != neqs)
    {
        printf("Dimension unmatched!\n");
        MY_FREE_ALL;
        return SPEX_INCORRECT_INPUT;
    }

    // -------------------------------------------------------------------------
    // find out if there are any fixed variables, and remove all fixed
    // variables if exist
    // -------------------------------------------------------------------------
    int64_t num_nz_removed = 0;
    int64_t num_var_removed = 0;
    for (int64_t j = 0; j < nvars; j++)
    {
        lb = (double) (MY_MAT_X(Mlb, j));
        ub = (double) (MY_MAT_X(Mub, j));
        //if (lb<0)printf("lb[%ld]=%f\n",j+1,lb);
        if (lb == ub)
        {
            //printf("x[%ld]=%f\n",j+1,lb);

            // remove all columns in A corresponding to these variables
            for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
            {
                int64_t i = MA_CSC->i[p];
                if (Mb->type == SPEX_FP64)
                {
                    Mb->x.fp64[i] -= (double) (MY_MAT_X(MA_CSC, p)) * lb;
                }
                else
                {
                    Mb->x.int64[i] -= (int64_t) (MY_MAT_X(MA_CSC, p) * lb);
                }
                num_nz_removed++;
            }
            (*z0_handle) += MY_MAT_X(Mc, j) * lb;
            num_var_removed++;
        }
        else if (num_var_removed > 0)
        {
            // shift all remaining entries
            int64_t shifted_j = j - num_var_removed;
            MY_SHIFT_X(Mc,  shifted_j, j);
            MY_SHIFT_X(Mlb, shifted_j, j);
            MY_SHIFT_X(Mub, shifted_j, j);
            MA_CSC->p[shifted_j] = MA_CSC->p[j] - num_nz_removed;
            for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
            {
                int64_t shifted_p = p - num_nz_removed;
                MA_CSC->i[shifted_p] = MA_CSC->i[p];
                MY_SHIFT_X(MA_CSC, shifted_p, p);
            }
        }
    }
    int64_t remain_vars = nvars - num_var_removed;
    int64_t remain_nz   = nz - num_nz_removed;
    Mc->nz  = remain_vars;
    Mc->m   = remain_vars;
    Mlb->nz = remain_vars;
    Mlb->m  = remain_vars;
    Mub->nz = remain_vars;
    Mub->m  = remain_vars;
    MA_CSC->p[remain_vars] = remain_nz;
    MA_CSC->n  = remain_vars;
    printf("m=%ld n=%ld fixed_var=%ld nz=%ld real_nnz=%ld\n", neqs,
        nvars, num_var_removed, nz, remain_nz);

    // -------------------------------------------------------------------------
    // check bounds for each variables to make sure lb[j] = 0
    // -------------------------------------------------------------------------
    // find out the number of new columns need by check if there is any bounds
    // for xi is (-inf,inf)
    int64_t new_vars = 0;
    int64_t new_nz = 0;
    for (int64_t j = 0; j < remain_vars; j++)
    {
        lb = (double) (MY_MAT_X(Mlb, j));
        ub = (double) (MY_MAT_X(Mub, j));
        if (lb == -1e308 && ub == 1e308)
        {
            new_nz += MA_CSC->p[j+1] - MA_CSC->p[j];
            new_vars++;
        }
    }
    if (new_vars > 0)
    {
        if (MA_CSC->nzmax < remain_nz+new_nz)
        {
            bool okx, oki;
            if (MA_CSC->type == SPEX_FP64)
            {
                MA_CSC->x.fp64 = (double*)SPEX_realloc(remain_nz+new_nz,
                    MA_CSC->nzmax, sizeof(double), MA_CSC->x.fp64, &okx);
            }
            else
            {
                MA_CSC->x.int64 = (int64_t*)SPEX_realloc(remain_nz+new_nz,
                    MA_CSC->nzmax, sizeof(int64_t), MA_CSC->x.int64, &okx);
            }
            MA_CSC->i = (int64_t*)SPEX_realloc(remain_nz+new_nz,
                MA_CSC->nzmax, sizeof(int64_t), MA_CSC->i, &oki);
            if (!okx || !oki)
            {
                MY_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }
        }
        MA_CSC->nzmax = remain_nz;
        if (nvars < remain_vars+new_vars)
        {
            bool okp, okc, okl, oku;
            MA_CSC->p = (int64_t*)SPEX_realloc(remain_vars+new_vars+1,
                nvars+1, sizeof(int64_t), MA_CSC->p, &okp);
            //expand Mc
            if (Mc->type == SPEX_FP64)
            {
                Mc->x.fp64 = (double*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(double), Mc->x.fp64, &okc);
            }
            else
            {
                Mc->x.int64 = (int64_t*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(int64_t), Mc->x.int64, &okc);
            }

            //expand Mlb
            if (Mlb->type == SPEX_FP64)
            {
                Mlb->x.fp64 = (double*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(double), Mlb->x.fp64, &okl);
            }
            else
            {
                Mlb->x.int64 = (int64_t*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(int64_t), Mlb->x.int64, &okl);
            }
            //expand Mub
            if (Mub->type == SPEX_FP64)
            {
                Mub->x.fp64 = (double*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(double), Mub->x.fp64, &oku);
            }
            else
            {
                Mub->x.int64 = (int64_t*)SPEX_realloc(remain_vars+new_vars,
                    nvars, sizeof(int64_t), Mub->x.int64, &oku);
            }

            if (!okp || !okc || !okl || !oku)
            {
                MY_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }
        }
        Mc->nzmax  = remain_vars;
        Mlb->nzmax = remain_vars;
        Mub->nzmax = remain_vars;

        // adding entries
        for (int64_t j = 0; j < remain_vars; j++)
        {
            lb = (double) (MY_MAT_X(Mlb, j));
            ub = (double) (MY_MAT_X(Mub, j));
            if (lb == -1e308)
            {
                if (ub == 1e308)
                {
                    // let xj = xj1 - xj2, where xj1 > 0 and xj2 > 0
                    for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
                    {
                        MA_CSC->i[remain_nz] = MA_CSC->i[p];
                        MY_SHIFT_NEG_X(MA_CSC, remain_nz, p);
                        remain_nz++;
                    }
                    MA_CSC->p[remain_vars] = remain_nz;
                    if (Mlb->type == SPEX_FP64)
                    {
                        Mlb->x.fp64[j] = 0;
                    }
                    else
                    {
                        Mlb->x.int64[j] = 0;
                    }
                    MY_SHIFT_NEG_X(Mc,  remain_vars, j);
                    MY_SHIFT_X(Mlb, remain_vars, j);
                    MY_SHIFT_X(Mub, remain_vars, j);
                    remain_vars++;
                }
                else
                {
                    // let xj = -(xj-ub[j]) so that xj > 0
                    for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
                    {
                        int64_t i = MA_CSC->i[p];
                        if (Mb->type == SPEX_FP64)
                        {
                            Mb->x.fp64[i] -= (double) (MY_MAT_X(MA_CSC, p)) *ub;
                        }
                        else
                        {
                            Mb->x.int64[i] -= (int64_t)(MY_MAT_X(MA_CSC, p)*ub);
                        }
                        if (MA_CSC->type == SPEX_FP64)
                        {
                            MA_CSC->x.fp64[p] *= -1;
                        }
                        else
                        {
                            MA_CSC->x.int64[p] *= -1;
                        }
                    }
                    (*z0_handle) += MY_MAT_X(Mc, j) * ub;
                    if (Mc->type == SPEX_FP64)
                    {
                        Mc->x.fp64[j] *= -1;
                    }
                    else
                    {
                        Mc->x.int64[j] *= -1;
                    }
                    //ub = - lb;
                    if (Mub->type == SPEX_FP64)
                    {
                        Mub->x.fp64[j] = -lb;
                    }
                    else
                    {
                        Mub->x.int64[j] = -1*(MY_MAT_X(Mlb, j));
                    }
                    //lb = 0.0;
                    if (Mlb->type == SPEX_FP64)
                    {
                        Mlb->x.fp64[j] = 0;
                    }
                    else
                    {
                        Mlb->x.int64[j] = 0;
                    }
                }
            }
            //------------------------------------------------------------------
            // find if there is any nonzero lower bound, shift the corresponding
            // variables such that the lower bound is 0.
            //------------------------------------------------------------------
         /*   else if (lb != 0)
            {
                //printf("%ld: %lf %lf\n", j, lb, ub);
                for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
                {
                    int64_t i = MA_CSC->i[p];
                    if (Mb->type == SPEX_FP64)
                    {
                        Mb->x.fp64[i] -= (double) (MY_MAT_X(MA_CSC, p)) * lb;
                    }
                    else
                    {
                        Mb->x.int64[i] -= (int64_t) (MY_MAT_X(MA_CSC, p) * lb);
                    }
                }
                ub = ub - lb;
                (*z0_handle) += MY_MAT_X(Mc, j) * lb;
                lb = 0.0;
                //printf("%ld: %lf %lf\n\n\n", j, lb, ub);
            }*/
        }
        Mc->nz  = remain_vars;
        Mc->m   = remain_vars;
        Mlb->nz = remain_vars;
        Mlb->m  = remain_vars;
        Mub->nz = remain_vars;
        Mub->m  = remain_vars;
        MA_CSC->p[remain_vars] = remain_nz;
        MA_CSC->n  = remain_vars;
        printf("m=%ld n=%ld new_vars=%ld real_nnz=%ld\n", neqs, remain_vars,
            new_vars, remain_nz);
    }

    // -------------------------------------------------------------------------
    // load matrices to LP
    // -------------------------------------------------------------------------
    I = (int*)    SPEX_malloc((remain_nz+1) * sizeof(int));
    J = (int*)    SPEX_malloc((remain_nz+1) * sizeof(int));
    X = (double*) SPEX_malloc((remain_nz+1) * sizeof(double));
    if (!I || !J || !X)
    {
        MY_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    I[0] = 0;
    J[0] = 0;
    X[0] = 0;

    for (int j = 0; j < remain_vars; j++)
    {
        for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
        {
            J[p+1] = j + 1;
            I[p+1] = (int) (MA_CSC->i[p]) + 1;
            X[p+1] = (double) (MY_MAT_X(MA_CSC, p));
        }
    }

    glp_add_rows(LP, neqs);
    glp_add_cols(LP, remain_vars);
    glp_load_matrix(LP, (int) remain_nz, I, J, X);

    for (int64_t j = 0; j < remain_vars; j++)
    {
        lb = (double) (MY_MAT_X(Mlb, j));
        ub = (double) (MY_MAT_X(Mub, j));
        if (lb < 0)
        {
            printf("%ld: lb=%lf(<0) ub=%lf coef=%lf\n", j, lb, ub,
                MY_MAT_X(Mc, j));
        }
        // find if there is any positive lower bound, shift the corresponding
        // variables such that the lower bound is 0.
        if (lb != 0)
        {
            //printf("%ld: %lf %lf\n", j, lb, ub);
            for (int64_t p = MA_CSC->p[j]; p < MA_CSC->p[j+1]; p++)
            {
                int64_t i = MA_CSC->i[p];
                if (Mb->type == SPEX_FP64)
                {
                    Mb->x.fp64[i] -= (double) (MY_MAT_X(MA_CSC, p)) * lb;
                }
                else
                {
                    Mb->x.int64[i] -= (int64_t) (MY_MAT_X(MA_CSC, p) * lb);
                }
            }
            ub = ub - lb;
            (*z0_handle) += MY_MAT_X(Mc, j) * lb;
            lb = 0.0;
            //printf("%ld: %lf %lf\n\n\n", j, lb, ub);
        }
        if (ub <= 0)
        {
            printf("%ld: %lf %lf\n", j, lb, ub);
        }
        glp_set_col_bnds(LP, j+1, GLP_DB, lb, ub);
        /*if (MY_MAT_X(Mc, j) < 0)
        {
            printf("c[%ld]=%lf<0\n",j,MY_MAT_X(Mc,j));
        }*/
        glp_set_obj_coef(LP, j+1, MY_MAT_X(Mc, j));
    }

    for (int64_t i = 0; i < neqs; i++)
    {
        glp_set_row_bnds(LP, i+1, GLP_FX, MY_MAT_X(Mb, i), 0.0);
    }

    // convert to mpz matrices
    OK1(SPEX_matrix_copy(A_handle, SPEX_CSC,   SPEX_MPZ, MA_CSC, option));
    OK1(SPEX_matrix_copy(c_handle, SPEX_DENSE, SPEX_MPZ, Mc, option));
    // convert to double matrix if needed
    if (Mb->type != SPEX_FP64)
    {
        OK1(SPEX_matrix_copy(b_handle, SPEX_DENSE, SPEX_FP64, Mb, option));
    }
    else
    {
        *b_handle = Mb;
        Mb = NULL;
    }

    // free the workspace
    MY_FREE_WORK;
    return SPEX_OK;
}


// function to compute A = A+v*v^T

SPEX_info SPEX_A_plus_vvT
(
    SPEX_matrix A0, // m-by-n SPEX_DYNAMIC_CSC matrix as A
    const SPEX_matrix M, // m-by-k SPEX_CSC matrix, whose j-th column is v
    const int64_t j // j-th column of M is used as v
)
{
    int64_t p, pj, pi, j1, i1, target, target_sym;
    SPEX_info info;

    for (pj = M->p[j]; pj < M->p[j+1]; pj++)
    {
        j1 = M->i[pj];

        target = -1;
        for (p = 0; p < A0->v[j1]->nz; p++)
        {
            if (A0->v[j1]->i[p] == j1) {target = p; break;}
        }
        if (target == -1)
        {
            target = A0->v[j1]->nz;
            A0->v[j1]->i[target] = j1;
            A0->v[j1]->nz++;
        }
        OK(SPEX_mpz_addmul(A0->v[j1]->x[target], M->x.mpz[pj], M->x.mpz[pj]));

        for (pi = pj+1; pi < M->p[j+1]; pi++)
        {
            i1 = M->i[pi];

            target = -1;
            for (p = 0; p < A0->v[j1]->nz; p++)
            {
                if (A0->v[j1]->i[p] == i1) {target = p; break;}
            }
            if (target == -1)
            {
                target = A0->v[j1]->nz;
                A0->v[j1]->i[target] = i1;
                A0->v[j1]->nz++;
                target_sym = A0->v[i1]->nz;
                A0->v[i1]->i[target_sym] = j1;
                A0->v[i1]->nz++;
            }
            else
            {
                target_sym = -1;
                for (p = 0; p < A0->v[i1]->nz; p++)
                {
                    if (A0->v[i1]->i[p] == j1) {target_sym = p; break;}
               }
                if (target_sym == -1)
                {
                    printf("this shouldn't happen, ");
                    printf("have fun finding the bug\n");
                    return 0;
                }
            }
            info = SPEX_mpz_addmul(A0->v[j1]->x[target],
                               M->x.mpz[pj], M->x.mpz[pi]);
            if (info != SPEX_OK) return info;
            info = SPEX_mpz_set(A0->v[i1]->x[target_sym],
                            A0->v[j1]->x[target]);
            if (info != SPEX_OK) return info;
        }
    }
    return SPEX_OK;
}

SPEX_info SPEX_matrix_equal
(
    bool *Isequal,
    const SPEX_matrix L1,
    const SPEX_matrix L_update,
    const int64_t *P_update
)
{
    *Isequal = false;
    int64_t i, j, p, pi, n = L1->m;
    int sgn;
    SPEX_info info;
    mpz_t tmpz, tmpz1;
    info = SPEX_mpz_init(tmpz);
    if (info != SPEX_OK) return info;
    info = SPEX_mpz_init(tmpz1);
    if (info != SPEX_OK)
    {
        mpz_clear(tmpz);
        return info;
    }

    for (j = 0; j < n; j++)
    {
        for (p = L1->p[j]; p < L1->p[j+1]; p++)
        {
            i = L1->i[p];
            //gmp_printf("%ld %Zd\n",i, L1->x.mpz[p]);
            info = SPEX_mpz_sgn(&sgn, L1->x.mpz[p]);
            if (info != SPEX_OK)
            {
                mpz_clear(tmpz);
                mpz_clear(tmpz1);
                return info;
            }
            if (sgn == 0) {continue;}

            bool found = false;
            for (pi = 0; pi < L_update->v[j]->nz; pi++)
            {
                //printf("P[%ld]=%ld %ld\n",i, P_update[i],
                //    L_update->v[j]->i[pi]);
                if (L_update->v[j]->i[pi] == P_update[i])
                {
                    mpq_get_num(tmpz1, L_update->v[j]->scale);
                    mpz_mul(tmpz, L_update->v[j]->x[pi], tmpz1);
                    mpq_get_den(tmpz1, L_update->v[j]->scale);
                    mpz_divexact(tmpz, tmpz, tmpz1);
                    info = SPEX_mpz_cmp(&sgn, tmpz, L1->x.mpz[p]);
                    if (info != SPEX_OK)
                    {
                        mpz_clear(tmpz);
                        mpz_clear(tmpz1);
                        return info;
                    }
                    if (sgn != 0)
                    {
                        printf("found different value for L(%ld, %ld)\n",
                            i, j);
                        gmp_printf("factor: %Zd\nupdate: %Zd\n",
                            L1->x.mpz[p], tmpz);
                        gmp_printf("update: %Zd\nscale:%Qd\n",
                            L_update->v[j]->x[pi],
                            L_update->v[j]->scale);
                        *Isequal = false;
                        mpz_clear(tmpz);
                        mpz_clear(tmpz1);
                        return SPEX_OK;
                    }
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                printf("failed to find corresponding L(%ld, %ld)\n", i, j);
                *Isequal = false;
                mpz_clear(tmpz);
                mpz_clear(tmpz1);
                return SPEX_OK;
            }
        }
    }
    *Isequal = true;
    mpz_clear(tmpz);
    mpz_clear(tmpz1);
    return SPEX_OK;
}

//------------------------------------------------------------------------------
// spex_update_verify.c: verify if A=LD^(-1)U after factorization update
//------------------------------------------------------------------------------

/* Purpose: This function is to verify if A=L(P,:)D^(-1)U(:,Q) after
 * factorization update. This is done by solving LD^(-1)U*x=b via the updated
 * factorization, and check if A*x=b holds rational-arthmetically. This
 * function is provided here only used for debugging purposes, as the routines
 * within SPEX are gauranteed to be exact.
 */


#undef MY_FREE_ALL
#define MY_FREE_ALL                   \
{                                     \
    SPEX_matrix_free(&b, option);     \
    SPEX_matrix_free(&x, option);     \
    SPEX_matrix_free(&b2, option);    \
    mpq_clear(temp);                  \
}


SPEX_info MY_update_verify
(
    bool *Is_correct,     // if the factorization is correct
    SPEX_factorization F,// LU factorization of A
    const SPEX_matrix A,     // Input matrix of CSC MPZ
    const SPEX_options option// command options
)
{
    SPEX_info info;
    int64_t tmp, i, n = F->L->n;
    int r;
    mpq_t temp;
    SPEX_matrix b = NULL; // the dense right-hand-side matrix to be generated
    SPEX_matrix x = NULL; // the dense solution matrix to be generated
    SPEX_matrix b2 = NULL; // the dense matrix to store the result of A*x

    OK1(SPEX_mpq_init(temp));
    OK1(SPEX_matrix_allocate(&b , SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));
    OK1(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPQ, n, 1, n, false,
        true, option));

    // -------------------------------------------------------------------------
    // generate random right-hand-size vector
    // -------------------------------------------------------------------------
    // initialize random number generator
    int seed = 10;
    srand(seed);
    for (i = 0; i < n; i++)
    {
        tmp = i+1;//rand(); //TODO
        OK1(SPEX_mpz_set_si(b->x.mpz[i], tmp));
    }

    // -------------------------------------------------------------------------
    // solve LD^(-1)Ux = b for x
    // -------------------------------------------------------------------------
    OK1(SPEX_update_solve(&x, F, b, option));

    // -------------------------------------------------------------------------
    // compute b2 = A*x
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        OK1(SPEX_mpq_sgn(&r, x->x.mpq[i]));
        if (r == 0) { continue;}

        for (int64_t p = A->p[i]; p < A->p[i+1]; p++)
        {
            int64_t j = A->i[p];
            OK1(SPEX_mpq_set_z(temp, A->x.mpz[p]));
            // b2[j] += x[i]*A(j,i)
            OK1(SPEX_mpq_mul(temp, temp, x->x.mpq[i]));
            OK1(SPEX_mpq_add(b2->x.mpq[j], b2->x.mpq[j], temp));
        }
    }
    //--------------------------------------------------------------------------
    // Apply scales of A and b to b2 before comparing the b2 with scaled b'
    //--------------------------------------------------------------------------
    OK1(SPEX_mpq_div(temp, b->scale, A->scale));

    // Apply scaling factor, but ONLY if it is not 1
    OK1(SPEX_mpq_cmp_ui(&r, temp, 1, 1));
    if (r != 0)
    {
        for (i = 0; i < n; i++)
        {
            OK1(SPEX_mpq_mul(b2->x.mpq[i], b2->x.mpq[i], temp));
        }
    }

    // -------------------------------------------------------------------------
    // check if b2 == b
    // -------------------------------------------------------------------------
    *Is_correct = true;
    for (i = 0; i < n; i++)
    {
        // temp = b[i] (correct b)
        OK1(SPEX_mpq_set_z(temp, b->x.mpz[i]));

        // set check false if b!=b2
        OK1(SPEX_mpq_equal(&r, temp, b2->x.mpq[i]));
        if (r == 0)
        {
            *Is_correct = false;
            break;
        }
    }

    //--------------------------------------------------------------------------
    // Print result
    //--------------------------------------------------------------------------

    if (!(*Is_correct))
    {
        // This can never happen.
        printf ("ERROR! Factorization is wrong.\n") ;
    }
    else
    {
        printf ("Factorization is verified to be correct and exact.\n") ;
    }

    MY_FREE_ALL;
    return SPEX_OK;
}
