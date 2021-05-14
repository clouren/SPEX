//------------------------------------------------------------------------------
// SPEX_Update/Test/test.c: support functions for the test programs
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

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
    SPEX_matrix *A,
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
    SPEX_matrix **A_handle,// handle of matrix to create
    FILE *f,             // file to read from, already open
    SPEX_options *option
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
    SPEX_matrix *A = NULL;

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


// load constraint matrix to the LP problem
#define MY_MATRIX_ENTRY(M, i) \
    (M->type == SPEX_FP64) ? M->x.fp64[i] : (double) M->x.int64[i]
SPEX_info SPEX_load_matrix_to_LP
(
    SPEX_matrix *A, 
    glp_prob *LP
)
{
    if (A->kind != SPEX_TRIPLET)
    {
        printf("input matrix is not triplet matrix\n");
        return SPEX_INCORRECT_INPUT;
    }
    printf("m=%ld n=%ld nz=%ld\n",A->m, A->n,A->nz);
    glp_add_rows(LP, A->m);
    glp_add_cols(LP, A->n);
    int64_t nz = A->nz;
    int *i = NULL, *j = NULL;
    double *x = NULL;
    i = (int*)    SPEX_malloc((nz+1) * sizeof(int));
    j = (int*)    SPEX_malloc((nz+1) * sizeof(int));
    x = (double*) SPEX_malloc((nz+1) * sizeof(double));
    if (!i || !j || !x)
    {
        SPEX_FREE(i);
        SPEX_FREE(j);
        SPEX_FREE(x);
        return SPEX_OUT_OF_MEMORY;
    }
    i[0] = 0;
    j[0] = 0;
    x[0] = 0;
    for (int64_t p = 0; p < nz; p++)
    {
        i[p+1] = (int) A->i[p] + 1;
        j[p+1] = (int) A->j[p] + 1;
        x[p+1] = MY_MATRIX_ENTRY(A, p);
    }
    glp_load_matrix(LP, (int) nz, i, j, x);
    SPEX_FREE(i);
    SPEX_FREE(j);
    SPEX_FREE(x);
    return SPEX_OK;
}

// read the matrices and construct LP obj
#define OK1(method)                \
{                                 \
    info = method;                \
    if (info != SPEX_OK)          \
    {                             \
        SPEX_matrix_free(&A, option);\
        SPEX_matrix_free(&b, option);\
        SPEX_matrix_free(&c, option);\
        SPEX_matrix_free(&M1, option);\
        SPEX_matrix_free(&M2, option);\
        if (File != NULL) fclose(File);\
        return;                    \
    }                              \
}

void SPEX_construct_LP
(
    glp_prob *LP,
    SPEX_matrix **A_handle,
    SPEX_matrix **b_handle,
    SPEX_matrix **c_handle,
    char *file_name,
    SPEX_options *option
)
{
    SPEX_info info;
    double lb, ub;
    int file_name_len = strlen(file_name);
    SPEX_matrix *A = NULL, *b = NULL, *c = NULL, *M1 = NULL, *M2 = NULL;
    FILE *File = NULL;
    char *suffix = "";

    // ---------------------------------------------------------------------------
    // read A matrix
    // ---------------------------------------------------------------------------
    suffix = ".mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file");
        return;
    }

    // Read matrix from given file as a triplet matrix
    OK1(SPEX_mmread(&M1, File, option));
    fclose(File); File = NULL;
    // load matrix to LP
    OK1(SPEX_load_matrix_to_LP(M1, LP));
    // convert matrix to a CSC matrix
    OK1(SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, M1, option));
    // free the triplet matrix
    OK1(SPEX_matrix_free(&M1, option));

    // ---------------------------------------------------------------------------
    // read b matrix
    // ---------------------------------------------------------------------------
    suffix = "_b.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        return;
    }
    OK1(SPEX_mmread(&M1, File, option));
    fclose(File); File = NULL;
    for (int64_t i = 0; i < M1->m; i++)
    {
        glp_set_row_bnds(LP, i+1, GLP_FX, MY_MATRIX_ENTRY(M1, i), 0.0);
    }
    // convert matrix to a CSC matrix
    OK1(SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, M1, option));
    // free the triplet matrix
    SPEX_matrix_free(&M1, option);

    // ---------------------------------------------------------------------------
    // read c matrix
    // ---------------------------------------------------------------------------
    suffix = "_c.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        return;
    }
    OK1(SPEX_mmread(&M1, File, option));
    fclose(File);File = NULL;
    for (int64_t i = 0; i < M1->m; i++)
    {
        glp_set_obj_coef(LP, i+1, MY_MATRIX_ENTRY(M1, i));
    }
    // convert matrix to a CSC matrix
    OK1(SPEX_matrix_copy(&c, SPEX_DENSE, SPEX_MPZ, M1, option));
    // free the triplet matrix
    SPEX_matrix_free(&M1, option);

    // ---------------------------------------------------------------------------
    // read lb and ub matrix
    // ---------------------------------------------------------------------------
    suffix = "_lo.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        return;
    }
    info = SPEX_mmread(&M1, File, option);
    fclose(File);
    if (info != SPEX_OK)
    {
        SPEX_matrix_free(&M1, option);
        SPEX_matrix_free(&M2, option);
        return;
    }
    suffix = "_hi.mtx";
    strcpy(file_name+file_name_len,suffix);
    printf("reading file %s\n",file_name);
    File = fopen(file_name, "r");
    if (File == NULL)
    {
        perror("Error while opening file\n");
        return;
    }
    OK1(SPEX_mmread(&M2, File, option));
    fclose(File);File = NULL;
    for (int64_t i = 0; i < M1->m; i++)
    {
        lb = MY_MATRIX_ENTRY(M1, i);
        ub = MY_MATRIX_ENTRY(M2, i);
        if (lb<0)printf("lb[%ld]=%f\n",i+1,lb);
        if (lb == ub)
        {
            glp_set_col_bnds(LP, i+1, GLP_FX, lb, 0.0);
            printf("x[%ld]=%f\n",i+1,lb);
        }
        else
        {
            glp_set_col_bnds(LP, i+1, GLP_DB, lb, ub);
        }
    }
    SPEX_matrix_free(&M1, option);
    SPEX_matrix_free(&M2, option);
    *A_handle = A;
    *b_handle = b;
    *c_handle = c;
    return;
}