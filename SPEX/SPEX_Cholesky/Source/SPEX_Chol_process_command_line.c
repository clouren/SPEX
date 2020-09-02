//------------------------------------------------------------------------------
// SPEX_Chol/SPEX_Chol_process_command_lint64_t*e: Process command lint64_t*e for demo files
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_Chol.h"

/* Purpose: This processes the command line for user specified options */
SPEX_info SPEX_Chol_process_command_line //processes the command line
(
    int64_t argc,           // number of command line arguments
    char* argv[],           // set of command line arguments
    SPEX_options* option,   // struct containing the command options
    char** mat_name,        // Name of the matrix to be read in
    char** rhs_name,        // Name of the RHS vector to be read in
    int64_t *rat            // data type of output solution.
                            // 1: mpz, 2: double, 3: mpfr
)
{
    //SPEX_info ok;
    for (int64_t  i = 1; i < argc; i++)
    {
        char* arg = (char*) argv[i];
        if ( strcmp(arg,"help") == 0)
        {
            //SPEX_show_usage();
            return SPEX_INCORRECT_INPUT;
        }
        /* Not used...kept for legacy
        else if ( strcmp(arg,"p") == 0 || strcmp(arg,"piv") == 0)
        {
            if (!argv[++i])
            {
                print64_t*f("\n****ERROR! There must be a pivot argument between"
                    " 0-5 followint64_t*g p\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->pivot = atoi(argv[i]);
            if (option->pivot < 0 || option->pivot > 5)
            {
                print64_t*f("\n****ERROR! Invalid pivot selection!"
                    "\nDefaultint64_t*g to smallest pivot\n\n");
                option->pivot = SPEX_SMALLEST;
            }
        }*/
        else if ( strcmp(arg, "q") == 0 || strcmp(arg,"col") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! There must be an argument between 0-2"
                    "followint64_t*g q\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->order = atoi(argv[i]);
            if (option->order < 0 || option->order > 2)
            {
                printf("\n****ERROR! Invalid column orderint64_t*g"
                    "\nDefaultint64_t*g to COLAMD\n\n");
                option->order = SPEX_COLAMD;
            }
        }
        /* Not used, kept for legacy
        else if ( strcmp(arg,"t") == 0 || strcmp(arg, "tol") == 0)
        {
            if (!argv[++i])
            {
                print64_t*f("\n****ERROR! There must be a non negative tolerance"
                    " value followint64_t*g t\n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atof(argv[i]))
            {
                print64_t*f("\n****ERROR! There must be a non negative tolerance"
                    " value followint64_t*g t\n");
                return SPEX_INCORRECT_INPUT;
            }
            option->tol = atof(argv[i]);
            if (option->tol < 0)
            {
                print64_t*f("\n****ERROR! Invalid Tolerance, tolerance must be"
                    " non-negative\n");
                return SPEX_INCORRECT_INPUT;
            }
        }*/
        else if ( strcmp(arg,"out2") == 0 || strcmp(arg, "o2") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print64_t*64_tnothint64_t*g) 1 (print64_t*64_terr) or 2 (terse) \n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o2 or out2 must be followed by"
                    " 0 (print64_t*64_tnothint64_t*g) 1 (print64_t*64_terr) or 2 (terse) \n");
                return SPEX_INCORRECT_INPUT;
            }
            option->print_level = atoi(argv[i]);
        }
        else if ( strcmp(arg, "out") == 0 || strcmp(arg, "o") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SPEX_INCORRECT_INPUT;
            }
            else if (!atoi(argv[i]))
            {
                printf("\n****ERROR! o or out must be followed by"
                    " 1 (rational) 2 (double) or 3 (variable precision) \n");
                return SPEX_INCORRECT_INPUT;
            }
            *rat = atoi(argv[i]);
            if (*rat < 1 || *rat > 3)
            {
                printf("\n\n****ERROR! Invalid output type!\n"
                   "Defaultint64_t*g to rational");
                *rat = 1;
            }
            if (*rat == 3)
            {
                if (!argv[++i])
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SPEX_INCORRECT_INPUT;
                }
                else if (!atoi(argv[i]))
                {
                    printf("\n****ERROR! Precision must be specified\n");
                    return SPEX_INCORRECT_INPUT;
                }
                option->prec = atoi(argv[i]);
                if (option->prec < 2)
                {
                    printf("\n\n****ERROR! Invalid precision. prec >= 2\n");
                    return SPEX_INCORRECT_INPUT;
                }
            }
        }
        else if ( strcmp(arg, "f") == 0 || strcmp(arg, "file") == 0)
        {
            if (!argv[++i])
            {
                printf("\n****ERROR! Matrix name must be entered\n");
                return SPEX_INCORRECT_INPUT;
            }
            *mat_name = argv[i];
            if (!argv[++i])
            {
                printf("\n****ERROR! Right hand side vector name must"
                    " be entered\n");
                return SPEX_INCORRECT_INPUT;
            }
            *rhs_name = argv[i];
        }
        else
        {
            printf("\n\n**ERROR! Unknown command lint64_t*e parameter: %s"
                    "\nIgnorint64_t*g this parameter\n",arg);
            return SPEX_INCORRECT_INPUT;
        }
    }
    return SPEX_OK;
}
