/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    ulong cflags = UWORD(0);

    FLINT_TEST_INIT(state);

    flint_printf("cosh_series....");
    fflush(stdout);    

    /* Check aliasing of a and c */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a, b;
        slong n = n_randint(state, 50) + 1;

        fmpq_poly_init(a);
        fmpq_poly_init(b);

        fmpq_poly_randtest_not_zero(a, state, n_randint(state, 50) + 1, 50);
        fmpq_poly_set_coeff_ui(a, 0, UWORD(0));

        fmpq_poly_canonicalise(a);

        fmpq_poly_cosh_series(b, a, n);
        fmpq_poly_cosh_series(a, a, n);

        cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(b) ? 0 : 2;
        result = (fmpq_poly_equal(a, b) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpq_poly_debug(a), flint_printf("\n\n");
            fmpq_poly_debug(b), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(a);
        fmpq_poly_clear(b);
    }

    /* Check cosh(A)^2-1 = sinh(A)^2 */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t A, coshA, sinhA, B, C, one;
        slong n = n_randint(state, 80) + 1;

        fmpq_poly_init(A);
        fmpq_poly_init(coshA);
        fmpq_poly_init(sinhA);
        fmpq_poly_init(B);
        fmpq_poly_init(C);
        fmpq_poly_init(one);

        fmpq_poly_randtest_not_zero(A, state, n_randint(state, 60) + 1, 80);
        fmpq_poly_set_coeff_ui(A, 0, UWORD(0));

        fmpq_poly_cosh_series(coshA, A, n);
        fmpq_poly_sinh_series(sinhA, A, n);
        fmpq_poly_mullow(B, coshA, coshA, n);
        fmpq_poly_set_coeff_ui(one, 0, UWORD(1));
        fmpq_poly_sub(B, B, one);
        fmpq_poly_mullow(C, sinhA, sinhA, n);

        cflags |= fmpq_poly_is_canonical(coshA) ? 0 : 1;
        cflags |= fmpq_poly_is_canonical(sinhA) ? 0 : 2;
        result = (fmpq_poly_equal(B, C) && !cflags);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("A = "), fmpq_poly_debug(A), flint_printf("\n\n");
            flint_printf("cosh(A) = "), fmpq_poly_debug(coshA), flint_printf("\n\n");
            flint_printf("sinh(A) = "), fmpq_poly_debug(sinhA), flint_printf("\n\n");
            flint_printf("cflags = %wu\n\n", cflags);
            fflush(stdout);
            flint_abort();
        }

        fmpq_poly_clear(A);
        fmpq_poly_clear(coshA);
        fmpq_poly_clear(sinhA);
        fmpq_poly_clear(B);
        fmpq_poly_clear(C);
        fmpq_poly_clear(one);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
