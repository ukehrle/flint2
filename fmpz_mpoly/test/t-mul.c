/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"


int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_add(t1, g, h, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul(k1, f, t1, ctx);
            fmpz_mpoly_assert_canonical(k1, ctx);
            fmpz_mpoly_mul(t1, f, g, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul(t2, f, h, ctx);
            fmpz_mpoly_assert_canonical(t2, ctx);
            fmpz_mpoly_add(k2, t1, t2, ctx);
            fmpz_mpoly_assert_canonical(k2, ctx);
            result = fmpz_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_mul(f, f, g, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_mul(g, f, g, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            result = fmpz_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

