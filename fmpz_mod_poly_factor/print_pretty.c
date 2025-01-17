/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "fmpz_mod_poly.h"

void fmpz_mod_poly_factor_print_pretty(const fmpz_mod_poly_factor_t fac,
                                    const char *var, const fmpz_mod_ctx_t ctx)
{
    slong i;
    for (i = 0; i < fac->num; i++)
    {
        fmpz_mod_poly_print_pretty(fac->poly + i, var, ctx);
        flint_printf(" ^ %wd\n", fac->exp[i]);
    }
}
