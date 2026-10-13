/* ========================================================================== */
/* === cholmod_extra ======================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * Copyright (C) 2012 Jaakko Luttinen
 *
 * cholmod_extra.h is licensed under Version 2 of the GNU General
 * Public License, or (at your option) any later version. See LICENSE
 * for a text of the license.
 * -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
 * This file is part of CHOLMOD Extra Module.
 *
 * CHOLDMOD Extra Module is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of
 * the License, or (at your option) any later version.
 *
 * CHOLMOD Extra Module is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHOLMOD Extra Module.  If not, see
 * <http://www.gnu.org/licenses/>.
 * -------------------------------------------------------------------------- */

/* -----------------------------------------------------------------------------
 * CHOLMOD Extra Module.
 *
 * Sparse matrix routines.
 *
 * cholmod_spinv		sparse inverse (from simplicial Cholesky)
 *
 * Requires the Core module, and three packages: CHOLMOD, AMD and COLAMD.
 * Optionally uses the Supernodal and Partition modules.
 * -------------------------------------------------------------------------- */

#ifndef CHOLMOD_EXTRA_H
#define CHOLMOD_EXTRA_H

// #include <SuiteSparse_config.h>
// #include <cholmod_config.h>
// #include <cholmod_core.h>

#ifndef NPARTITION
// #include <cholmod_partition.h>
#endif

#ifndef NSUPERNODAL
// #include <cholmod_supernodal.h>
#endif

/* -------------------------------------------------------------------------- */
/* cholmod_spinv:  compute the sparse inverse of a sparse matrix              */
/* -------------------------------------------------------------------------- */

cholmod_sparse *cholmod_SIS
(
    /* ---- input ---- */
    cholmod_factor *L,	/* factorization to use */
    /* --------------- */
    cholmod_common *Common
) ;

cholmod_sparse *cholmod_l_SIS( cholmod_factor *L, cholmod_common *Common ) ;


#endif
