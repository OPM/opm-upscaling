/*****************************************************************************
 *   Copyright (C) 2009-2010 by Melanie Darcis                               *
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 *
 * \brief Provides the class with the tabulated values of CO2 for the
 *        benchmark3 problem
 */
#ifndef OPM_BENCHMARK3_CO2TABLES_HH
#define OPM_BENCHMARK3_CO2TABLES_HH

#include <dune/porsol/blackoil/co2fluid/opm/old_material/tabulatedmaterial2.hh>
#include <dune/porsol/blackoil/co2fluid/opm/old_material/tabulatedmaterial2hires.hh>

#include <assert.h>

namespace Opm
{
namespace Benchmark3
{
// the real work is done by some external program which provides
// ready-to-use tables.
#include "benchmark3co2values.inc"
}
}

#endif
