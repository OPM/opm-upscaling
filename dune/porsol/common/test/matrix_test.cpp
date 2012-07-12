/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <opm/core/utility/have_boost_redef.hpp>

#if defined(HAVE_DYNAMIC_BOOST_TEST)
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing


#define BOOST_TEST_MODULE MatrixTests
#include <boost/test/unit_test.hpp>

#include "../Matrix.hpp"



BOOST_AUTO_TEST_CASE(copy_assignment_tests)
{
    using namespace Dune;
    OwnCMatrix m1(1,1,(const double*)0);
    m1(0,0) = 3.14;
    OwnCMatrix m2(m1);
    BOOST_CHECK(m1 == m2);
    BOOST_CHECK_EQUAL(m1(0,0), m2(0,0));
    double storage_m3[1];
    SharedCMatrix m3(1,1,storage_m3);
    m3 = m2;
    BOOST_CHECK(m2 == m3);
    double storage_m4[1];
    SharedCMatrix m4(1,1,storage_m4);
    m4 = m3;
    BOOST_CHECK(m3 == m4);
    BOOST_CHECK_EQUAL(m4.data(), storage_m3); // Note this behaviour (identical types)...
    double storage_m5[1];
    SharedFortranMatrix m5(1,1,storage_m5);
    m5 = m3;
    BOOST_CHECK_EQUAL(m5.data(), storage_m5); // ... compared to this (different types).
}
