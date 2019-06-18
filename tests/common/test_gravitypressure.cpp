/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

  This file is part of The Open Porous Media project (OPM).

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

#include <config.h>

#define BOOST_TEST_MODULE UpscaleHelperGravityPressure

// Enable intrusive unit-testing of private property
//
//   ::Opm::RelPermUpscaleHelper::dP
//
#define UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP 1

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif

#include <opm/upscaling/RelPermUtils.hpp>
#undef  UNITTEST_TRESPASS_PRIVATE_PROPERTY_DP

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <opm/parser/eclipse/Units/Units.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

namespace {                     // Case utilities
    namespace GridInput {
        namespace detail {
            std::string RUNSPEC()
            {
                return R"~~(
RUNSPEC

DIMENS
  1 1 10 /

OIL
WATER

METRIC
)~~";
            }

            // Note: Must be COORD + ZCORN because of details of
            // implementation of Opm::UpscalerBase--notably the
            // representation of the rock properties.
            std::string GRID()
            {
                return R"~~(
GRID

COORD
  0 0 0   0 0 0
  1 0 0   1 0 0
  0 1 0   0 1 0
  1 1 0   1 1 0
/

ZCORN
  4*0.0
  4*0.1

  4*0.1
  4*0.2

  4*0.2
  4*0.3

  4*0.3
  4*0.4

  4*0.4
  4*0.5

  4*0.5
  4*0.6

  4*0.6
  4*0.7

  4*0.7
  4*0.8

  4*0.8
  4*0.9

  4*0.9
  4*1.0
/

PERMX
  10*100
/

PERMY
  10*100
/

PERMZ
  10*10
/

PORO
  10*0.3
/
)~~";
            }
        } // detail

        std::string basic()
        {
            return detail::RUNSPEC() + detail::GRID();
        }
    } // GridInput

    namespace Options
    {
        using Opt = std::map<std::string, std::string>;

        namespace detail {
            Opt tesselateBasic()
            {
                return
                {
                    {"linsolver_tolerance",       "1e-12"},
                    {"linsolver_verbosity",           "0"},
                    {"linsolver_type",                "3"},
                    {"linsolver_max_iterations",      "0"},
                    {"linsolver_smooth_steps",        "1"},
                    {"linsolver_prolongate_factor", "1.0"},
                    {"minPerm",                   "1e-12"}, // mD
                };
            }

            Opt fluidSystem()
            {
                return
                {
                    { "fluids", "ow" },
                };
            }

            Opt boundaryCondition()
            {
                return
                {
                    { "bc", "fixed" },
                };
            }

            Opt gravityOn()
            {
                return
                {
                    { "gravity", "10.0" },
                };
            }

            Opt gravityOff()
            {
                return
                {
                    { "gravity", "Zero" },
                };
            }

            Opt gravityEquator()
            {
                return
                {
                    { "gravity", "9.80665" },
                };
            }

            Opt densityNormal()
            {
                return
                {
                    { "waterDensity", "1.0" }, // g/cm3
                    { "oilDensity"  , "0.8" }, // g/cm3
                };
            }

            Opt densityWrong()
            {
                // String->double conversion internally effected by
                // std::strtod() which returns zero (0.0) when unable to
                // convert an input string.  No error checking at time of
                // writing.
                //
                // These options therefore behave as if set to "0.0".
                return
                {
                    { "waterDensity", "One" }, // g/cm3
                    { "oilDensity"  , "Zero Point Eight" }, // g/cm3
                };
            }

            Opt mergeOptions(Opt&&                      into,
                             std::initializer_list<Opt> others)
            {
                for (const auto& other : others) {
                    into.insert(other.begin(), other.end());
                }

                return into;
            }
        } // detail

        Opt baseCase()
        {
            return detail::mergeOptions(detail::tesselateBasic(),
                                        { detail::boundaryCondition(),
                                          detail::fluidSystem(),
                                          detail::gravityOn(),
                                          detail::densityNormal() });
        }

        Opt noGravity()
        {
            return detail::mergeOptions(detail::tesselateBasic(),
                                        { detail::boundaryCondition(),
                                          detail::fluidSystem(),
                                          detail::gravityOff(),
                                          detail::densityNormal() });
        }

        Opt noGravityWrong()
        {
            return detail::mergeOptions(detail::tesselateBasic(),
                                        { detail::boundaryCondition(),
                                          detail::fluidSystem(),
                                          detail::gravityOn(),
                                          detail::densityWrong() });
        }

        Opt equator()
        {
            return detail::mergeOptions(detail::tesselateBasic(),
                                        { detail::boundaryCondition(),
                                          detail::fluidSystem(),
                                          detail::gravityEquator(),
                                          detail::densityNormal() });
        }
    } // Options

    namespace TestCase
    {
        namespace Results
        {
            struct Computed
            {
                explicit Computed(std::vector<double>&& x);

                std::vector<double> data;
            };

            struct Expected
            {
                explicit Expected(std::vector<double>&& x);

                std::vector<double> data;
            };

            Computed::Computed(std::vector<double>&& x)
                : data(std::move(x))
            {}

            Expected::Expected(std::vector<double>&& x)
                : data(std::move(x))
            {}

            void compare(const Computed& computed,
                         const Expected& expected)
            {
                const auto& c = computed.data;
                const auto& e = expected.data;

                BOOST_REQUIRE_EQUAL(c.size(), e.size());

                for (decltype(c.size())
                         i = 0, n = c.size(); i < n; ++i)
                {
                    BOOST_CHECK_CLOSE(c[i], e[i], 100.0e-6);
                }
            }
        } // Results

        namespace Expect
        {
            namespace detail
            {
                std::vector<double>::size_type
                numCells()
                {
                    return 10;
                }

                double offset()
                {
                    return 5.0;
                }

                double g10()
                {
                    return 10.0;
                }

                double gEquator()
                {
                    return 9.80665;
                }

                double dRho()
                {
                    using namespace Opm;

                    return unit::convert::
                        from(1.0 - 0.8,
                             prefix::milli*unit::kilogram
                             / unit::cubic(prefix::centi*unit::meter));
                }

                std::vector<double>
                gravityPressure(const double grav)
                {
                    auto dp = std::vector<double>(numCells(), 0.0);

                    if (std::abs(grav) > 0.0) {
                        const auto off  = offset();
                        const auto drho = dRho();

                        decltype(dp.size()) i = 0;

                        for (auto& x : dp) {
                            x = grav * drho * ((i + 0.5 - off) / numCells());

                            ++i;
                        }
                    }

                    return dp;
                }
            } // detail

            std::vector<double> noGravity()
            {
                return detail::gravityPressure(0.0);
            }

            std::vector<double> gravity10()
            {
                return detail::gravityPressure(detail::g10());
            }

            std::vector<double> gravityEquator()
            {
                return detail::gravityPressure(detail::gEquator());
            }
        } // Expect

        namespace Compute
        {
            std::vector<double>
            modelDP(const std::string& input,
                    Options::Opt&      options)
            {
                // We're NOT master (rank == 0) in this context.
                //
                // This is to suppress informational/diagnostic output from
                // RelPermUpscaleHelper during grid construction (i.e., when
                // calling tesselateGrid()).
                const auto mpi_rank = 1;

                int m_argc = boost::unit_test::framework::master_test_suite().argc;
                char** m_argv = boost::unit_test::framework::master_test_suite().argv;
                Dune::MPIHelper::instance(m_argc, m_argv);

                Opm::RelPermUpscaleHelper helper{mpi_rank, options};
                {
                    auto parse = Opm::Parser{};
                    auto deck  = parse.parseString(input);

                    const auto minPerm = 0;      // mD
                    const auto maxPerm = 10.0e3; // mD
                    const auto minPoro = 0.1;

                    // Note: sanityCheckInput() extracts and stores private
                    // copies of ZCORN and PERMX from the 'deck'.
                    //
                    // This function *must* be called before tesselateGrid().
                    helper.sanityCheckInput(deck, minPerm, maxPerm, minPoro);

                    // Function tesselateGrid() cannot be called unless the
                    // upscaler's boundary conditions are well defined.
                    helper.setupBoundaryConditions();

                    helper.tesselateGrid(deck);
                }

                helper.calculateCellPressureGradients();

                // Note: 'private' member of RelPermUpscaleHelper.
                return helper.dP;
            }
        } // Compute

        void run(std::function<std::string()>         input,
                 std::function<Options::Opt()>        options,
                 std::function<std::vector<double>()> expect)
        {
            auto computed = std::vector<double>{};
            {
                auto opt = options(); // Live for duration of modelDP().
                computed = Compute::modelDP(input(), opt);
            }

            Results::compare(Results::Computed{std::move(computed)},
                             Results::Expected{expect()});
        }
    } // TestCase
} // Anonymous

// =====================================================================
// Actual test suite below
// =====================================================================

BOOST_AUTO_TEST_SUITE()

BOOST_AUTO_TEST_CASE (BaseCase)
{
    TestCase::run(GridInput::basic,
                  Options::baseCase,
                  TestCase::Expect::gravity10);
}

BOOST_AUTO_TEST_CASE (NoGravity)
{
    TestCase::run(GridInput::basic,
                  Options::noGravity,
                  TestCase::Expect::noGravity);
}

BOOST_AUTO_TEST_CASE (NoGravityWrong)
{
    TestCase::run(GridInput::basic,
                  Options::noGravityWrong,
                  TestCase::Expect::noGravity);
}

BOOST_AUTO_TEST_CASE (GravityEquator)
{
    TestCase::run(GridInput::basic,
                  Options::equator,
                  TestCase::Expect::gravityEquator);
}

BOOST_AUTO_TEST_SUITE_END()
