/*
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

// This test is meaningful when run on several MPI ranks, e.g.,
//
//   mpirun -np 4 test_collectivecheck
//
// In that case a failure on a subset of the ranks must make *every* rank
// throw Opm::CollectiveFailure, and none of the ranks may hang.  It also
// runs (and passes) in a serial run.  There, the original exception must
// propagate unchanged.

#include <config.h>

#define BOOST_TEST_MODULE CollectiveCheck

#include <boost/test/unit_test.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <opm/upscaling/CollectiveCheck.hpp>

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace {

    struct MPIFixture
    {
        MPIFixture()
        {
            auto& suite = boost::unit_test::framework::master_test_suite();
            int argc = suite.argc;
            char** argv = suite.argv;

            helper = &Dune::MPIHelper::instance(argc, argv);
        }

        Dune::MPIHelper* helper{nullptr};
    };

    int rank()
    {
        return Dune::MPIHelper::getCommunication().rank();
    }

    int size()
    {
        return Dune::MPIHelper::getCommunication().size();
    }

    // Check that the communicator is still usable and that all ranks are
    // still in step with each other.
    void checkStillInSync()
    {
        BOOST_CHECK_EQUAL(Dune::MPIHelper::getCommunication().sum(1), size());
    }

} // Anonymous namespace

BOOST_GLOBAL_FIXTURE(MPIFixture);

BOOST_AUTO_TEST_CASE(Success_On_All_Ranks)
{
    int count = 0;
    BOOST_CHECK_NO_THROW(Opm::collectiveCheck([&count]() { ++count; }));
    BOOST_CHECK_EQUAL(count, 1);

    checkStillInSync();
}

BOOST_AUTO_TEST_CASE(Failure_On_Last_Rank)
{
    const auto failingRank = size() - 1;

    auto work = [failingRank]() {
        if (rank() == failingRank) {
            throw std::runtime_error { "Injected failure" };
        }
    };

    if (size() > 1) {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), Opm::CollectiveFailure);
    }
    else {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), std::runtime_error);
    }

    checkStillInSync();
}

BOOST_AUTO_TEST_CASE(Failure_On_Rank_Zero_Only)
{
    auto work = []() {
        if (rank() == 0) {
            throw std::runtime_error { "Injected failure" };
        }
    };

    if (size() > 1) {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), Opm::CollectiveFailure);
    }
    else {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), std::runtime_error);
    }

    checkStillInSync();
}

BOOST_AUTO_TEST_CASE(Failure_On_All_Ranks)
{
    auto work = []() { throw std::logic_error { "Injected failure" }; };

    if (size() > 1) {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), Opm::CollectiveFailure);
    }
    else {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), std::logic_error);
    }

    checkStillInSync();
}

BOOST_AUTO_TEST_CASE(Failure_Not_Derived_From_Std_Exception)
{
    // Legacy code throws C strings and other non-exception objects.
    auto work = []() {
        if (rank() == size() - 1) {
            throw "C string exception";
        }
    };

    if (size() > 1) {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), Opm::CollectiveFailure);
    }
    else {
        BOOST_CHECK_THROW(Opm::collectiveCheck(work), const char*);
    }

    checkStillInSync();
}

BOOST_AUTO_TEST_CASE(Abort_On_Unilateral_Failure)
{
    // This case terminates the entire run, so it only does anything when
    // explicitly requested through the test program's own command line
    // option "--abort-test".  The CTest driver runs it as an expected
    // failure on several ranks.
    const auto& suite = boost::unit_test::framework::master_test_suite();

    bool requested = false;
    for (int i = 1; i < suite.argc; ++i) {
        requested = requested || (std::string { suite.argv[i] } == "--abort-test");
    }

    if (! requested) {
        return;
    }

    if (rank() == size() - 1) {
        // A single rank gives up on its own, outside of collectiveCheck().
        // This must take down the other ranks too.  In a run with a single
        // rank, the call returns and the test case is reported as a failure.
        Opm::abortParallelRun(EXIT_FAILURE);

        BOOST_FAIL("abortParallelRun() must not return in a parallel run");
    }
    else {
        // The other ranks wait for the failed rank.  If they are not taken
        // down, then the test hangs and times out.
        Dune::MPIHelper::getCommunication().barrier();

        BOOST_FAIL("Should have been aborted while waiting");
    }
}

BOOST_AUTO_TEST_CASE(Several_Phases_In_Sequence)
{
    // A failure in one phase must not leave the ranks out of step for
    // later phases.
    BOOST_CHECK_NO_THROW(Opm::collectiveCheck([]() {}));

    BOOST_CHECK_THROW(Opm::collectiveCheck([]() {
        if (rank() == size() / 2) {
            throw std::runtime_error { "Injected failure" };
        }
    }), std::runtime_error);   // CollectiveFailure is-a runtime_error

    BOOST_CHECK_NO_THROW(Opm::collectiveCheck([]() {}));

    checkStillInSync();
}
