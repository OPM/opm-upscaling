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

#include <config.h>

#include <opm/upscaling/CollectiveCheck.hpp>

#include <exception>
#include <functional>
#include <iostream>
#include <string_view>

#if HAVE_MPI
#include <mpi.h>

#include <sstream>
#endif  // HAVE_MPI

#if HAVE_MPI

namespace {

    void report(const int rank, std::string_view what)
    {
        std::cerr << "Rank " << rank << ": " << what << '\n';
    }

    int numberOfRanks()
    {
        int initialised = 0;
        int finalised = 0;
        MPI_Initialized(&initialised);
        MPI_Finalized(&finalised);

        int size = 1;
        if (initialised && !finalised) {
            MPI_Comm_size(MPI_COMM_WORLD, &size);
        }

        return size;
    }

} // Anonymous namespace

void Opm::collectiveCheck(const std::function<void()>& work)
{
    const auto size = numberOfRanks();

    if (size == 1) {
        work();
        return;
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int failed = 0;
    try {
        work();
    }
    catch (const std::exception& e) {
        failed = 1;
        report(rank, e.what());
    }
    catch (const char* msg) {
        // Some legacy code throws C strings.
        failed = 1;
        report(rank, msg);
    }
    catch (...) {
        failed = 1;
        report(rank, "Unknown exception");
    }

    int numFailed = 0;
    MPI_Allreduce(&failed, &numFailed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (numFailed > 0) {
        std::ostringstream msg;
        msg << numFailed << " of " << size
            << " MPI ranks failed.  Shutting down.";

        if (rank == 0) {
            std::cerr << msg.str() << '\n';
        }

        throw CollectiveFailure { msg.str() };
    }
}

void Opm::abortParallelRun(const int exitCode)
{
    if (numberOfRanks() == 1) {
        return;
    }

    // MPI_Abort() does not necessarily flush the C++ streams.
    std::cout.flush();
    std::cerr.flush();

    MPI_Abort(MPI_COMM_WORLD, exitCode);
}

#else  // !HAVE_MPI

void Opm::collectiveCheck(const std::function<void()>& work)
{
    work();
}

void Opm::abortParallelRun(const int)
{}

#endif  // HAVE_MPI
