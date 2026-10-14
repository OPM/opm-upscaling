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

/** @file CollectiveCheck.hpp
    @brief Make all MPI ranks agree on whether a piece of work failed.
 */

#ifndef OPM_UPSCALING_COLLECTIVE_CHECK_HPP
#define OPM_UPSCALING_COLLECTIVE_CHECK_HPP

#include <functional>
#include <stdexcept>

namespace Opm {

    //! \brief Exception thrown on every rank when at least one rank failed.
    //!
    //! The original error messages have already been reported (with rank
    //! number) on standard error by the time this exception is thrown, so
    //! handlers typically just need to exit with a failure code.
    class CollectiveFailure : public std::runtime_error
    {
    public:
        using std::runtime_error::runtime_error;
    };

    //! \brief Run a unit of work on all ranks and agree on its outcome.
    //!
    //! Must be called collectively by all ranks in MPI_COMM_WORLD.
    //!
    //! Each rank runs \p work, catching any exception (including non-
    //! std::exception ones).  The ranks then exchange a success/failure
    //! flag in a single reduction.  If any rank failed, then the failing
    //! ranks report their error message on standard error and *every* rank
    //! throws CollectiveFailure.  This lets all ranks leave the parallel
    //! region together instead of some ranks exiting while others are
    //! blocked in a barrier, receive or reduction.
    //!
    //! The reduction is the first point at which a failed rank can be told
    //! apart from a slow one, so ranks that finish early will wait here
    //! for the slowest rank before the failure is known.
    //!
    //! In a run with a single rank (or a build without MPI) \p work is
    //! simply called and its exceptions propagate unchanged.
    //!
    //! \param[in] work Callable to run on this rank.  It must not itself
    //!    make collective MPI calls that could be skipped by an exception.
    void collectiveCheck(const std::function<void()>& work);

    //! \brief Terminate all ranks because this rank cannot continue.
    //!
    //! This is a last-resort safety net for failures that happen on one
    //! rank outside of collectiveCheck().  A rank that simply leaves the
    //! program in that situation will leave the other ranks waiting for it
    //! in some barrier, receive, or reduction.  It is not a substitute for
    //! collectiveCheck(), since an abort gives no chance for an orderly
    //! shutdown and the MPI runtime typically prints diagnostics of its own.
    //!
    //! In a run with more than one rank, this flushes the standard output
    //! and error streams and calls MPI_Abort() on MPI_COMM_WORLD with the
    //! given exit code.  The function does not return in that case.  In a
    //! run with a single rank (or a build without MPI) the function returns
    //! immediately, and the caller should exit in the usual way.
    //!
    //! Call this *after* reporting the reason for the failure.
    //!
    //! \param[in] exitCode Exit code for the failed run.
    void abortParallelRun(int exitCode);

} // namespace Opm

#endif // OPM_UPSCALING_COLLECTIVE_CHECK_HPP
