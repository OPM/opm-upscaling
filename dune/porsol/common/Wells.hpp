/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLS_HEADER_INCLUDED
#define OPM_WELLS_HEADER_INCLUDED


// Forward declaration.
namespace Dune
{
    class EclipseGridParser;
}


namespace Opm
{



    /// A class designed to encapsulate a set of rate- or
    /// pressure-controlled wells. Current implementation is a
    /// skeleton only.
    class Wells
    {
    public:
        void init(const Dune::EclipseGridParser& parser);

        // Well-centric interface.
        int numWells() const;
        enum WellType { Injector, Producer };
        WellType type(int wellnum) const;
        enum WellControl { Rate, Pressure };
        WellControl control(int wellnum) const;
        double target(int wellnum) const;
        int numPerforations(int wellnum) const;
        int wellCell(int wellnum, int perfnum) const;
        double wellIndex(int wellnum, int perfnum) const;
        double pressureDelta(int wellnum, int perfnum) const;

        // Cell-centric interface.
        double wellOutflowRate(int cell) const;
    private:
    };


    // ------------ Method implementations --------------

    inline void Wells::init(const Dune::EclipseGridParser& parser)
    {
    }

    inline int Wells::numWells() const
    {
        return 0;
    }

    inline Wells::WellType Wells::type(int wellnum) const
    {
        return Injector;
    }

    inline Wells::WellControl Wells::control(int wellnum) const
    {
        return Rate;
    }

    inline double Wells::target(int wellnum) const
    {
        return 0.0;
    }

    inline int Wells::numPerforations(int wellnum) const
    {
        return 0;
    }

    inline int Wells::wellCell(int wellnum, int perfnum) const
    {
        return -1;
    }

    inline double Wells::wellIndex(int wellnum, int perfnum) const
    {
        return 0.0;
    }

    inline double Wells::pressureDelta(int wellnum, int perfnum) const
    {
        return 0.0;
    }

    inline double Wells::wellOutflowRate(int cell) const
    {
        return 0.0;
    }


} // namespace Opm


#endif // OPM_WELLS_HEADER_INCLUDED
