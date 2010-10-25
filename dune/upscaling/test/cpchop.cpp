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

#include <dune/common/EclipseGridParser.hpp>
#include <dune/common/param/ParameterGroup.hpp>
#include <cstdlib>

using namespace Dune;


template <typename T>
void outputFilteredImpl(std::ostream& os,
                        const std::vector<int>& nto,
                        const std::vector<T>& field,
                        const std::string& keyword)
{
    os << keyword << '\n';
    int sz = nto.size();
    T last = std::numeric_limits<T>::max();
    int repeats = 0;
    for (int i = 0; i < sz; ++i) {
        T val = field[nto[i]];
        if (val == last) {
            ++repeats;
        } else {
            if (repeats == 1) {
                os << last << '\n';
            } else if (repeats > 1) {
                os << repeats << '*' << last << '\n';
            }
            last = val;
            repeats = 1;
        }
    }
    if (repeats == 1) {
        os << last << '\n';
    } else if (repeats > 1) {
        os << repeats << '*' << last << '\n';
    }
    os << "/\n\n";
}

void outputFilteredInt(std::ostream& os,
                       const EclipseGridParser& parser,
                       const std::vector<int>& nto,
                       const std::string& keyword)
{
    if (parser.hasField(keyword)) {
        const std::vector<int>& field = parser.getIntegerValue(keyword);
        outputFilteredImpl(os, nto, field, keyword);
    }
}

void outputFilteredDouble(std::ostream& os,
                          const EclipseGridParser& parser,
                          const std::vector<int>& nto,
                          const std::string& keyword)
{
    if (parser.hasField(keyword)) {
        const std::vector<double>& field = parser.getFloatingPointValue(keyword);
        outputFilteredImpl(os, nto, field, keyword);
    }
}



int main(int argc, char** argv)
{
    parameter::ParameterGroup param(argc, argv);
    std::string gridfilename = param.get<std::string>("gridfilename");
    std::cout << "Parsing grdecl file." << std::endl;
    EclipseGridParser parser(gridfilename);
    const int* dims = &parser.getSPECGRID().dimensions[0];
    std::cout << "Parsed grdecl file with dimensions (" << dims[0] << ", " << dims[1] << ", " << dims[2] << ")" << std::endl;

    // The cells with i coordinate in [imin, imax) are included, similar for j.
    int imin = param.getDefault("imin", 0);
    int imax = param.getDefault("imax", dims[0]);
    int jmin = param.getDefault("jmin", 0);
    int jmax = param.getDefault("jmax", dims[1]);
    int new_dims[3] = { imax - imin, jmax - jmin, -1 };

    // Filter the coord field
    const std::vector<double>& COORD = parser.getFloatingPointValue("COORD");
    int num_coord = COORD.size();
    if (num_coord != 6*(dims[0] + 1)*(dims[1] + 1)) {
        std::cerr << "Error! COORD size (" << COORD.size() << ") not consistent with SPECGRID\n";
        return EXIT_FAILURE;
    }
    int num_new_coord = 6*(new_dims[0] + 1)*(new_dims[1] + 1);
    std::vector<double> new_COORD(num_new_coord, 1e100);
    for (int j = jmin; j < jmax + 1; ++j) {
        for (int i = imin; i < imax + 1; ++i) {
            int pos = (dims[0] + 1)*j + i;
            int new_pos = (new_dims[0] + 1)*(j-jmin) + (i-imin);
            std::copy(COORD.begin() + 6*pos, COORD.begin() + 6*(pos + 1), new_COORD.begin() + 6*new_pos);
        }
    }

    // Get the z limits, check if they must be changed to make a shoe-box.
    // This means that zmin must be greater than or equal to the highest
    // coordinate of the bottom surface, while zmax must be less than or
    // equal to the lowest coordinate of the top surface.
    double zmin = param.getDefault("zmin", -1e100);
    double zmax = param.getDefault("zmax", 1e100);
    int layersz = 8*dims[0]*dims[1];
    const std::vector<double>& ZCORN = parser.getFloatingPointValue("ZCORN");
    int num_zcorn = ZCORN.size();
    if (num_zcorn != layersz*dims[2]) {
        std::cerr << "Error! ZCORN size (" << ZCORN.size() << ") not consistent with SPECGRID\n";
        return EXIT_FAILURE;
    }
    double botmax = *std::max_element(ZCORN.begin(), ZCORN.begin() + layersz/2);
    double topmin = *std::min_element(ZCORN.begin() + dims[2]*layersz - layersz/2,
                                 ZCORN.begin() + dims[2]*layersz);
    zmin = std::max(zmin, botmax);
    zmax = std::min(zmax, topmin);

    // We must find the maximum and minimum k value for the given z limits.
    // First, find the first layer with a z-coordinate strictly above zmin.
    int kmin = -1;
    for (int k = 0; k < dims[2]; ++k) {
        double layer_max = *std::max_element(ZCORN.begin() + k*layersz, ZCORN.begin() + (k + 1)*layersz);
        if (layer_max > zmin) {
            kmin = k;
            break;
        }
    }
    // Then, find the last layer with a z-coordinate strictly below zmax.
    int kmax = -1;
    for (int k = dims[2]; k > 0; --k) {
        double layer_min = *std::min_element(ZCORN.begin() + (k - 1)*layersz, ZCORN.begin() + k*layersz);
        if (layer_min < zmax) {
            kmax = k;
            break;
        }
    }
    new_dims[2] = kmax - kmin;

    // Filter the ZCORN field, build mapping from new to old cells.
    std::vector<double> new_ZCORN(8*new_dims[0]*new_dims[1]*new_dims[2], 1e100);
    std::vector<int> new_to_old_cell(new_dims[0]*new_dims[1]*new_dims[2], -1);
    int cellcount = 0;
    int delta[3] = { 1, 2*dims[0], 4*dims[0]*dims[1] };
    int new_delta[3] = { 1, 2*new_dims[0], 4*new_dims[0]*new_dims[1] };
    for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {
                new_to_old_cell[cellcount++] = dims[0]*dims[1]*k + dims[0]*j + i;
                int old_ix = 2*(i*delta[0] + j*delta[1] + k*delta[2]);
                int new_ix = 2*((i-imin)*new_delta[0] + (j-jmin)*new_delta[1] + (k-kmin)*new_delta[2]);
                int old_indices[8] = { old_ix, old_ix + delta[0],
                                       old_ix + delta[1], old_ix + delta[1] + delta[0],
                                       old_ix + delta[2], old_ix + delta[2] + delta[0],
                                       old_ix + delta[2] + delta[1], old_ix + delta[2] + delta[1] + delta[0] };
                int new_indices[8] = { new_ix, new_ix + new_delta[0],
                                       new_ix + new_delta[1], new_ix + new_delta[1] + new_delta[0],
                                       new_ix + new_delta[2], new_ix + new_delta[2] + new_delta[0],
                                       new_ix + new_delta[2] + new_delta[1], new_ix + new_delta[2] + new_delta[1] + new_delta[0] };
                for (int cc = 0; cc < 8; ++cc) {
                    new_ZCORN[new_indices[cc]] = std::min(zmax, std::max(zmin, ZCORN[old_indices[cc]]));
                }
            }
        }
    }

    // Output new versions of SPECGRID, COORD, ZCORN, ACTNUM, PERMX, PORO, SATNUM.
    std::string filebase = param.get<std::string>("filebase");
    std::ofstream out(filebase.c_str());
    if (!out) {
        std::cerr << "Could not open file " << filebase << "\n";
        return EXIT_FAILURE;
    }
    out << "SPECGRID\n" << new_dims[0] << ' ' << new_dims[1] << ' ' << new_dims[2]
        << " 1 F\n/\n\n";
    out << "COORD\n";
    for (int i = 0; i < num_new_coord/6; ++i) {
        for (int j = 0; j < 6; ++j) {
            out << "  " << new_COORD[6*i + j];
        }
        out << '\n';
    }
    out << "/\n\n";
    out << "ZCORN\n";
    int num_new_zcorn = new_ZCORN.size();
    assert(num_new_zcorn%8 == 0);
    for (int i = 0; i < num_new_zcorn/8; ++i) {
        for (int j = 0; j < 8; ++j) {
            out << "  " << new_ZCORN[8*i + j];
        }
        out << '\n';
    }
    out << "/\n\n";

    outputFilteredInt(out, parser, new_to_old_cell, "ACTNUM");
    outputFilteredDouble(out, parser, new_to_old_cell, "PERMX");
    outputFilteredDouble(out, parser, new_to_old_cell, "PERMY");
    outputFilteredDouble(out, parser, new_to_old_cell, "PERMY");
    outputFilteredDouble(out, parser, new_to_old_cell, "PORO");
    outputFilteredInt(out, parser, new_to_old_cell, "SATNUM");
}
