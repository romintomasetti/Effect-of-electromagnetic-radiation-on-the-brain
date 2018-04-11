#include "vtlSPoints_romin.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace vtl_romin;

SPoints::SPoints() : id(-1), o(), np1(), np2(), dx()
{
}

VTL_API_ROMIN std::ostream &
vtl_romin::operator<<(std::ostream &out, vtl_romin::SPoints const &obj)
{
    out << "SPoints: np=" << obj.np() << "\n";
    out << "\to = " << obj.o << '\n';
    out << "\tL = " << obj.L() << '\n';
    out << "\trange = ["
              << obj.np1[0] << '-' << obj.np2[0] << ", "
              << obj.np1[1] << '-' << obj.np2[1] << ", "
              << obj.np1[2] << '-' << obj.np2[2] << "]\n";
    out << "\tdx = " << obj.dx << '\n';
    return out;
}

SPoints 
SPoints::split(int numprocs, int myid)
{
    SPoints subgrid = *this; // copy the whole grid
    subgrid.id = myid; // assign rank# to the sub-grid

    int sdir = 0; // split direction

    // split grid along 'sdir'
    int nfloor = np()[sdir] / numprocs;
    int rem = np()[sdir] % numprocs;
    assert(np()[sdir] == nfloor * numprocs + rem);

    // build the list if idx separating each subdomain
    std::vector<int> nps(numprocs + 1);
    nps[0] = np1[sdir];
    for (int i = 0; i < numprocs + 1; ++i)
    {
        nps[i + 1] = nps[i] + nfloor;
        if (i < rem)
            nps[i + 1]++;
    }

    for (int i = 0; i < 3; ++i)
    {
        subgrid.np1[i] = np1[i];
        subgrid.np2[i] = np2[i];
    }
    subgrid.np1[sdir] = nps[myid];
    subgrid.np2[sdir] = nps[myid + 1] - 1;

    // add extra layers from direct neighborhood
    if (myid != 0)
        subgrid.np1[sdir] -= 1;

    return subgrid;
}
