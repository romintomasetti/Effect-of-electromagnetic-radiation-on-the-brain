#ifndef VTLSPOINTS_ROMIN_H
#define VTLSPOINTS_ROMIN_H

#include "vtl_romin.h"
#include "vtlVec3_romin.h"
#include <map>
#include <string>
#include <vector>


namespace vtl_romin
{
    class SPoints;
    VTL_API_ROMIN std::ostream &operator<<(std::ostream &out, SPoints const &obj);

// a structured grid

class VTL_API_ROMIN SPoints
{
  public:
    int id;     ///< rank of the grid
    Vec3d o;    ///< origin
    Vec3i np1;  ///< starting indices
    Vec3i np2;  ///< ending indices
    Vec3d dx;   ///< spacing
    
    std::map<std::string, std::vector<double> *> scalars;
    std::map<std::string, std::vector<double> *> vectors;

    SPoints();
	SPoints split(int numprocs, int myid);
    Vec3d L() const { return Vec3d(np2-np1)*dx; }
	Vec3i np() const { return np2-np1+1; }
	int nbp() const { Vec3i a = np(); return a[0]*a[1]*a[2]; }
    friend VTL_API_ROMIN std::ostream &operator<<(std::ostream &out, SPoints const &obj);
};
}

#endif //VTLSPOINTS_ROMIN_H