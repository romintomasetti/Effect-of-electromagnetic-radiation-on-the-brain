#ifndef VTLSPOINTS_H
#define VTLSPOINTS_H

#include "vtl.h"
#include "vtlVec3.h"
#include <map>
#include <string>
#include <vector>


namespace vtl
{
    class SPoints;
    VTL_API std::ostream &operator<<(std::ostream &out, SPoints const &obj);

/**
 * @brief A structured set of points (vtkStructuredPoints in VTK)
 */
class VTL_API SPoints
{
  public:
    int id;     ///< MPI rank of the grid
    Vec3d o;    ///< origin
    Vec3i np1;  ///< starting indices
    Vec3i np2;  ///< ending indices
    Vec3d dx;   ///< spacing
    
    /// scalar values stored at the points
    std::map<std::string, std::vector<double> *> scalars; 
    /// vector values stored at the points
    std::map<std::string, std::vector<double> *> vectors;

    SPoints();
	SPoints split(int numprocs, int myid);
    Vec3d L() const { return Vec3d(np2-np1)*dx; }
	Vec3i np() const { return np2-np1+1; }
	int nbp() const { Vec3i a = np(); return a[0]*a[1]*a[2]; }
    friend VTL_API std::ostream &operator<<(std::ostream &out, SPoints const &obj);
};
}

#endif //VTLSPOINTS_H