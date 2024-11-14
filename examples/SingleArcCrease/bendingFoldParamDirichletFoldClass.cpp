/* 
 Simulating a fold that results from bending a plate with clamped 
 boundary conditions. Only NonlinearMembraneEnergy and SimpleBendingEnergy are used.
 The gravitational energy is not used in this example.
 We set the weight of the contribution of the middle edge to the bending energy to 0.
 We thus expect a lot of the bending to happen in the middle of the plate, at the fold.
 Imagine a piece of paper that is folded in the middle in both directions.
 So that the fold does not point into any direction, but the edge just doesnt
 resist bending.
*/

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Core.h>
#include <goast/DiscreteShells.h>
#include <typeinfo>
#include <goast/Folds.h>
//#include <goast/Folds/Fold.h>

//==============================================================================================================
typedef DefaultConfigurator::VectorType VectorType;
typedef DefaultConfigurator::VecType VecType;
typedef DefaultConfigurator::RealType RealType;

/** 
 * \brief Simulation of a simple fold 
 * \author Johannssen
 *
 * We optimize this energy by direct optimization via gradient descent or BFGS.
 * 
 */


int main(int argc, char *argv[])
{
    try{

        std::cerr << "=================================================================================" << std::endl;
        std::cerr << "SIMPLE ARC FOLD SIMULATION WITH BENDING AND MEMBRANE ENERGIES" << std::endl;
        std::cerr << "=================================================================================" << std::endl << std::endl;
        
    // load flat plate and prepare the arc crease
        TriMesh plate;
        OpenMesh::IO::read_mesh(plate, "../../data/plate/rectangle_criss_cross.ply");
        MeshTopologySaver plateTopol( plate );
        VectorType plateGeomRef; // plateGeomRef is the geometry without the arc crease in the mesh
        VectorType plateGeomInitial;
        getGeometry(plate,plateGeomRef);
        getGeometry(plate,plateGeomInitial);

        // Setup the spline and control points
        std::vector<tsReal> ctrlp;
        std::vector<VecType> foldLocations;
        for(int i = 0; i < plateTopol.getNumVertices(); i++){
            VecType coords;
            getXYZCoord<VectorType, VecType>( plateGeomRef, coords, i);

            if(coords[0] == 1.0){
                VecType temp;
                temp[0] = 1.25 - std::pow(coords[1] - 0.5, 2.0);
                temp[1] = coords[1];
                temp[2] = 0.0;
                ctrlp.push_back(temp[0]);
                ctrlp.push_back(temp[1]);
                ctrlp.push_back(temp[2]);
            }
        }

        size_t numControlPoints = ctrlp.size() / 3.0;

        tsReal *ctrlp_arr = new tsReal[ctrlp.size()];

        for(int i = 0; i < ctrlp.size(); i++){
            ctrlp_arr[i] = ctrlp[i];
        }

        tsBSpline tspline;
        tsStatus status;

        ts_bspline_new(numControlPoints,3,3,TS_CLAMPED, &tspline, &status);
        ts_bspline_interpolate_cubic_natural(ctrlp_arr, numControlPoints,3,&tspline,&status);

        tsDeBoorNet net;
        ts_bspline_eval(&tspline, 0.0, &net, &status);
        tsReal *result;
        //ts_deboornet_result(&net, &result, &status);
        std::cout << "Result: " << result[0] << " " << result[1] << " " << result[2] << std::endl;
        // fold( tspline , foldLocations);

        //fold.applyFoldToMeshGeometry(plateTopol, plateGeomRef);
        }
        catch ( BasicException &el ){
            std::cerr << std::endl << "ERROR!! CAUGHT FOLLOWING EXECEPTION: " << std::endl << el.getMessage() << std::endl << std::flush;
    }
    return 0;
}