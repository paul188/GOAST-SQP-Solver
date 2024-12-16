#include <iostream>
#include <chrono>
#include <ctime>
#include <string>

#include <goast/Quads.h>
#include <goast/Core.h>
#include <goast/DiscreteShells.h>

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

/**/
int main(int argc, char *argv[])
{

try{
    std::cerr << "=================================================================================" << std::endl;
    std::cerr << "SIMPLE BENDING MEMBRANE SIMULATION FOR QUADS" << std::endl;
    std::cerr << "=================================================================================" << std::endl << std::endl;

    QuadMesh qm;
    qm.read_overwrite_quad_mesh("../../data/plate/quadMesh2.ply");

    TriMesh plate;
    plate = qm.makeTriMeshCentroid();

    int num_barycenters = qm.num_faces();

    OpenMesh::IO::write_mesh(plate, "plateCentroid.ply");

    VectorType weights = VectorType::Ones(plate.n_edges());

    MeshTopologySaver plateTopol( plate );

    // Only one initial reference and deformed geometry -> only care for vertex positions
    VectorType plateGeomRef;
    VectorType plateGeomDef;

    getGeometry( plate, plateGeomRef );
    getGeometry( plate, plateGeomDef );

    std::vector<int> bdryMask;
    for(int i = 0; i < plateTopol.getNumVertices(); i++)
    {
        VecType coords;
        getXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);

        if(coords[0] < 0.05 || coords[0] > 0.95){
            bdryMask.push_back(i);
        }

        if(coords[0] == 0.5){
            coords[2] -= 0.1;
        }

        coords[0]-=0.5;
        coords[0]*=0.8;
        coords[0]+=0.4;

        setXYZCoord<VectorType, VecType>( plateGeomDef, coords, i);
    }

    setGeometry(plate, plateGeomDef);

    OpenMesh::IO::write_mesh(plate, "plate_def.ply");

    auto def_geometries = std::make_shared<std::vector<VectorType>>();
    
    SimpleBendingEnergy<DefaultConfigurator> E_bend(plateTopol, plateGeomRef, true, weights);
    SimpleBendingGradientDef<DefaultConfigurator> DE_bend(plateTopol, plateGeomRef, weights);
    SimpleBendingHessianDef<DefaultConfigurator> D2E_bend(plateTopol, plateGeomRef, weights);

    NonlinearMembraneEnergy<DefaultConfigurator> E_membrane(plateTopol, plateGeomRef, true);
    NonlinearMembraneGradientDef<DefaultConfigurator> DE_membrane(plateTopol, plateGeomRef);
    NonlinearMembraneHessianDef<DefaultConfigurator> D2E_membrane(plateTopol, plateGeomRef);

    RealType factor_membrane = 10000.0;
    RealType factor_bending = 1.0;

    VectorType factors(2);
    factors[0] = factor_membrane;
    factors[1] = factor_bending;

    AdditionOp<DefaultConfigurator> E_tot( factors, E_membrane, E_bend);
    AdditionGradient<DefaultConfigurator> DE_tot( factors, DE_membrane, DE_bend);
    AdditionHessian<DefaultConfigurator> D2E_tot(factors, D2E_membrane, D2E_bend);

    // set outer optimization parameters
    OptimizationParameters<DefaultConfigurator> optPars;
    optPars.setNewtonIterations( 10000 );
    optPars.setQuietMode( SHOW_ALL );
    VectorType initialization = plateGeomDef;

    std::cerr<< "Startting Newton Linesearch..."<<std::endl;
    LineSearchNewton<DefaultConfigurator> NLS( E_tot, DE_tot, D2E_tot, optPars);
    NLS.setBoundaryMask( bdryMask );
    NLS.solve( initialization, plateGeomDef);

    // saving
    setGeometry( plate, plateGeomDef );
    OpenMesh::IO::write_mesh(plate, "bendingFoldSol_withNewton2.ply");

}catch(std::exception &e){
    std::cerr << "Exception caught: " << e.what() << std::endl;
}
return 0;
}