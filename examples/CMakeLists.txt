project(examples)

add_executable(multires_example multires_example.cpp)
add_executable(interpolationShells geodesic_interpolation_shells.cpp)
add_executable(extrapolationShells geodesic_extrapolation_shells.cpp)
add_executable(dissipationModes dissipationModes.cpp)
add_executable(decimation decimation.cpp)
add_executable(elasticMean elasticMean.cpp)
add_executable(parallelTransport parallelTransport.cpp)
add_executable(deformationTransfer deformationTransfer.cpp)
add_executable(nonlinearSpline nonlinearSpline.cpp)
add_executable(sectionalCurvature sectionalCurvature.cpp)
add_executable(registration registration.cpp)
add_executable(DirichletEnergy DirichletEnergy.cpp)
add_executable(parametricMCM parametricMCM.cpp)
add_executable(gravitationBending GravitationBending/gravitationBending.cpp)
add_executable(bendingFold SimpleFold/bendingFold.cpp)
add_executable(bendingFoldParamInterpolation SingleArcCrease/bendingFoldParamInterpolation.cpp )
add_executable(bendingFoldParamDirichlet SingleArcCrease/bendingFoldParamDirichlet.cpp)
add_executable(bendingSkewedCompression SkewedCompression/bendingSkewedCompression.cpp)
add_executable(bendingSmallArc SmallArc/bendingSmallArc.cpp)
add_executable(bendingSmallArcInterpolation SmallArc/bendingSmallArcInterpolation.cpp)
#add_executable(SimpleFoldOptimization2 FoldOptimization/SimpleFoldOptimization2.cpp)
add_executable(QuadMeshTest QuadMeshExamples/QuadMeshTest.cpp)
#add_executable(SimpleFoldOptimization FoldOptimization/SimpleFoldOptimization.cpp)
add_executable(QuadMeshTest2 QuadMeshExamples/QuadMeshTest2.cpp)
add_executable(ToyExample FoldOptimization/ToyExample.cpp)
add_executable(ToyExampleNewton FoldOptimization/ToyExampleNewton.cpp)
# link libraries
target_link_libraries(ToyExampleNewton PRIVATE GOAST::GOAST)

target_link_libraries(ToyExample PRIVATE GOAST::GOAST)

target_link_libraries(QuadMeshTest PRIVATE GOAST::GOAST)

target_link_libraries(QuadMeshTest2 PRIVATE GOAST::GOAST)

target_link_libraries(multires_example PRIVATE GOAST::GOAST)

target_link_libraries(interpolationShells PRIVATE GOAST::GOAST)

target_link_libraries(extrapolationShells PRIVATE GOAST::GOAST)

target_link_libraries(dissipationModes PRIVATE GOAST::GOAST)

target_link_libraries(decimation PRIVATE GOAST::GOAST)

target_link_libraries(elasticMean PRIVATE GOAST::GOAST)

target_link_libraries(parallelTransport PRIVATE GOAST::GOAST)

target_link_libraries(deformationTransfer PRIVATE GOAST::GOAST)

target_link_libraries(nonlinearSpline PRIVATE GOAST::GOAST)

target_link_libraries(sectionalCurvature PRIVATE GOAST::GOAST)

target_link_libraries(registration PRIVATE GOAST::GOAST)

target_link_libraries(DirichletEnergy PRIVATE GOAST::GOAST)

target_link_libraries(parametricMCM PRIVATE GOAST::GOAST)

target_link_libraries(gravitationBending PRIVATE GOAST::GOAST)

target_link_libraries(bendingFold PRIVATE GOAST::GOAST)

target_link_libraries(bendingFoldParamInterpolation PRIVATE GOAST::GOAST)

target_link_libraries(bendingSkewedCompression PRIVATE GOAST::GOAST)

target_link_libraries(bendingSmallArc PRIVATE GOAST::GOAST)

target_link_libraries(bendingFoldParamDirichlet PRIVATE GOAST::GOAST)

target_link_libraries(bendingSmallArcInterpolation PRIVATE GOAST::GOAST)

#target_link_libraries(SimpleFoldOptimization PRIVATE GOAST::GOAST)

#target_link_libraries(SimpleFoldOptimization2 PRIVATE GOAST::GOAST)

# Phasefield example currently requires Ipopt
if (GOAST_WITH_IPOPT)
    add_executable(Phasefield Phasefield.cpp)
    target_link_libraries(Phasefield PRIVATE GOAST)
    target_link_libraries(Phasefield PRIVATE GOAST::All)
endif ()

# NRIC Examples
add_executable(NRICGeodesics NRICGeodesics.cpp)
target_link_libraries(NRICGeodesics PRIVATE GOAST::GOAST)

add_executable(NRICConstrictions NRICConstrictions.cpp)
target_link_libraries(NRICConstrictions PRIVATE GOAST::GOAST)

add_executable(NRICFolding NRICFolding.cpp)
target_link_libraries(NRICFolding PRIVATE GOAST::GOAST)

add_executable(NRICMean NRICMean.cpp)
target_link_libraries(NRICMean PRIVATE GOAST::GOAST)

add_executable(NRICIsometricMode NRICIsometricMode.cpp)
target_link_libraries(NRICIsometricMode PRIVATE GOAST::All)

add_executable(nonlinearProjection nonlinearProjection.cpp)
target_link_libraries(nonlinearProjection PRIVATE GOAST::GOAST)