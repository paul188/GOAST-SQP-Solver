#pragma once
#include <goast/Core.h>

template <typename ConfiguratorType>
class MyObjectFactory {
    public:
        MyObjectFactory(){
            factors.resize(2);
            factors << 1.0, 1.0;
        }
        ~MyObjectFactory(){}

        AdditionOp<ConfiguratorType> produceE(const MeshTopologySaver &plateTopol, const VectorType &plateGeomRef, const bool ActiveShellIsDeformed, const VectorType &edge_weights){
            auto *E_bend = new SimpleBendingEnergy<ConfiguratorType>(plateTopol, plateGeomRef, ActiveShellIsDeformed, edge_weights);
            auto *E_membrane = new NonlinearMembraneEnergy<ConfiguratorType>(plateTopol, plateGeomRef, ActiveShellIsDeformed);
            AdditionOp<ConfiguratorType> temp(factors, *E_bend, *E_membrane);
            // Move E_bend and E_membrane to the AdditionOp constructor, dereferencing them
            return temp;
        }

        AdditionGradient<ConfiguratorType> produceDE(const MeshTopologySaver &plateTopol, const VectorType &plateGeomRef, const VectorType &edge_weights){
            auto *DE_bend = new SimpleBendingGradientDef<ConfiguratorType>(plateTopol, plateGeomRef, edge_weights);
            auto *DE_membrane = new NonlinearMembraneGradientDef<ConfiguratorType>(plateTopol, plateGeomRef);
            AdditionGradient<ConfiguratorType> temp(factors, *DE_bend, *DE_membrane);
            return temp;
        }

        AdditionHessian<ConfiguratorType> produceD2E(const MeshTopologySaver &plateTopol, const VectorType &plateGeomRef, const VectorType &edge_weights){
            auto *D2E_bend = new SimpleBendingHessianDef<ConfiguratorType>(plateTopol, plateGeomRef, edge_weights);
            auto *D2E_membrane = new NonlinearMembraneHessianDef<ConfiguratorType>(plateTopol, plateGeomRef);
            AdditionHessian<ConfiguratorType> temp(factors, *D2E_bend, *D2E_membrane);
            return temp;
        }

        AdditionHessian<ConfiguratorType> produceD2E_mix(const MeshTopologySaver &plateTopol, const VectorType &plateGeomRef, const bool ActiveShellIsDeformed, const bool firstDerivWRTDef, const VectorType &edge_weights){
            auto *D2E_bend_mix = new SimpleBendingHessianMixed<ConfiguratorType>(plateTopol, plateGeomRef, true, false, edge_weights);
            auto *D2E_membrane_mix = new NonlinearMembraneHessianMixed<ConfiguratorType>(plateTopol, plateGeomRef, true, false);
            AdditionHessian<ConfiguratorType> temp(factors, *D2E_bend_mix, *D2E_membrane_mix);
            return temp;
        }
    private:
        VectorType factors;
};
