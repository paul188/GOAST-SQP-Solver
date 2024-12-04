#pragma once
#include <goast/Core.h>
#include <goast/SQP/DOFHandling/ProblemDOFs.h>

template <typename ConfiguratorType>
class MyObjectFactory {
    protected:   
        typedef typename ConfiguratorType::RealType    RealType;
        typedef typename ConfiguratorType::VectorType  VectorType;

    public:
        MyObjectFactory(const VectorType &factors, 
                        const MeshTopologySaver &plateTopol, 
                        const VectorType &edge_weights) : 
                        _factors(factors), 
                        _plateTopol(plateTopol), 
                        _edge_weights(edge_weights){}
        ~MyObjectFactory(){}

        void produceE(const ProblemDOFs<ConfiguratorType> &problemDOFs, RealType &Dest) const {
            // first, generate plateGeomRef
            VectorType plateGeomRef;
            problemDOFs.getCurrentReferenceGeometry(plateGeomRef);

            VectorType plateGeomDef = problemDOFs.getVertexDOFs();

            // then, generate the energy
            auto *E_bend = new SimpleBendingEnergy<ConfiguratorType>(_plateTopol, plateGeomRef, true, _edge_weights);
            auto *E_membrane = new NonlinearMembraneEnergy<ConfiguratorType>(_plateTopol, plateGeomRef, true);
            AdditionOp<ConfiguratorType> temp(_factors, *E_bend, *E_membrane);
            temp.apply(plateGeomDef,Dest);
        }

        void produceDE_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, VectorType& Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();

            VectorType plateGeomDef = problemDOFs.getVertexDOFs();

            // then, generate the gradient
            auto *DE_bend = new SimpleBendingGradientDef<ConfiguratorType>(_plateTopol, plateGeomRef, _edge_weights);
            auto *DE_membrane = new NonlinearMembraneGradientDef<ConfiguratorType>(_plateTopol, plateGeomRef);
            AdditionGradient<ConfiguratorType> temp(_factors, *DE_bend, *DE_membrane);
            temp.apply(plateGeomDef, Dest);
        }

        void produceD2E_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType& Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getVertexDOFs();

            // then, generate the hessian
            auto *D2E_bend = new SimpleBendingHessianDef<ConfiguratorType>(_plateTopol, plateGeomRef, _edge_weights);
            auto *D2E_membrane = new NonlinearMembraneHessianDef<ConfiguratorType>(_plateTopol, plateGeomRef);
            AdditionHessian<ConfiguratorType> temp(_factors, *D2E_bend, *D2E_membrane);
            temp.apply(plateGeomDef, Dest);
        }

        // produce the Hessian w.r.t. vertex and fold dofs
        void produceD2E_mix(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType &Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getVertexDOFs();
            
            // generate mixed hessians w.r.t. deformed and undeformed geometry
            auto *D2E_bend_def_undef = new SimpleBendingHessianMixed<ConfiguratorType>(_plateTopol, plateGeomRef, true, true, _edge_weights);
            auto *D2E_membrane_def_undef = new NonlinearMembraneHessianMixed<ConfiguratorType>(_plateTopol, plateGeomRef, true, true);
            AdditionHessian<ConfiguratorType> temp(_factors, *D2E_bend_def_undef, *D2E_membrane_def_undef);
            MatrixType D2E_def_undef_val;
            temp.apply(plateGeomDef, D2E_def_undef_val);

            // gradient of the reference geometry w.r.t. the fold DOFs
            MatrixType DFoldDofs_val = problemDOFs.getReferenceGeometryGradient();

            Dest = D2E_def_undef_val*DFoldDofs_val;
        }

    private:
        const VectorType &_factors;
        const MeshTopologySaver &_plateTopol;
        const VectorType &_edge_weights;
};
