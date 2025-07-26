#pragma once
#include <goast/Core.h>
#include <goast/SQP/DOFHandling/ProblemDOFs.h>
#include <goast/DiscreteShells.h>

template <typename ConfiguratorType>
class EnergyFactory
{
    protected:
        typedef typename ConfiguratorType::RealType    RealType;
        typedef typename ConfiguratorType::VectorType  VectorType;
    public:
        EnergyFactory(const VectorType &factors, 
                      const MeshTopologySaver &plateTopol,
                      const VectorType &edge_weights):
                        _factors(factors),
                        _plateTopol(plateTopol),
                        _edge_weights(edge_weights)
                    {
                        assert(edge_weights.size() == plateTopol.getNumEdges());
                    }
        ~EnergyFactory(){}

        virtual void produceE(const ProblemDOFs<ConfiguratorType> &problemDOFs, RealType &Dest) const = 0;
        virtual void produceDE_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, VectorType& Dest) const = 0;
        virtual void produceD2E_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType& Dest) const = 0;
        virtual void produceD2E_mix(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType &Dest) const = 0;

    protected:
        const VectorType &_factors;
        const MeshTopologySaver &_plateTopol;
        const VectorType &_edge_weights;
};

template <typename ConfiguratorType>
class ElasticEnergyFactory : public EnergyFactory<ConfiguratorType>{
    protected:   
        typedef typename ConfiguratorType::RealType    RealType;
        typedef typename ConfiguratorType::VectorType  VectorType;

    public:
        ElasticEnergyFactory(const VectorType &factors, 
                        const MeshTopologySaver &plateTopol, 
                        const VectorType &edge_weights) : 
                        EnergyFactory<ConfiguratorType>(factors, plateTopol, edge_weights){
                            assert(factors.size() == 2);
                        }
        ~ElasticEnergyFactory(){}

        void produceE(const ProblemDOFs<ConfiguratorType> &problemDOFs, RealType &Dest) const override{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // then, generate the energy
            auto *E_bend = new SimpleBendingEnergy<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, this->_edge_weights);
            auto *E_membrane = new NonlinearMembraneEnergy<ConfiguratorType>(this->_plateTopol, plateGeomRef, true);
            AdditionOp<ConfiguratorType> temp(this->_factors, *E_membrane,*E_bend);
            temp.apply(plateGeomDef,Dest);
        }

        void produceDE_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, VectorType& Dest) const override{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // then, generate the gradient
            auto *DE_bend = new SimpleBendingGradientDef<ConfiguratorType>(this->_plateTopol, plateGeomRef, this->_edge_weights);
            auto *DE_membrane = new NonlinearMembraneGradientDef<ConfiguratorType>(this->_plateTopol, plateGeomRef);
            AdditionGradient<ConfiguratorType> temp(this->_factors, *DE_membrane, *DE_bend);
            temp.apply(plateGeomDef, Dest);
        }

        void produceD2E_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType& Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // then, generate the hessian
            auto *D2E_bend = new SimpleBendingHessianDef<ConfiguratorType>(this->_plateTopol, plateGeomRef, this->_edge_weights);
            auto *D2E_membrane = new NonlinearMembraneHessianDef<ConfiguratorType>(this->_plateTopol, plateGeomRef);
            AdditionHessian<ConfiguratorType> temp(this->_factors,*D2E_membrane, *D2E_bend);
            temp.apply(plateGeomDef, Dest);
        }

        // produce the Hessian w.r.t. vertex and fold dofs
        void produceD2E_mix(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType &Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();
            
            // generate mixed hessians w.r.t. deformed and undeformed geometry
            auto *D2E_bend_def_undef = new SimpleBendingHessianMixed<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, true, this->_edge_weights);
            auto *D2E_membrane_def_undef = new NonlinearMembraneHessianMixed<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, true);
            AdditionHessian<ConfiguratorType> temp(this->_factors, *D2E_membrane_def_undef, *D2E_bend_def_undef);
            MatrixType D2E_def_undef_val;
            temp.apply(plateGeomDef, D2E_def_undef_val);

            // gradient of the reference geometry w.r.t. the fold DOFs
            MatrixType DFoldDofs_val = problemDOFs.getReferenceGeometryGradient();

            Dest = D2E_def_undef_val*DFoldDofs_val;
        }
};

// used in the arc fold optimization
/*
    Same as elastic energy factory, but with additional gravitational force acting on force_vertices
*/
template <typename ConfiguratorType>
class ElasticGravitationalFactory : public EnergyFactory<ConfiguratorType>{
    protected:   
        typedef typename ConfiguratorType::RealType    RealType;
        typedef typename ConfiguratorType::VectorType  VectorType;

    public:
        ElasticGravitationalFactory(const VectorType &factors, 
                        const MeshTopologySaver &plateTopol, 
                        const VectorType &edge_weights,
                        const VectorType &gravity_direction,
                        const VectorType &mass_distribution) :
                        EnergyFactory<ConfiguratorType>(factors,plateTopol,edge_weights),
                        _gravity_direction(gravity_direction),
                        _mass_distribution(mass_distribution)
                        {
                            assert(gravity_direction.size() == 3);
                            assert(mass_distribution.size() == plateTopol.getNumVertices());
                            assert(factors.size() == 3);
                        }
        ~ElasticGravitationalFactory(){}

        void produceE(const ProblemDOFs<ConfiguratorType> &problemDOFs, RealType &Dest) const {
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // then, generate the energy
            auto *E_bend = new SimpleBendingEnergy<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, this->_edge_weights);
            auto *E_membrane = new NonlinearMembraneEnergy<ConfiguratorType>(this->_plateTopol, plateGeomRef, true);
            auto *E_grav = new GravitationalEnergy<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, _mass_distribution, _gravity_direction);
            AdditionOp<ConfiguratorType> temp(this->_factors, *E_membrane,*E_bend, *E_grav);
            temp.apply(plateGeomDef,Dest);
        }

        void produceDE_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, VectorType& Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();

            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // then, generate the gradient
            auto *DE_bend = new SimpleBendingGradientDef<ConfiguratorType>(this->_plateTopol, plateGeomRef, this->_edge_weights);
            auto *DE_membrane = new NonlinearMembraneGradientDef<ConfiguratorType>(this->_plateTopol, plateGeomRef);
            auto *DE_grav = new GravitationalEnergyGradientDef<ConfiguratorType>(this->_plateTopol, plateGeomRef, _mass_distribution, _gravity_direction);
            AdditionGradient<ConfiguratorType> temp(this->_factors, *DE_membrane, *DE_bend, *DE_grav);
            temp.apply(plateGeomDef, Dest);
        }

        void produceD2E_vertex(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType& Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();

            // second derivative of gravitational energy is zero
            auto *D2E_bend = new SimpleBendingHessianDef<ConfiguratorType>(this->_plateTopol, plateGeomRef, this->_edge_weights);
            auto *D2E_membrane = new NonlinearMembraneHessianDef<ConfiguratorType>(this->_plateTopol, plateGeomRef);
            VectorType first_factors(2);
            first_factors[0] = this->_factors[0];
            first_factors[1] = this->_factors[1];
            AdditionHessian<ConfiguratorType> temp(first_factors,*D2E_membrane, *D2E_bend);
            temp.apply(plateGeomDef, Dest);
        }

        // produce the Hessian w.r.t. vertex and fold dofs
        void produceD2E_mix(const ProblemDOFs<ConfiguratorType> &problemDOFs, MatrixType &Dest) const{
            // first, generate plateGeomRef
            VectorType plateGeomRef = problemDOFs.getReferenceGeometry();
            VectorType plateGeomDef = problemDOFs.getDeformedGeometry();
            
            // generate mixed hessians w.r.t. deformed and undeformed geometry
            auto *D2E_bend_def_undef = new SimpleBendingHessianMixed<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, true, this->_edge_weights);
            auto *D2E_membrane_def_undef = new NonlinearMembraneHessianMixed<ConfiguratorType>(this->_plateTopol, plateGeomRef, true, true);
            VectorType first_factors(2);
            first_factors[0] = this->_factors[0];
            first_factors[1] = this->_factors[1];
            AdditionHessian<ConfiguratorType> temp(first_factors, *D2E_membrane_def_undef, *D2E_bend_def_undef);
            MatrixType D2E_def_undef_val;
            temp.apply(plateGeomDef, D2E_def_undef_val);

            // gradient of the reference geometry w.r.t. the fold DOFs
            MatrixType DFoldDofs_val = problemDOFs.getReferenceGeometryGradient();

            Dest = D2E_def_undef_val*DFoldDofs_val;
        }

    private:
        const VectorType &_mass_distribution;
        const VectorType &_gravity_direction;
};