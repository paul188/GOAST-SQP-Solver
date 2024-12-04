#pragma once

#include <goast/Core.h>
#include "../Headers/Utils/ObjectFactory.h"
#include "../Headers/DOFHandling/FoldDofs.h"

// This is the first derivative of the energy w.r.t. the deformed geometry,
// but takes as input the foldDof Vector
// This is first and foremost for testing the derivative
// Deformed geometry is held constant
template <typename ConfiguratorType>
class EnergyFoldDof : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::VectorType> {
    protected:
        typedef typename ConfiguratorType::RealType    RealType;
        typedef typename ConfiguratorType::VectorType  VectorType;
        typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

        const MeshTopologySaver& _topology;
        const VectorType& _plateGeomDef;
        const VectorType& _edge_weights;
        const MyObjectFactory<ConfiguratorType>& _objectFactory;
        std::unique_ptr<FoldDofs<ConfiguratorType>> _foldDofsPtr;
        const BoundaryDOFS<ConfiguratorType>& _boundaryDOFs;

    public:
        EnergyFoldDof(const MeshTopologySaver &topology, 
              const VectorType &plateGeomDef, 
              const VectorType &edge_weights, 
              const MyObjectFactory<ConfiguratorType> &objectFactory, 
              std::unique_ptr<FoldDofs<ConfiguratorType>> foldDofsPtr,
              const BoundaryDOFS<ConfiguratorType>& boundaryDOFs)
                : _plateGeomDef(plateGeomDef), 
                _objectFactory(objectFactory), 
                _topology(topology), 
                _edge_weights(edge_weights), 
                _foldDofsPtr(std::move(foldDofsPtr)),
                _boundaryDOFs( boundaryDOFs ) {}

        void apply(const VectorType &t, VectorType &dest) const override{
            VectorType plateGeomRef_temp;
            _foldDofsPtr -> apply(t, plateGeomRef_temp);
            auto DE = _objectFactory.produceDE(_topology, plateGeomRef_temp,_edge_weights);
            DE.apply(_plateGeomDef, dest);
            _boundaryDOFs.transformToReducedSpace(dest);
        }
};


template <typename ConfiguratorType>
class EnergyGradientFoldDof : public BaseOp<typename ConfiguratorType::VectorType, typename ConfiguratorType::SparseMatrixType> {
    protected:
        protected:
            typedef typename ConfiguratorType::RealType    RealType;
            typedef typename ConfiguratorType::VectorType  VectorType;
            typedef typename ConfiguratorType::SparseMatrixType  MatrixType;

            const MeshTopologySaver& _topology;
            const VectorType& _plateGeomDef;
            const VectorType& _edge_weights;
            const MyObjectFactory<ConfiguratorType>& _objectFactory;
            std::unique_ptr<FoldDofs<ConfiguratorType>> _foldDofsPtr;
            std::unique_ptr<FoldDofsGradient<ConfiguratorType>> _DfoldDofsPtr;
            const BoundaryDOFS<ConfiguratorType>& _boundaryDOFs;

        public:
            EnergyGradientFoldDof(const MeshTopologySaver &topology,
                                  const VectorType &plateGeomDef, 
                                  const VectorType &edge_weights,
                                  std::unique_ptr<FoldDofs<ConfiguratorType>> foldDofsPtr, 
                                  std::unique_ptr<FoldDofsGradient<ConfiguratorType>> DFoldDofsPtr,
                                  const MyObjectFactory<ConfiguratorType>& objectFactory,
                                  const BoundaryDOFS<ConfiguratorType>& boundaryDOFs)
                                : _topology(topology),
                                  _plateGeomDef(plateGeomDef),
                                  _edge_weights(edge_weights),
                                  _foldDofsPtr(std::move(foldDofsPtr)),
                                  _DfoldDofsPtr(std::move(DFoldDofsPtr)),
                                  _objectFactory(objectFactory),
                                  _boundaryDOFs(boundaryDOFs) {}

            void apply(const VectorType &t, MatrixType &dest) const override{
                VectorType plateGeomRef_temp;
                _foldDofsPtr -> apply(t, plateGeomRef_temp);
                auto D2E_mix = _objectFactory.produceD2E_mix(_topology, plateGeomRef_temp, true, false, _edge_weights);
                MatrixType D2E_mix_val;
                D2E_mix.apply(_plateGeomDef, D2E_mix_val);

                MatrixType DFoldDofs_val;
                _DfoldDofsPtr -> apply(t,DFoldDofs_val); 
                dest = D2E_mix_val.transpose()*DFoldDofs_val;
                _boundaryDOFs.transformRowToReducedSpace(dest);
            }
};