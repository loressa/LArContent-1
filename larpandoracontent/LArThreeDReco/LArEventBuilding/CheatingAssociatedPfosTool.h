/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/CheatingAssociatedPfosTool.h
 * 
 *  @brief  Header file for cheating associated pfos tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_ASSOCIATED_PFOS_TOOL_H
#define LAR_CHEATING_ASSOCIATED_PFOS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

#include "Objects/MCParticle.h"

namespace lar_content
{

/**
 *  @brief  CheatingAssociatedPfosTool class
 */
class CheatingAssociatedPfosTool : public PfoRelationTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    CheatingAssociatedPfosTool();

    void Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    /**
     *  @brief  ...
     * 
     *  @param pPfo ...
     *
     *  @param vtxPos ...
     *
     *  @return ...
     */
    float GetClosestDistance(const pandora::Pfo * const pPfo, const pandora::CartesianVector vtxPos) const;

    typedef std::multimap<const pandora::MCParticle * const, NeutrinoHierarchyAlgorithm::PfoInfo * const> MCParticlePfoInfoMultiMap; ///< Multimap from each MCParticle to any corresponding Pfos
    typedef std::list<std::pair<NeutrinoHierarchyAlgorithm::PfoInfo * const, float>>                      PfoInfoDistanceList;       ///< ...
    typedef std::map<const pandora::MCParticle * const, NeutrinoHierarchyAlgorithm::PfoInfo * const>      MCParticlePfoInfoMap;      ///< ...

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *CheatingAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new CheatingAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_ASSOCIATED_PFOS_TOOL_H
