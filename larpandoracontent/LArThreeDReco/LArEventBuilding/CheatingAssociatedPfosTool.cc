/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/CheatingAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the cheating associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/CheatingAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef NeutrinoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

CheatingAssociatedPfosTool::CheatingAssociatedPfosTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingAssociatedPfosTool::Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const Vertex *const pNeutrinoVertex, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
    
    std::cout << "\033[1;31m******************************************************************************\033[0m" << std::endl;
    
    // Make a multimap from the MC particles to all the corresponding PfoInfo objects.
    MCParticlePfoInfoMultiMap mcParticlePfoInfoMultiMap;
    MCParticleSet allMcParticles;
        
    for (const PfoInfoMap::value_type &mapEntry : pfoInfoMap)
    {
        try
        {
            const MCParticle * const pMainMCParticle(LArMCParticleHelper::GetMainMCParticle(mapEntry.first));
            allMcParticles.insert(pMainMCParticle);
            mcParticlePfoInfoMultiMap.emplace(pMainMCParticle, mapEntry.second);
        }
        catch (...)
        {
        }
    }
        
    MCParticlePfoInfoMap mcParticleToClosestPfoInfoMap;

    for (const MCParticle * const pMCParticle : allMcParticles)
    {
        const CartesianVector vtxPos(pMCParticle->GetVertex());
        
        // Find all the PFParticles corresponding to this MC particle and store them in a list.
        PfoInfoDistanceList pfoInfoDistanceList;
        
        const std::pair<MCParticlePfoInfoMultiMap::const_iterator, MCParticlePfoInfoMultiMap::const_iterator> ret = mcParticlePfoInfoMultiMap.equal_range(pMCParticle);
        if (ret.first == ret.second) continue; // if the MCParticle is not in the multimap
        
        for (MCParticlePfoInfoMultiMap::const_iterator iter = ret.first; iter != ret.second; ++iter)
        {
            PfoInfo * const pPfoInfo(iter->second);
            const float closestDistance(this->GetClosestDistance(pPfoInfo->GetThisPfo(), vtxPos));
            pfoInfoDistanceList.emplace_back(pPfoInfo, closestDistance);
        }
        
        pfoInfoDistanceList.sort([](const PfoInfoDistanceList::value_type &lhs, const PfoInfoDistanceList::value_type &rhs) { return lhs.second > rhs.second; });
        mcParticleToClosestPfoInfoMap.emplace(pMCParticle, pfoInfoDistanceList.front().first);
        
        // Make the internal p-d links.
        for (PfoInfoDistanceList::const_iterator iter = pfoInfoDistanceList.begin(); std::next(iter, 1) != pfoInfoDistanceList.end(); ++iter)
        {
            const PfoInfoDistanceList::const_iterator nextIter(std::next(iter, 1));
            
            PfoInfo * const pParentPfoInfo(iter->first);
            PfoInfo * const pDaughterPfoInfo(nextIter->first);
            
            pParentPfoInfo->AddDaughterPfo(pDaughterPfoInfo->GetThisPfo());
            pDaughterPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
            
            const LArPointingCluster pointingCluster(*(pDaughterPfoInfo->GetSlidingFitResult3D()));
            const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - vtxPos).GetMagnitudeSquared() <
                                (pointingCluster.GetOuterVertex().GetPosition() - vtxPos).GetMagnitudeSquared());
            
            pDaughterPfoInfo->SetInnerLayerAssociation(useInner);
            
            std::cout << "\033[1;31m ---> Forged internal p-d link \033[0m" << std::endl;
        }
    }
    
    // Form the vertex associations.
    for (const MCParticlePfoInfoMap::value_type &mapEntry : mcParticleToClosestPfoInfoMap)
    {
        const MCParticle * const pMCParticle(mapEntry.first);
        PfoInfo * const pPfoInfo(mapEntry.second);
                
        // Work out vertex association.        
        const CartesianVector &neutrinoVertex(pNeutrinoVertex->GetPosition());
        const LArPointingCluster pointingCluster(*(pPfoInfo->GetSlidingFitResult3D()));
        
        const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - neutrinoVertex).GetMagnitudeSquared() <
                            (pointingCluster.GetOuterVertex().GetPosition() - neutrinoVertex).GetMagnitudeSquared());
        
        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle) && LArMCParticleHelper::IsPrimary(pMCParticle))
        {
            pPfoInfo->SetNeutrinoVertexAssociation(true);
            pPfoInfo->SetInnerLayerAssociation(useInner);
            std::cout << "\033[1;31m ---> Made vtx association \033[0m" << std::endl;
        }
    }
    
    // Make the real p-d links.
    for (const MCParticle * pMCParticle : allMcParticles)
    {
        for (const MCParticle * pDaughterMCParticle : pMCParticle->GetDaughterList())
        {
            const MCParticlePfoInfoMap::const_iterator findIter(mcParticleToClosestPfoInfoMap.find(pDaughterMCParticle));
            if (findIter == mcParticleToClosestPfoInfoMap.end()) continue;
            
            PfoInfo * const pClosestDaughterPfoInfo(findIter->second);
            const CartesianVector vtxPos(pDaughterMCParticle->GetVertex());
            
            const std::pair<MCParticlePfoInfoMultiMap::const_iterator, MCParticlePfoInfoMultiMap::const_iterator> ret = mcParticlePfoInfoMultiMap.equal_range(pMCParticle);
            float bestDistance(std::numeric_limits<float>::max());
            PfoInfo * pBestParentPfoInfo(nullptr);
            
            for (MCParticlePfoInfoMultiMap::const_iterator iter = ret.first; iter != ret.second; ++iter)
            {
                const float distance(this->GetClosestDistance(iter->second->GetThisPfo(), vtxPos));
                
                if (distance < bestDistance)
                {
                    bestDistance = distance;
                    pBestParentPfoInfo = iter->second;
                }
            }
            
            if (pBestParentPfoInfo)
            {
                pBestParentPfoInfo->AddDaughterPfo(pClosestDaughterPfoInfo->GetThisPfo());
                pClosestDaughterPfoInfo->SetParentPfo(pBestParentPfoInfo->GetThisPfo());
                
                const LArPointingCluster pointingCluster(*(pClosestDaughterPfoInfo->GetSlidingFitResult3D()));
                const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - vtxPos).GetMagnitudeSquared() <
                                    (pointingCluster.GetOuterVertex().GetPosition() - vtxPos).GetMagnitudeSquared());
                
                pClosestDaughterPfoInfo->SetInnerLayerAssociation(useInner);
                std::cout << "\033[1;31m ---> Forged real p-d link \033[0m" << std::endl;
            }
        }
    }
    
    std::cout << "\033[1;31m******************************************************************************\033[0m" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CheatingAssociatedPfosTool::GetClosestDistance(const ParticleFlowObject * const pPfo, const CartesianVector vtxPos) const
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, HitType::TPC_3D, clusterList);

    if (clusterList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float bestDistance(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pPfoCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(vtxPos, pPfoCluster));

        if (thisDistance < bestDistance)
            bestDistance = thisDistance;
    }

    return bestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingAssociatedPfosTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    (void) xmlHandle;
    /*
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNeutrinoVertexDistance", m_minNeutrinoVertexDistance));
    */
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
