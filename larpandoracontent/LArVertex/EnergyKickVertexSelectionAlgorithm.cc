/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the energy kick vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EnergyKickVertexSelectionAlgorithm::EnergyKickVertexSelectionAlgorithm() :
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_rOffset(10.f),
    m_xOffset(0.06),
    m_epsilon(0.06),
    m_asymmetryConstant(3.f),
    m_maxAsymmetryDistance(5.f),
    m_minAsymmetryCosAngle(0.9962),
    m_maxAsymmetryNClusters(2),
    
    m_beamDeweightingConstant(1.f),
    m_showerDeweightingConstant(1.f),
    m_showerCollapsingConstant(1.f),
    m_minShowerSpineLength(15.f),
    m_showerClusteringDistance(3.f),
    m_showerAngleConstant(0.f),
    m_showerDistanceConstant(0.f),
    m_vertexClusterDistance(4.f),
    m_tempShowerLikeStrength(0.f),
    m_globalAsymmetryConstant(1.f),
    m_useGlobalEnergyAsymmetry(false),
    m_showerAsymmetryConstant(1.f),
    m_useShowerEnergyAsymmetry(false),
    m_showerClusterNumberConstant(1.f),
    m_minShowerInwardnessDistance(1.5),
    m_showerInwardnessConstant(0.f),
    m_minShowerClusterHits(12),
    m_closestSlidingFitCanBeShowers(false),
    m_noLocalShowerAsymmetry(false),
    m_useShowerClusteringApproximation(false),
    m_useShowerClusterNumber(false)
    
    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &/*kdTreeU*/, HitKDTree2D &/*kdTreeV*/, HitKDTree2D &/*kdTreeW*/, VertexScoreList &vertexScoreList) const
{
    ClusterList clustersU, clustersV, clustersW;

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
             if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                 std::cout << "EnergyKickVertexSelectionAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
            
        ClusterList &clusterList((TPC_VIEW_U == hitType) ? clustersU : (TPC_VIEW_V == hitType) ? clustersV : clustersW);
        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }
    
    ShowerClusterMap showerClusterMapU, showerClusterMapV, showerClusterMapW;
    ShowerClusterList showerClusterListU, showerClusterListV, showerClusterListW;
    
    this->CalculateShowerClusterMap(clustersU, showerClusterListU, showerClusterMapU);
    this->CalculateShowerClusterMap(clustersV, showerClusterListV, showerClusterMapV);
    this->CalculateShowerClusterMap(clustersW, showerClusterListW, showerClusterMapW);
    
    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW, singleClusterSlidingFitDataListU, singleClusterSlidingFitDataListV, singleClusterSlidingFitDataListW;
    
    this->CalculateClusterSlidingFits(slidingFitDataListU, clustersU, showerClusterListU, singleClusterSlidingFitDataListU);
    this->CalculateClusterSlidingFits(slidingFitDataListV, clustersV, showerClusterListV, singleClusterSlidingFitDataListV);
    this->CalculateClusterSlidingFits(slidingFitDataListW, clustersW, showerClusterListW, singleClusterSlidingFitDataListW);
    
    for (const Vertex *const pVertex : vertexVector)
    {
        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float beamDeweightingScore(this->IsBeamModeOn() ? (-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant() / m_beamDeweightingConstant) : 0.f);

        float energyKick(0.f), energyAsymmetry(0.f), globalEnergyAsymmetry(0.f), showerEnergyAsymmetry(0.f);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), energyKick, energyAsymmetry, globalEnergyAsymmetry, slidingFitDataListU, showerClusterMapU, showerEnergyAsymmetry, singleClusterSlidingFitDataListU);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), energyKick, energyAsymmetry, globalEnergyAsymmetry, slidingFitDataListV, showerClusterMapV, showerEnergyAsymmetry, singleClusterSlidingFitDataListV);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), energyKick, energyAsymmetry, globalEnergyAsymmetry, slidingFitDataListW, showerClusterMapW, showerEnergyAsymmetry, singleClusterSlidingFitDataListW);

        const float energyKickScore(-energyKick / m_epsilon);
        const float energyAsymmetryScore(energyAsymmetry / m_asymmetryConstant);
        const float globalEnergyAsymmetryScore(m_useGlobalEnergyAsymmetry ? globalEnergyAsymmetry / m_globalAsymmetryConstant : 0.f);
        const float showerEnergyAsymmetryScore(m_useShowerEnergyAsymmetry ? showerEnergyAsymmetry / m_showerAsymmetryConstant : 0.f);

        vertexScoreList.emplace_back(pVertex, beamDeweightingScore + energyKickScore + energyAsymmetryScore + globalEnergyAsymmetryScore + showerEnergyAsymmetryScore);
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateShowerClusterMap(const ClusterList &inputClusterList, ShowerClusterList &showerClusterList, ShowerClusterMap &showerClusterMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    std::map<const Cluster *const, std::pair<CartesianVector, CartesianVector>> hitPositionMap;
    
    ClusterList showerLikeClusters;
    
    for (const Cluster *const pCluster : inputClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minShowerClusterHits)
            continue;
        
        if (this->IsClusterShowerLike(pCluster))
            showerLikeClusters.push_back(pCluster);
            
        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
        
        CaloHitVector clusterCaloHitVector(clusterCaloHitList.begin(), clusterCaloHitList.end());
        std::sort(clusterCaloHitVector.begin(), clusterCaloHitVector.end(), LArClusterHelper::SortHitsByPosition);
        
        if (clusterCaloHitVector.empty())
            continue;
        
        hitPositionMap.emplace(pCluster, std::make_pair(clusterCaloHitVector.front()->GetPositionVector(), clusterCaloHitVector.back()->GetPositionVector()));
    }
    
    ClusterVector showerLikeClusterVector(showerLikeClusters.begin(), showerLikeClusters.end());
    std::sort(showerLikeClusterVector.begin(), showerLikeClusterVector.end(), LArClusterHelper::SortByNHits);
        
    while (!showerLikeClusterVector.empty())
    {
        ClusterList clusterList;
        clusterList.push_back(showerLikeClusterVector.back());
        showerLikeClusterVector.pop_back();
        
        
        bool addedCluster = true;
        while (addedCluster && !showerLikeClusterVector.empty())
        {
            addedCluster = false;
            for (const Cluster * const pCluster : clusterList)
            {
                const auto &clusterVectorPairIter = hitPositionMap.find(pCluster);
                if (clusterVectorPairIter == hitPositionMap.end())
                    continue;
                    
                const auto &clusterVectorPair = clusterVectorPairIter->second;
                for (auto iter = showerLikeClusterVector.begin(); iter != showerLikeClusterVector.end(); /* no increment */)
                {
                    const auto &newVectorPairIter = hitPositionMap.find(*iter);
                    if (newVectorPairIter == hitPositionMap.end())
                        continue;
                        
                    const auto &newVectorPair = newVectorPairIter->second;
                    
                    const bool satisfiesClosestDistance = m_useShowerClusteringApproximation ? ((newVectorPair.first - clusterVectorPair.first).GetMagnitude() < m_showerClusteringDistance ||
                        (newVectorPair.first - clusterVectorPair.second).GetMagnitude() < m_showerClusteringDistance ||
                        (newVectorPair.second - clusterVectorPair.first).GetMagnitude() < m_showerClusteringDistance ||
                        (newVectorPair.second - clusterVectorPair.second).GetMagnitude() < m_showerClusteringDistance) : (LArClusterHelper::GetClosestDistance(pCluster, *iter) < m_showerClusteringDistance);
                    
                    if (satisfiesClosestDistance)
                    {
                        clusterList.push_back(*iter);
                        showerLikeClusterVector.erase(iter);
                        addedCluster = true;
                        
                        break;
                    }
                    
                    ++iter;
                }
                
                if (addedCluster)
                    break;
            }
        }
        
        ClusterVector clusterVector(clusterList.begin(), clusterList.end());
        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
        
        int totHits(0);
        for (const Cluster * const pCluster : clusterVector)
            totHits += pCluster->GetNCaloHits();
        
        if (totHits < m_minClusterCaloHits)
            continue;
        
        showerClusterList.emplace_back(clusterVector, slidingFitPitch, m_slidingFitWindow);
    }
    
    for (const ShowerCluster &showerCluster : showerClusterList)
    {
        for (const Cluster * const pCluster : showerCluster.GetClusters())
            showerClusterMap.emplace(pCluster, showerCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateClusterSlidingFits(SlidingFitDataList &slidingFitDataList, const ClusterList &inputClusterList, const ShowerClusterList &showerClusterList, SlidingFitDataList &singleClusterSlidingFitDataList) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    ClusterVector sortedClusters(inputClusterList.begin(), inputClusterList.end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster * const pCluster : sortedClusters)
    {
        if (m_closestSlidingFitCanBeShowers && this->IsClusterShowerLike(pCluster))
            continue;
        
        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;
            
        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
        slidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
    }

    for (const Cluster * const pCluster : sortedClusters)
    {
        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;
            
        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
        singleClusterSlidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
    }
    
    if (m_closestSlidingFitCanBeShowers)
    {
        for (const ShowerCluster &showerCluster : showerClusterList)
        {        
            int totHits(0);
            for (const Cluster * const pCluster : showerCluster.GetClusters())
                totHits += pCluster->GetNCaloHits();
                
            if (totHits < m_minClusterCaloHits)
                continue;
            
            // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
            const int slidingFitWindow(std::min(static_cast<int>(totHits), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
            slidingFitDataList.emplace_back(showerCluster.GetClusters(), slidingFitWindow, slidingFitPitch, m_minClusterCaloHits);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyScoresForView(const CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, float &globalEnergyAsymmetry, 
    const SlidingFitDataList &slidingFitDataList, const ShowerClusterMap &showerClusterMap, float &showerEnergyAsymmetry, const SlidingFitDataList &singleClusterSlidingFitDataList) const
{
    (void) showerEnergyAsymmetry;
    
    unsigned int totHits(0);
    bool useEnergy(true), useAsymmetry(true);
    float totEnergy(0.f), totEnergyKick(0.f), totHitKick(0.f);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);
    ClusterVector asymmetryClusters;

    // Find the closest cluster to the vertex candidate.
    SlidingFitDataList closeSlidingFitData;
    
    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        for (const Cluster *const pCluster : slidingFitData.GetClusterVector())
        {
            if (this->IsClusterShowerLike(pCluster))
                continue;
            
            const float distance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
            
            if (distance < m_vertexClusterDistance)
            {
                closeSlidingFitData.push_back(slidingFitData);
                break;
            }
        }
    }
    
    SlidingFitData closestSlidingFitData;
    float largestFitEnergy(0.f);
    
    for (const SlidingFitData &slidingFitData : closeSlidingFitData)
    {
        float totalEnergy(0.f);
        
        for (const Cluster * const pCluster : slidingFitData.GetClusterVector())
            totalEnergy += pCluster->GetElectromagneticEnergy();
            
        if (totalEnergy > largestFitEnergy)
        {
            closestSlidingFitData = slidingFitData;
            largestFitEnergy = totalEnergy;
        }
    }

    for (const SlidingFitData &slidingFitData : singleClusterSlidingFitDataList)
    {
        for (const Cluster *const pCluster : slidingFitData.GetClusterVector())
        {
            if (pCluster->GetElectromagneticEnergy() < std::numeric_limits<float>::epsilon())
                useEnergy = false;

            const CartesianVector vertexToMinLayer(slidingFitData.GetMinLayerPosition() - vertexPosition2D);
            const CartesianVector vertexToMaxLayer(slidingFitData.GetMaxLayerPosition() - vertexPosition2D);

            const bool minLayerClosest(vertexToMinLayer.GetMagnitudeSquared() < vertexToMaxLayer.GetMagnitudeSquared());
            const CartesianVector &clusterDisplacement((minLayerClosest) ? vertexToMinLayer : vertexToMaxLayer);
            const CartesianVector &clusterDirection((minLayerClosest) ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection());

            this->IncrementEnergyKickParameters(pCluster, clusterDisplacement, clusterDirection, totEnergyKick, totEnergy, totHitKick, totHits, closestSlidingFitData, showerClusterMap);

            if (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance)
            {
                if (m_noLocalShowerAsymmetry && this->IsClusterShowerLike(pCluster))
                    continue;
                    
                useAsymmetry &= this->IncrementEnergyAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);
                useAsymmetry &= this->IncrementEnergyAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);
                asymmetryClusters.push_back(pCluster);
            }
        }
    }
    
    float newShowerEnergyAsymmetry(1.f);
    for (const auto &mapEntry : showerClusterMap)
    {
        bool useShowerCluster = false;
        const auto &showerCluster = mapEntry.second;
        
        for (const Cluster * const pCluster : showerCluster.GetClusters())
        {
            if (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) > m_vertexClusterDistance)
                continue;
            
            useShowerCluster = true;
            break;
        }
        
        if (useShowerCluster)
        {
            const auto &showerFit = showerCluster.GetFit();
            
            float rL(0.f), rT(0.f);
            showerFit.GetLocalPosition(vertexPosition2D, rL, rT);
            
            CartesianVector showerDirection(0.f, 0.f, 0.f);
            showerFit.GetGlobalFitDirection(rL, showerDirection);
            
            const float projectedVtxPosition = vertexPosition2D.GetDotProduct(showerDirection);
            float beforeVtxEnergy(0.f), afterVtxEnergy(0.f);
            
            for (const Cluster * const pCluster : showerCluster.GetClusters())
            {
                CaloHitList caloHitList;
                pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

                CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
                std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

                for (const CaloHit *const pCaloHit : caloHitVector)
                {
                    if (pCaloHit->GetPositionVector().GetDotProduct(showerDirection) < projectedVtxPosition)
                        beforeVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                        
                    else if (pCaloHit->GetPositionVector().GetDotProduct(showerDirection) > projectedVtxPosition)
                        afterVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                }
            }
            
            if (beforeVtxEnergy + afterVtxEnergy > 0.f)
                newShowerEnergyAsymmetry = std::fabs(afterVtxEnergy - beforeVtxEnergy) / (afterVtxEnergy + beforeVtxEnergy);
                
            break;
        }
    }
        
    showerEnergyAsymmetry += newShowerEnergyAsymmetry;
    
    // Default: maximum asymmetry (i.e. not suppressed), zero for energy kick (i.e. not suppressed)
    if ((0 == totHits) || (useEnergy && (totEnergy < std::numeric_limits<float>::epsilon())))
    {
        energyAsymmetry += 1.f;
        return;
    }

    energyKick += useEnergy ? (totEnergyKick / totEnergy) : (totHitKick / static_cast<float>(totHits));
    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);
    energyAsymmetry += useAsymmetry ? this->CalculateEnergyAsymmetry(useEnergy, vertexPosition2D, asymmetryClusters, localWeightedDirectionSum) : 1.f;
    
    if ((useEnergy && energyWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)) || (!useEnergy && hitWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)))
        globalEnergyAsymmetry += 0.f;
    
    else
        this->IncrementGlobalEnergyAsymmetry(globalEnergyAsymmetry, useEnergy, vertexPosition2D, singleClusterSlidingFitDataList, localWeightedDirectionSum); //change me
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementGlobalEnergyAsymmetry(float &globalEnergyAsymmetry, const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const SlidingFitDataList &slidingFitDataList, const CartesianVector &localWeightedDirectionSum) const
{
    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxHitEnergy(0.f), afterVtxHitEnergy(0.f);
    unsigned int beforeVtxHitCount(0), afterVtxHitCount(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        for (const Cluster *const pCluster : slidingFitData.GetClusterVector())
        {
        
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
            std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

            for (const CaloHit *const pCaloHit : caloHitVector)
            {
                if (pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection) < evtProjectedVtxPos)
                {
                    beforeVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++beforeVtxHitCount;
                }
                else
                {
                    afterVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++afterVtxHitCount;
                }
            }
        }
    }
    
    // Use energy metrics if possible, otherwise fall back on hit counting.
    const float totHitEnergy(beforeVtxHitEnergy + afterVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);
    
    if (useEnergyMetrics)
        globalEnergyAsymmetry += std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;
    
    else
    {
        if (0 == totHitCount)
            throw StatusCodeException(STATUS_CODE_FAILURE);
            
        globalEnergyAsymmetry += std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyKickParameters(const Cluster *const pCluster, const CartesianVector &clusterDisplacement,
    const CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits, const SlidingFitData &closestSlidingFitData, const ShowerClusterMap &showerClusterMap) const
{
    float impactParameter(clusterDisplacement.GetCrossProduct(clusterDirection).GetMagnitude());   
    const float displacement(clusterDisplacement.GetMagnitude());
        
    if (this->IsClusterShowerLike(pCluster))
    {       
        const float showerCollapsingFactor = m_showerCollapsingConstant * (1.f - m_tempShowerLikeStrength * this->IsClusterShowerLike(pCluster, closestSlidingFitData, showerClusterMap));
 
    //    std::cout << "SHOWER COLLAPSING FACTOR: " << showerCollapsingFactor << std::endl;
 
        totEnergyKick += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy() * 
                        (showerCollapsingFactor * impactParameter + m_xOffset) / (displacement + m_rOffset);
                        
        totEnergy += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy();
        
        totHitKick += static_cast<float>(m_showerDeweightingConstant) * static_cast<float>(pCluster->GetNCaloHits()) * 
                      (showerCollapsingFactor * impactParameter + m_xOffset) / (displacement + m_rOffset);
                      
        totHits += static_cast<int>(std::round(m_showerDeweightingConstant * static_cast<float>(pCluster->GetNCaloHits())));
    }

    else
    {
        totEnergyKick += pCluster->GetElectromagneticEnergy() * (impactParameter + m_xOffset) / (displacement + m_rOffset);
        totEnergy += pCluster->GetElectromagneticEnergy();
        
        totHitKick += static_cast<float>(pCluster->GetNCaloHits()) * (impactParameter + m_xOffset) / (displacement + m_rOffset);
        totHits += pCluster->GetNCaloHits();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster, const SlidingFitData &closestSlidingFitData, const ShowerClusterMap &showerClusterMap) const
{
    if (pCluster->GetParticleId() != E_MINUS || LArClusterHelper::GetLength(pCluster) >= m_minShowerSpineLength)
        return 1.f;
    
    auto findIter = showerClusterMap.find(pCluster);
    if (findIter == showerClusterMap.end()) // it's a tiny bit of shower that was too crap to incorporate into a shower, so we definitely want to collapse it
        return 1.f;
     
    if (closestSlidingFitData.GetClusterVector().empty())
        return 1.f;
        
    const ShowerCluster &showerCluster = findIter->second;
    
    auto showerDistanceMapIter = m_showerDataMap.find(std::make_pair(&showerCluster, &closestSlidingFitData));
    if (showerDistanceMapIter != m_showerDataMap.end())
        return showerDistanceMapIter->second;
    
    float closestShowerDistance = std::numeric_limits<float>::max();
    
    for (const Cluster * const pClosestCluster : closestSlidingFitData.GetClusterVector())
    {
        ClusterList showerClusterList;
        
        for (const Cluster *const pShowerCluster : showerCluster.GetClusters())
            showerClusterList.push_back(pShowerCluster);
        
        const float distance = LArClusterHelper::GetClosestDistance(pClosestCluster, showerClusterList);
        
        if (distance < closestShowerDistance)
            closestShowerDistance = distance;
    }
                             
    CartesianVector closestShowerFitDirection(0.f, 0.f, 0.f), closestFitDirection(0.f, 0.f, 0.f), closestFitPosition(0.f, 0.f, 0.f), farthestFitPosition(0.f, 0.f, 0.f);
    const float maxMaxDistance = (showerCluster.GetFit().GetGlobalMaxLayerPosition() - closestSlidingFitData.GetMaxLayerPosition()).GetMagnitude();
    const float maxMinDistance = (showerCluster.GetFit().GetGlobalMaxLayerPosition() - closestSlidingFitData.GetMinLayerPosition()).GetMagnitude();
    const float minMaxDistance = (showerCluster.GetFit().GetGlobalMinLayerPosition() - closestSlidingFitData.GetMaxLayerPosition()).GetMagnitude();
    const float minMinDistance = (showerCluster.GetFit().GetGlobalMinLayerPosition() - closestSlidingFitData.GetMinLayerPosition()).GetMagnitude();
                            
    if (maxMaxDistance <= maxMinDistance && maxMaxDistance <= minMaxDistance && maxMaxDistance <= minMinDistance)
    {
        closestShowerFitDirection = showerCluster.GetFit().GetGlobalMaxLayerDirection();
        closestFitDirection = closestSlidingFitData.GetMaxLayerDirection();
        
        closestFitPosition = closestSlidingFitData.GetMaxLayerPosition();
        farthestFitPosition = closestSlidingFitData.GetMinLayerPosition();
    }
        
    else if (maxMinDistance <= maxMaxDistance && maxMinDistance <= minMaxDistance && maxMinDistance <= minMinDistance)
    {
        closestShowerFitDirection = showerCluster.GetFit().GetGlobalMaxLayerDirection();
        closestFitDirection = closestSlidingFitData.GetMinLayerDirection();
        
        closestFitPosition = closestSlidingFitData.GetMinLayerPosition();
        farthestFitPosition = closestSlidingFitData.GetMaxLayerPosition();
    }
    
    else if (minMaxDistance <= maxMaxDistance && minMaxDistance <= maxMinDistance && minMaxDistance <= minMinDistance)
    {
        closestShowerFitDirection = showerCluster.GetFit().GetGlobalMinLayerDirection();
        closestFitDirection = closestSlidingFitData.GetMaxLayerDirection();
        
        closestFitPosition = closestSlidingFitData.GetMaxLayerPosition();
        farthestFitPosition = closestSlidingFitData.GetMinLayerPosition();
    }
    
    else
    {
        closestShowerFitDirection = showerCluster.GetFit().GetGlobalMinLayerDirection();
        closestFitDirection = closestSlidingFitData.GetMinLayerDirection();
        
        closestFitPosition = closestSlidingFitData.GetMinLayerPosition();
        farthestFitPosition = closestSlidingFitData.GetMaxLayerPosition();
    }

    const float cosTheta = std::fabs(closestShowerFitDirection.GetDotProduct(closestFitDirection));
    
    
    int outsideHitCount(0), insideHitCount(0);
    
    if (maxMaxDistance < m_minShowerInwardnessDistance || maxMinDistance < m_minShowerInwardnessDistance || minMaxDistance < m_minShowerInwardnessDistance || minMinDistance < m_minShowerInwardnessDistance)
    {
        const float fitLength = (farthestFitPosition - closestFitPosition).GetDotProduct(closestFitDirection);
        
        for (const Cluster * const pShowerCluster : showerCluster.GetClusters())
        {
            CaloHitList caloHitList;
            pShowerCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList); 
            
            CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
            std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);
            
            for (const CaloHit * const pCaloHit : caloHitVector)
            {
                const float caloHitProjPos = (pCaloHit->GetPositionVector() - closestFitPosition).GetDotProduct(closestFitDirection);
                if (caloHitProjPos < 0.f || caloHitProjPos > fitLength)
                    ++outsideHitCount;
                    
                else
                    ++insideHitCount;
            }
        }
    }
    
    const float showerInwardness = ((insideHitCount + outsideHitCount) == 0) ? 0.5 : std::fabs(static_cast<float>(insideHitCount - outsideHitCount)) / static_cast<float>(insideHitCount + outsideHitCount);
    const float showerClusterNumberFactor = m_useShowerClusterNumber ? (1.f - std::exp(-m_showerClusterNumberConstant * static_cast<float>(showerCluster.GetClusters().size()))) : 1.f;
    
    const float metric = std::exp(-closestShowerDistance * m_showerDistanceConstant - m_showerAngleConstant * cosTheta) * showerClusterNumberFactor * std::exp(-m_showerInwardnessConstant * showerInwardness);
    m_showerDataMap.emplace(std::make_pair(&showerCluster, &closestSlidingFitData), metric);
    
    return metric;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnergyKickVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    auto findIter = m_showerLikeClusterMap.find(pCluster);
    if (findIter != m_showerLikeClusterMap.end())
        return findIter->second;
    
    bool isClusterShowerLike = (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
    m_showerLikeClusterMap.emplace(pCluster, isClusterShowerLike);
    
    return isClusterShowerLike;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnergyKickVertexSelectionAlgorithm::IncrementEnergyAsymmetryParameters(const float weight, const CartesianVector &clusterDirection,
    CartesianVector &localWeightedDirectionSum) const
{
    // If the new axis direction is at an angle of greater than 90 deg to the current axis direction, flip it 180 degs.
    CartesianVector newDirection(clusterDirection);

    if (localWeightedDirectionSum.GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
    {
        const float cosOpeningAngle(localWeightedDirectionSum.GetCosOpeningAngle(clusterDirection));

        if (std::fabs(cosOpeningAngle) > m_minAsymmetryCosAngle)
        {
            if (cosOpeningAngle < 0.f)
                newDirection *= -1.f;
        }
        else
        {
            return false;
        }
    }

    localWeightedDirectionSum += newDirection * weight;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::CalculateEnergyAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const ClusterVector &asymmetryClusters, const CartesianVector &localWeightedDirectionSum) const
{
    if (asymmetryClusters.empty() || (asymmetryClusters.size() > m_maxAsymmetryNClusters))
        return 1.f;

    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxHitEnergy(0.f), afterVtxHitEnergy(0.f);
    unsigned int beforeVtxHitCount(0), afterVtxHitCount(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    for (const Cluster *const pCluster : asymmetryClusters)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            if (pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection) < evtProjectedVtxPos)
            {
                beforeVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                ++beforeVtxHitCount;
            }
            else
            {
                afterVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                ++afterVtxHitCount;
            }
        }
    }

    // Use energy metrics if possible, otherwise fall back on hit counting.
    const float totHitEnergy(beforeVtxHitEnergy + beforeVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);
    
    if (useEnergyMetrics && (totHitEnergy > std::numeric_limits<float>::max())) // ATTN!
        return std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::SlidingFitData::SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, 
        const float slidingFitPitch) :
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f)
{
    m_clusterVector.push_back(pCluster);
    
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
    m_minLayerDirection = slidingFitResult.GetGlobalMinLayerDirection();
    m_maxLayerDirection = slidingFitResult.GetGlobalMaxLayerDirection();
    m_minLayerPosition = slidingFitResult.GetGlobalMinLayerPosition();
    m_maxLayerPosition = slidingFitResult.GetGlobalMaxLayerPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::SlidingFitData::SlidingFitData(const pandora::ClusterVector &clusterVector, const int slidingFitWindow, 
        const float slidingFitPitch, int minClusterCaloHits) :
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f)
{
    const TwoDSlidingFitResult slidingFitResult(clusterVector, slidingFitWindow, slidingFitPitch);
    m_minLayerDirection = slidingFitResult.GetGlobalMinLayerDirection();
    m_maxLayerDirection = slidingFitResult.GetGlobalMaxLayerDirection();
    m_minLayerPosition = slidingFitResult.GetGlobalMinLayerPosition();
    m_maxLayerPosition = slidingFitResult.GetGlobalMaxLayerPosition();
    
    for (const Cluster * const pCluster : clusterVector)
    {
        if (pCluster->GetNCaloHits() < minClusterCaloHits)
                continue;
                
        m_clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::ShowerCluster::ShowerCluster(const pandora::ClusterVector &clusterVector, const int slidingFitWindow, 
        const float slidingFitPitch) :
    m_clusterVector(clusterVector),
    m_twoDSlidingFitResult(clusterVector, slidingFitWindow, slidingFitPitch)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyKickVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ROffset", m_rOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "XOffset", m_xOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Epsilon", m_epsilon));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AsymmetryConstant", m_asymmetryConstant));

    if ((m_rOffset < std::numeric_limits<float>::epsilon()) || (m_epsilon < std::numeric_limits<float>::epsilon()) || (m_asymmetryConstant < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "EnergyKickVertexSelection: Invalid parameter(s), ROffset " << m_rOffset << ", Epsilon " << m_epsilon << ", AsymmetryConstant " << m_asymmetryConstant << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAsymmetryCosAngle", m_minAsymmetryCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryNClusters", m_maxAsymmetryNClusters));
    
     PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamDeweightingConstant", m_beamDeweightingConstant));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerDeweightingConstant", m_showerDeweightingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerCollapsingConstant", m_showerCollapsingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAngleConstant", m_showerAngleConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerDistanceConstant", m_showerDistanceConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexClusterDistance", m_vertexClusterDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TempShowerLikeStrength", m_tempShowerLikeStrength));
       
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GlobalAsymmetryConstant", m_globalAsymmetryConstant));
        
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseGlobalEnergyAsymmetry", m_useGlobalEnergyAsymmetry));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAsymmetryConstant", m_showerAsymmetryConstant));
        
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerEnergyAsymmetry", m_useShowerEnergyAsymmetry));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusterNumberConstant", m_showerClusterNumberConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerInwardnessDistance", m_minShowerInwardnessDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerInwardnessConstant", m_showerInwardnessConstant));
      
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerClusterHits", m_minShowerClusterHits));
     
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClosestSlidingFitCanBeShowers", m_closestSlidingFitCanBeShowers));   
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NoLocalShowerAsymmetry", m_noLocalShowerAsymmetry));   
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusteringApproximation", m_useShowerClusteringApproximation));   
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusterNumber", m_useShowerClusterNumber));  
        
    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
