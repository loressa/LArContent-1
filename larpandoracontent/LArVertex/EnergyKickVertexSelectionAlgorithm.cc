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
    m_localAsymmetryConstant(0.4),
    m_minLocalAsymmetryClusterHits(10),
    m_useLocalAsymmetry(true),
    m_beamDeweightingConstant(1.f),
    m_showerDeweightingConstant(1.f),
    m_showerCollapsingConstant(0.f),
    m_minShowerSpineLength(15.f),
    m_useClusterCharacterisation(true),
    m_maxTrackLikeDeviation(2.f),
    m_showerClusteringDistance(3.f),
    m_showerAngleConstant(1.f),
    m_showerDistanceConstant(1.f),
    m_vertexClusterDistance(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &/*kdTreeU*/, HitKDTree2D &/*kdTreeV*/, HitKDTree2D &/*kdTreeW*/, VertexScoreList &vertexScoreList) const
{    
    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(slidingFitDataListU, slidingFitDataListV, slidingFitDataListW);

    ShowerClusterMap showerClusterMapU = this->CalculateShowerClusterMap(slidingFitDataListU);
    ShowerClusterMap showerClusterMapV = this->CalculateShowerClusterMap(slidingFitDataListV);
    ShowerClusterMap showerClusterMapW = this->CalculateShowerClusterMap(slidingFitDataListW);

    for (const Vertex *const pVertex : vertexVector)
    {
        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float beamDeweightingScore(this->IsBeamModeOn() ? std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant() / m_beamDeweightingConstant) : 1.f);

        float energyKick(0.f), energyAsymmetry(0.f), localEnergyAsymmetry(0.f);
        
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), energyKick, energyAsymmetry, localEnergyAsymmetry, slidingFitDataListU, showerClusterMapU, TPC_VIEW_U);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), energyKick, energyAsymmetry, localEnergyAsymmetry, slidingFitDataListV, showerClusterMapV, TPC_VIEW_V);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), energyKick, energyAsymmetry, localEnergyAsymmetry, slidingFitDataListW, showerClusterMapW, TPC_VIEW_W);

        const float energyKickScore(std::exp(-energyKick / m_epsilon));
        const float energyAsymmetryScore(std::exp(energyAsymmetry / m_asymmetryConstant));
        const float localEnergyAsymmetryScore(m_useLocalAsymmetry ? std::exp(localEnergyAsymmetry / m_localAsymmetryConstant) : 1.f);
        
        vertexScoreList.push_back(VertexScore(pVertex, beamDeweightingScore * energyKickScore * energyAsymmetryScore * localEnergyAsymmetryScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::ShowerClusterMap EnergyKickVertexSelectionAlgorithm::CalculateShowerClusterMap(const SlidingFitDataList &slidingFitDataList) const
{
    ClusterList showerLikeClusters;
    
    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());
        
        if (this->IsClusterShowerLike(pCluster))
            showerLikeClusters.push_back(pCluster);
    }
    
    ShowerClusterList showerClusterList;
    
    while (!showerLikeClusters.empty())
    {
        ClusterList clusterList;
        clusterList.push_back(showerLikeClusters.back());
        showerLikeClusters.pop_back();
        
        bool addedCluster = true;
        while (addedCluster && !showerLikeClusters.empty())
        {
            addedCluster = false;
            for (const Cluster * const pCluster : clusterList)
            {
                for (auto iter = showerLikeClusters.begin(); iter != showerLikeClusters.end(); /* no increment */)
                {
                    if (LArClusterHelper::GetClosestDistance(pCluster, *iter) < m_showerClusteringDistance)
                    {
                        clusterList.push_back(*iter);
                        showerLikeClusters.erase(iter);
                        addedCluster = true;
                        break;
                    }
                    
                    ++iter;
                }
                
                if (addedCluster)
                    break;
            }
        }
        
        showerClusterList.emplace_back(clusterList);
    }
    
    std::cout << "********** NUMBER OF SHOWERS: " << showerClusterList.size() << std::endl;
    
    //PandoraMonitoringApi::Pause(this->GetPandora());
    
    ShowerClusterMap showerClusterMap;
    for (const ShowerCluster &showerCluster : showerClusterList)
    {
        for (const Cluster * const pCluster : showerCluster.GetClusters())
            showerClusterMap.emplace(pCluster, showerCluster);
    }
    
    std::cout << "DONE" << std::endl;
    return showerClusterMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateClusterSlidingFits(SlidingFitDataList &slidingFitDataListU, SlidingFitDataList &slidingFitDataListV,
    SlidingFitDataList &slidingFitDataListW) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

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

        SlidingFitDataList &slidingFitDataList((TPC_VIEW_U == hitType) ? slidingFitDataListU : (TPC_VIEW_V == hitType) ? slidingFitDataListV : slidingFitDataListW);

        if (!slidingFitDataListW.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster * const pCluster : sortedClusters)
        {
            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;
            
            // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
            const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
            slidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyScoresForView(const CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, float &localEnergyAsymmetry,
    const SlidingFitDataList &slidingFitDataList, const ShowerClusterMap &showerClusterMap, HitType hitType) const
{
    //const CartesianVector vtxPosProj(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertexPosition2D, TPC_VIEW_W));
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vtxPosProj, "vertex candidate", GREEN, 1);
    
    unsigned int totHits(0);
    bool useEnergy(true);
    float totEnergy(0.f), totEnergyKick(0.f), totHitKick(0.f);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);
    ClusterList localAsymmetryClusterList;

    // Find the closest cluster to the vertex candidate.
    SlidingFitDataList closeSlidingFitData;
    
    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());
        const float distance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
        
        if (distance < m_vertexClusterDistance)
            closeSlidingFitData.push_back(slidingFitData);
    }

    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());

        if (pCluster->GetElectromagneticEnergy() < std::numeric_limits<float>::epsilon())
            useEnergy = false;

        const CartesianVector vertexToMinLayer(slidingFitData.GetMinLayerPosition() - vertexPosition2D);
        const CartesianVector vertexToMaxLayer(slidingFitData.GetMaxLayerPosition() - vertexPosition2D);

        const bool minLayerClosest(vertexToMinLayer.GetMagnitudeSquared() < vertexToMaxLayer.GetMagnitudeSquared());
        const CartesianVector &clusterDisplacement((minLayerClosest) ? vertexToMinLayer : vertexToMaxLayer);
        const CartesianVector &clusterDirection((minLayerClosest) ? slidingFitData.GetMinLayerDirection() : slidingFitData.GetMaxLayerDirection());

        this->IncrementEnergyKickParameters(pCluster, clusterDisplacement, clusterDirection, totEnergyKick, totEnergy, totHitKick, totHits, closeSlidingFitData, showerClusterMap, hitType);

        // Update the local event axis.
        if (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance)
        {
            this->IncrementEnergyAsymmetryParameters(pCluster->GetElectromagneticEnergy(), clusterDirection, energyWeightedDirectionSum);
            this->IncrementEnergyAsymmetryParameters(static_cast<float>(pCluster->GetNCaloHits()), clusterDirection, hitWeightedDirectionSum);
            localAsymmetryClusterList.push_back(pCluster);
        }
    }
    
    // ATTN: Look into the logic here...
    if ((0 == totHits) || (useEnergy && (totEnergy < std::numeric_limits<float>::epsilon())))
        energyKick += 0.f;
        
    else
        energyKick += useEnergy ? (totEnergyKick / totEnergy) : (totHitKick / static_cast<float>(totHits));
        
    const CartesianVector &localWeightedDirectionSum(useEnergy ? energyWeightedDirectionSum : hitWeightedDirectionSum);
        
    // Default: maximum asymmetry (i.e. not suppressed), zero for energy kick (i.e. not suppressed)
    if ((useEnergy && energyWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)) || (!useEnergy && hitWeightedDirectionSum == CartesianVector(0.f, 0.f, 0.f)))
        energyAsymmetry += 0.f;
    
    else
        this->IncrementEnergyAsymmetry(energyAsymmetry, localEnergyAsymmetry, localAsymmetryClusterList, useEnergy, vertexPosition2D, slidingFitDataList, localWeightedDirectionSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyKickParameters(const Cluster *const pCluster, const CartesianVector &clusterDisplacement,
    const CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits, const SlidingFitDataList &closeSlidingFitData, const ShowerClusterMap &showerClusterMap, HitType hitType) const
{
    float impactParameter(clusterDisplacement.GetCrossProduct(clusterDirection).GetMagnitude());   
    const float displacement(clusterDisplacement.GetMagnitude());
        
    if (this->IsClusterShowerLike(pCluster))
    {       
        const float showerCollapsingFactor = 1.f - this->IsClusterShowerLike(pCluster, closeSlidingFitData, showerClusterMap, hitType);
 
        std::cout << "SHOWER COLLAPSING FACTOR: " << showerCollapsingFactor << std::endl;
 
        totEnergyKick += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy() * 
                        (showerCollapsingFactor * impactParameter + m_xOffset) / (displacement + m_rOffset);
                        
        totEnergy += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy();
        
        totHitKick += static_cast<float>(m_showerDeweightingConstant * pCluster->GetNCaloHits()) * 
                      (showerCollapsingFactor * impactParameter + m_xOffset) / (displacement + m_rOffset);
                      
        totHits += static_cast<int>(std::round(m_showerDeweightingConstant * pCluster->GetNCaloHits()));
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

float EnergyKickVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster, const SlidingFitDataList &closeSlidingFitData, const ShowerClusterMap &showerClusterMap, HitType hitType) const
{
    if (pCluster->GetParticleId() != E_MINUS || LArClusterHelper::GetLength(pCluster) >= m_minShowerSpineLength)
        return 0.f;
    
    auto findIter = showerClusterMap.find(pCluster);
    if (findIter == showerClusterMap.end())
        throw;
        
    if (closeSlidingFitData.empty())
        return 0.f;
        
    ShowerCluster showerCluster = findIter->second;
    
    unsigned int hitCount(0);
    CartesianVector avgShowerPosition(0.f, 0.f, 0.f);
    
    float closestShowerDistance = std::numeric_limits<float>::max();
    SlidingFitData closestSlidingFitData;
    
    for (const Cluster *const pShowerCluster : showerCluster.GetClusters())
    {
        // Loop over Hits.
        CaloHitList caloHitList;
        pShowerCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            avgShowerPosition += pCaloHit->GetPositionVector();
            ++hitCount;
            
            for (const SlidingFitData &slidingFitData : closeSlidingFitData)
            {
                if (slidingFitData.GetCluster() == pShowerCluster)
                    continue;
                
                const float distance = LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), slidingFitData.GetCluster());
                if (distance < closestShowerDistance)
                {
                    closestShowerDistance = distance;
                    closestSlidingFitData = slidingFitData;
                }
            }
        }
    }
    
    if (closestShowerDistance == std::numeric_limits<float>::max())
        return 0.f;
    
    avgShowerPosition = avgShowerPosition * (1.f/static_cast<float>(hitCount));
    
    const CartesianVector closestFitEnd = (avgShowerPosition - closestSlidingFitData.GetMaxLayerPosition()).GetMagnitude() < (avgShowerPosition - closestSlidingFitData.GetMinLayerPosition()).GetMagnitude() 
                                          ? closestSlidingFitData.GetMaxLayerPosition() : closestSlidingFitData.GetMinLayerPosition();
                                          
    const CartesianVector closestFitDirection = (avgShowerPosition - closestSlidingFitData.GetMaxLayerPosition()).GetMagnitude() < (avgShowerPosition - closestSlidingFitData.GetMinLayerPosition()).GetMagnitude() 
                                          ? closestSlidingFitData.GetMaxLayerDirection() : closestSlidingFitData.GetMinLayerDirection();
    
    const float cosTheta = std::fabs((avgShowerPosition - closestFitEnd).GetDotProduct(closestFitDirection)) / ((avgShowerPosition - closestFitEnd).GetMagnitude());
    
    std::cout << "THETA IS " << std::acos(cosTheta) * 180.f / 3.141 << std::endl;
    std::cout << "CLOSEST DISTANCE IS " << closestShowerDistance << std::endl;
    
    const float metric = m_showerCollapsingConstant * std::pow(cosTheta, m_showerAngleConstant) * std::exp(-closestShowerDistance * m_showerDistanceConstant);
    
    (void) hitType;
    /*
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
    const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
    
    const CartesianVector showerStartU(LArGeometryHelper::ProjectPosition(this->GetPandora(), slidingFitResult.GetGlobalMinLayerPosition(), TPC_VIEW_W));
    const CartesianVector showerEndU(LArGeometryHelper::ProjectPosition(this->GetPandora(), slidingFitResult.GetGlobalMaxLayerPosition(), TPC_VIEW_W));
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStartU, "shower candidate", MAGENTA, 1);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerEndU, "shower candidate", MAGENTA, 1);
    const CartesianVector vtxClusterStartU(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestSlidingFitData.GetMinLayerPosition(), TPC_VIEW_W));
    const CartesianVector vtxClusterEndU(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestSlidingFitData.GetMaxLayerPosition(), TPC_VIEW_W));
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vtxClusterStartU, "vtx cluster", RED, 1);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vtxClusterEndU, "vtx cluster", RED, 1);
    
    const CartesianVector projecPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(),avgShowerPosition, TPC_VIEW_W));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projecPosition, "avg shower position", BLUE, 1);
    
    std::cout << "METRIC: " << metric << std::endl;
    */
    
    
    
    //PandoraMonitoringApi::Pause(this->GetPandora());
    
    return metric;
        
    /*
    // ATTN: use the already-calculated sliding fits instead!
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
    const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
        
    // Loop over Hits.
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    float avgDeviation(0.f);
    for (const CaloHit *const pCaloHit : caloHitVector)
    {
        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), rL, rT);
        
        CartesianVector fitPosition(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalPosition(rL, 0.f, fitPosition);
        
        avgDeviation += (fitPosition - pCaloHit->GetPositionVector()).GetMagnitude();
    }
    
    avgDeviation /= static_cast<float>(caloHitVector.size());
    
    return (avgDeviation > m_maxTrackLikeDeviation);
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EnergyKickVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    return (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyAsymmetryParameters(const float weight, const CartesianVector &clusterDirection,
     CartesianVector &localWeightedDirectionSum) const
{
    // If the new axis direction is at an angle of greater than 90 deg to the current axis direction, flip it 180 degs.
    CartesianVector newDirection(clusterDirection);

    if (localWeightedDirectionSum.GetMagnitudeSquared() > std::numeric_limits<float>::epsilon())
    {
        if (localWeightedDirectionSum.GetCosOpeningAngle(clusterDirection) < 0.f)
            newDirection *= -1.f;
    }

    localWeightedDirectionSum += newDirection * weight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyAsymmetry(float &energyAsymmetry, float &localEnergyAsymmetry, ClusterList &localAsymmetryClusterList, const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const SlidingFitDataList &slidingFitDataList, const CartesianVector &localWeightedDirectionSum) const
{
    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxHitEnergy(0.f), afterVtxHitEnergy(0.f);
    unsigned int beforeVtxHitCount(0), afterVtxHitCount(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    for (const SlidingFitData &slidingFitData : slidingFitDataList)
    {
        const Cluster *const pCluster(slidingFitData.GetCluster());
        
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
    
    for (auto iter = localAsymmetryClusterList.begin(); iter != localAsymmetryClusterList.end(); /* no increment */)
    {
        if ((*iter)->GetNCaloHits() < m_minLocalAsymmetryClusterHits)
            iter = localAsymmetryClusterList.erase(iter);
            
        else
            ++iter;
    }
    
    if (localAsymmetryClusterList.size() == 1)
    {
        //std::vector<std::pair<float, float>> energyAsymmetryVec;
        
        float localBeforeVtxHitEnergy(0.f), localAfterVtxHitEnergy(0.f);
        unsigned int localBeforeVtxHitCount(0), localAfterVtxHitCount(0);
        
        for (const Cluster *const pCluster : localAsymmetryClusterList)
        {            
            float thisLocalBeforeVtxHitEnergy(0.f), thisLocalAfterVtxHitEnergy(0.f);
            unsigned int thisLocalBeforeVtxHitCount(0), thisLocalAfterVtxHitCount(0);
            
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
            std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

            for (const CaloHit *const pCaloHit : caloHitVector)
            {
                const float projectedPos = pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection);
                
                if (projectedPos < evtProjectedVtxPos)
                {
                    thisLocalBeforeVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++thisLocalBeforeVtxHitCount;
                }
                
                else if (projectedPos > evtProjectedVtxPos)
                {
                    thisLocalAfterVtxHitEnergy += pCaloHit->GetElectromagneticEnergy();
                    ++thisLocalAfterVtxHitCount;
                }
            }
            /*
            if ((useEnergyMetrics && (localAfterVtxHitEnergy < std::numeric_limits<float>::epsilon() || localBeforeVtxHitEnergy < std::numeric_limits<float>::epsilon())) ||
                (!useEnergyMetrics && (localAfterVtxHitCount == 0 || localBeforeVtxHitCount == 0)))
             {
                 continue;
             }
            
            if (useEnergyMetrics)
                energyAsymmetryVec.emplace_back(std::fabs((localAfterVtxHitEnergy - localBeforeVtxHitEnergy)), localAfterVtxHitEnergy + localBeforeVtxHitEnergy);
                
            else
            {
                if (0 == (localAfterVtxHitCount + localBeforeVtxHitCount))
                    throw StatusCodeException(STATUS_CODE_FAILURE);
                    
                energyAsymmetryVec.emplace_back(std::fabs((static_cast<float>(localAfterVtxHitCount) - static_cast<float>(localBeforeVtxHitCount))), 
                                                static_cast<float>(localAfterVtxHitCount) + static_cast<float>(localBeforeVtxHitCount));
            }
            */
            
            if ((useEnergyMetrics && (thisLocalAfterVtxHitEnergy < std::numeric_limits<float>::epsilon() || thisLocalBeforeVtxHitEnergy < std::numeric_limits<float>::epsilon())) ||
                (!useEnergyMetrics && (thisLocalAfterVtxHitCount == 0 || thisLocalBeforeVtxHitCount == 0)))
             {
                 continue;
             }
             
             localAfterVtxHitEnergy += thisLocalAfterVtxHitEnergy;
             localBeforeVtxHitEnergy += thisLocalBeforeVtxHitEnergy;
             localAfterVtxHitCount += thisLocalAfterVtxHitCount;
             localBeforeVtxHitCount += thisLocalBeforeVtxHitCount;
        }
        
        if (useEnergyMetrics)
        {
            if ((localAfterVtxHitEnergy + localBeforeVtxHitEnergy) > std::numeric_limits<float>::epsilon() && localAfterVtxHitEnergy > std::numeric_limits<float>::epsilon() && localBeforeVtxHitEnergy > std::numeric_limits<float>::epsilon())
                localEnergyAsymmetry += std::fabs((localAfterVtxHitEnergy - localBeforeVtxHitEnergy)) / localAfterVtxHitEnergy + localBeforeVtxHitEnergy;
        }
        
        else
        {
            if (0 == (localAfterVtxHitCount + localBeforeVtxHitCount))
                throw StatusCodeException(STATUS_CODE_FAILURE);
            
            if (localAfterVtxHitCount > 0 && localBeforeVtxHitCount > 0)
            {
                localEnergyAsymmetry += std::fabs((static_cast<float>(localAfterVtxHitCount) - static_cast<float>(localBeforeVtxHitCount))) 
                                                    / static_cast<float>(localAfterVtxHitCount) + static_cast<float>(localBeforeVtxHitCount);
            }
        }
        
        /*
        if (!energyAsymmetryVec.empty())
        {
            std::sort(energyAsymmetryVec.begin(), energyAsymmetryVec.end(), [](const std::pair<float, float> &lhs, const std::pair<float, float> &rhs){ return lhs.first > rhs.first; });
            localEnergyAsymmetry += energyAsymmetryVec.front().first / energyAsymmetryVec.front().second;
        }
        */
    }
    
    // Use energy metrics if possible, otherwise fall back on hit counting.
    const float totHitEnergy(beforeVtxHitEnergy + afterVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);
    
    if (useEnergyMetrics)
        energyAsymmetry += std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;
    
    else
    {
        if (0 == totHitCount)
            throw StatusCodeException(STATUS_CODE_FAILURE);
            
        energyAsymmetry += std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::SlidingFitData::SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, 
        const float slidingFitPitch) :
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f),
    m_pCluster(pCluster)
{
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
    m_minLayerDirection = slidingFitResult.GetGlobalMinLayerDirection();
    m_maxLayerDirection = slidingFitResult.GetGlobalMaxLayerDirection();
    m_minLayerPosition = slidingFitResult.GetGlobalMinLayerPosition();
    m_maxLayerPosition = slidingFitResult.GetGlobalMaxLayerPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EnergyKickVertexSelectionAlgorithm::ShowerCluster::ShowerCluster(const pandora::ClusterList &clusterList) :
    m_clusterList(clusterList)
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
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LocalAsymmetryConstant", m_localAsymmetryConstant));

    if ((m_rOffset < std::numeric_limits<float>::epsilon()) || (m_epsilon < std::numeric_limits<float>::epsilon()) || (m_asymmetryConstant < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "EnergyKickVertexSelection: Invalid parameter(s), ROffset " << m_rOffset << ", Epsilon " << m_epsilon << ", AsymmetryConstant " << m_asymmetryConstant << ", LocalAsymmetryConstant " << m_localAsymmetryConstant << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseLocalAsymmetry", m_useLocalAsymmetry));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLocalAsymmetryClusterHits", m_minLocalAsymmetryClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamDeweightingConstant", m_beamDeweightingConstant));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerDeweightingConstant", m_showerDeweightingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerCollapsingConstant", m_showerCollapsingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseClusterCharacterisation", m_useClusterCharacterisation));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackLikeDeviation", m_maxTrackLikeDeviation));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAngleConstant", m_showerAngleConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerDistanceConstant", m_showerDistanceConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexClusterDistance", m_vertexClusterDistance));
        
        
        
        

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
