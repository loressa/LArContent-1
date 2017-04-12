/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the SVM vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArSVMHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h" 

#include "larpandoracontent/LArVertex/EnergyKickFeatureTool.h"
#include "larpandoracontent/LArVertex/LocalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/GlobalAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/ShowerAsymmetryFeatureTool.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"

#include "larpandoracontent/LArVertex/SVMVertexSelectionAlgorithm.h"

#include <chrono>

#define COUT(a) std::cout << "\033[1;32m" << a << "\033[0m" << std::endl; // ATTN temporary

using namespace pandora;

namespace lar_content
{

SVMVertexSelectionAlgorithm::SVMVertexSelectionAlgorithm() : VertexSelectionBaseAlgorithm(),
    m_classifyUsingPermutations(true),
    m_trainingSetMode(false),
    m_allowClassifyDuringTraining(false),
    m_produceAllTrainingPermutations(false),
    m_mcVertexXCorrection(0.f),
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_minShowerSpineLength(15.f),
    m_beamDeweightingConstant(0.4),
    m_localAsymmetryConstant(3.f),
    m_globalAsymmetryConstant(1.f),
    m_showerAsymmetryConstant(1.f),
    m_energyKickConstant(0.06),
    m_minTopNSeparation(3.f),
    m_topNSize(5),
    m_showerClusteringDistance(3.f),
    m_minShowerClusterHits(1),
    m_useShowerClusteringApproximation(false),
    m_drawThings(false),         // ATTN temporary 
    m_cheatTrackShowerId(false), // ATTN temporary
    m_cheatTheVertex(false)      // ATTN temporary
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{    
    ClusterList clustersU, clustersV, clustersW;
    this->GetClusterLists(m_inputClusterListNames, clustersU, clustersV, clustersW);
        
    SlidingFitDataList slidingFitDataListU, slidingFitDataListV, slidingFitDataListW;
    this->CalculateClusterSlidingFits(clustersU, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, m_minClusterCaloHits, m_slidingFitWindow, slidingFitDataListW);
    
    ShowerClusterList showerClusterListU, showerClusterListV, showerClusterListW;
    this->CalculateShowerClusterList(clustersU, showerClusterListU);
    this->CalculateShowerClusterList(clustersV, showerClusterListV);
    this->CalculateShowerClusterList(clustersW, showerClusterListW);
    
    // Create maps from hit types to objects for passing to feature tools
    ClusterListMap clusterListMap{{TPC_VIEW_U, clustersU}, 
                                  {TPC_VIEW_V, clustersV},
                                  {TPC_VIEW_W, clustersV}};
                                  
    SlidingFitDataListMap slidingFitDataListMap{{TPC_VIEW_U, slidingFitDataListU}, 
                                                {TPC_VIEW_V, slidingFitDataListV},
                                                {TPC_VIEW_W, slidingFitDataListW}};

    ShowerClusterListMap showerClusterListMap{{TPC_VIEW_U, showerClusterListU}, 
                                              {TPC_VIEW_V, showerClusterListV},
                                              {TPC_VIEW_W, showerClusterListW}};

    KDTreeMap kdTreeMap{{TPC_VIEW_U, kdTreeU},
                        {TPC_VIEW_V, kdTreeV},
                        {TPC_VIEW_W, kdTreeW}};
                        
    EventFeatureInfo eventFeatureInfo = this->CalculateEventFeatures(clustersU, clustersV, clustersW, vertexVector);
    
    VertexFeatureInfoMap vertexFeatureInfoMap;
    for (const Vertex * const pVertex : vertexVector)
        this->PopulateVertexFeatureInfoMap(beamConstants, clusterListMap, slidingFitDataListMap, showerClusterListMap, kdTreeMap, pVertex, vertexFeatureInfoMap);
    
    const bool allowedToClassify(!m_trainingSetMode || m_allowClassifyDuringTraining);
    
    if (m_cheatTheVertex) // ATTN temporary
    {
        const Vertex *pCheatedVertex(NULL);
        this->GetCheatedVertex(vertexVector, pCheatedVertex);
        
        if (pCheatedVertex)
            vertexScoreList.emplace_back(pCheatedVertex, 1.f);
            
        return;
    }
    
    VertexScoreList bestVertexScoreList;
    if (m_trainingSetMode || !m_classifyUsingPermutations)
    {
        // Get the top-N vertices and find the best one using MC info if in training set mode
        VertexScoreList initialScoreList;
        for (const Vertex * const pVertex : vertexVector)
            this->PopulateInitialScoreList(vertexFeatureInfoMap, pVertex, initialScoreList);
        
        VertexList topNVertices;
        this->GetTopNVertices(initialScoreList, topNVertices);
        
        const Vertex *pBestVertex(NULL);
        float bestVertexDr(std::numeric_limits<float>::max());
        
        if (m_trainingSetMode)
            this->GetBestVertex(topNVertices, pBestVertex, bestVertexDr);
        
        for (const Vertex * const pVertex : topNVertices)
        {
            if (m_drawThings) // ATTN temporary
            {
                static unsigned int vertexCounter(0);
                const CartesianVector vertexPositionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));
                const std::string label = "Top-" + std::to_string(m_topNSize) + " vtx " + std::to_string(vertexCounter++);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertexPositionW, label, RED, 2);
            }
            
            FloatVector featureList = this->GenerateFeatureList(eventFeatureInfo, vertexFeatureInfoMap, pVertex, topNVertices);
            if (allowedToClassify && featureList.size() != m_svMachine.GetNumFeatures())
            {
                std::cout << "SVMVertexSelectionAlgorithm: the number of features (" << featureList.size() << ") did not match the expected " << 
                             "number of features (" << m_svMachine.GetNumFeatures() << ")" << std::endl;
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            }
            
            if (allowedToClassify)
            {
                // ATTN temporary
                if (false)
                {
                    std::cout << "Features: ";
                    for (const float feature : featureList)
                        std::cout << feature << " ";
                        
                    std::cout << std::endl;
                    
                    std::cout << " => Class score: " << SVMHelper::CalculateClassificationScore(m_svMachine, featureList) << std::endl;
                    std::cout << std::endl;
                    
                    PandoraMonitoringApi::Pause(this->GetPandora());
                }
                
                bestVertexScoreList.emplace_back(pVertex, SVMHelper::CalculateClassificationScore(m_svMachine, featureList));
            }
                
            if (m_trainingSetMode && pBestVertex && (topNVertices.size() == m_topNSize) && (bestVertexDr < m_minTopNSeparation))
            {
                if (m_produceAllTrainingPermutations)
                {
                    for (const FloatVector &permutedFeatureList : this->GeneratePermutedFeatureLists(eventFeatureInfo, vertexFeatureInfoMap, pVertex, topNVertices))
                        SVMHelper::ProduceTrainingExample(m_trainingOutputFile, (pVertex == pBestVertex), permutedFeatureList);
                }
                
                else
                    SVMHelper::ProduceTrainingExample(m_trainingOutputFile, (pVertex == pBestVertex), featureList);
            }
        }
    }
    
    else if (m_classifyUsingPermutations && allowedToClassify && !vertexVector.empty())
    {
        VertexVector currentVertexSet(vertexVector);
        while (currentVertexSet.size() > 1)
        {
            VectorOfVertexVectors setsOfNVertices;
            auto iter = currentVertexSet.begin();
            
            while (true)
            {
                VertexVector nVertexSet;
                while (iter != currentVertexSet.end())
                {
                    if (nVertexSet.size() < m_topNSize)
                    {
                        nVertexSet.push_back(*iter);
                        ++iter;
                    }
                        
                    else
                        break;
                }
                
                if (nVertexSet.empty())
                    break;
                
                const bool endOfVector = (nVertexSet.size() < m_topNSize);
                
                while (nVertexSet.size() < m_topNSize)
                    nVertexSet.push_back(nVertexSet.back());
                    
                setsOfNVertices.push_back(nVertexSet);
                    
                if (endOfVector)
                    break;
            }
            
            currentVertexSet.clear();
            
            for (VertexVector &vertexSet : setsOfNVertices)
            {
                VertexScoreList nSetScoreList;
                for (int i = 0; i < m_topNSize; ++i)
                {
                    const Vertex *const pVertex(vertexSet.front());
                    const VertexList vertexList(vertexSet.begin(), vertexSet.end());
                    
                    FloatVector featureList(this->GenerateFeatureList(eventFeatureInfo, vertexFeatureInfoMap, pVertex, vertexList));
                    if (featureList.size() != m_svMachine.GetNumFeatures())
                    {
                        std::cout << "SVMVertexSelectionAlgorithm: the number of features (" << featureList.size() << ") did not match the expected " << 
                                     "number of features (" << m_svMachine.GetNumFeatures() << ")" << std::endl;
                        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
                    }
                    
                    nSetScoreList.emplace_back(pVertex, SVMHelper::CalculateClassificationScore(m_svMachine, featureList));
                    std::rotate(vertexSet.begin(), std::next(vertexSet.begin(), 1), vertexSet.end());
                }
                
                std::sort(nSetScoreList.begin(), nSetScoreList.end());
                currentVertexSet.push_back(nSetScoreList.front().GetVertex());
            }
        }
        
        if (currentVertexSet.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);
        
        bestVertexScoreList.emplace_back(currentVertexSet.front(), 0.f);
    }
    
    if (!allowedToClassify)
        throw StatusCodeException(STATUS_CODE_SUCCESS);
    
    // Now some fine-tuning using only the RPhi score.
    this->PopulateFinalVertexScoreList(vertexFeatureInfoMap, bestVertexScoreList, vertexVector, vertexScoreList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::CalculateShowerClusterList(const ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const
{
    ClusterEndPointsMap clusterEndPointsMap;
    ClusterList showerLikeClusters;
    this->GetShowerLikeClusterEndPoints(inputClusterList, showerLikeClusters, clusterEndPointsMap);
    
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    ClusterList availableShowerLikeClusters(showerLikeClusters.begin(), showerLikeClusters.end());
        
    while (!availableShowerLikeClusters.empty())
    {
        ClusterList showerCluster;
        showerCluster.push_back(availableShowerLikeClusters.back());
        availableShowerLikeClusters.pop_back();
        
        bool addedCluster(true);
        while (addedCluster && !availableShowerLikeClusters.empty())
        {
            addedCluster = false;
            for (const Cluster * const pCluster : showerCluster)
            {
                addedCluster = this->AddClusterToShower(clusterEndPointsMap, availableShowerLikeClusters, pCluster, showerCluster);
                
                if (addedCluster)
                    break;
            }
        }
        
        int totHits(0);
        for (const Cluster * const pCluster : showerCluster)
            totHits += pCluster->GetNCaloHits();
        
        if (totHits < m_minClusterCaloHits)
            continue;
        
        showerClusterList.emplace_back(showerCluster, slidingFitPitch, m_slidingFitWindow);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::GetShowerLikeClusterEndPoints(const ClusterList &clusterList, ClusterList &showerLikeClusters, 
    ClusterEndPointsMap &clusterEndPointsMap) const
{    
    for (const Cluster *const pCluster : clusterList)
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
        
        ClusterEndPoints clusterEndPoints(clusterCaloHitVector.front()->GetPositionVector(), clusterCaloHitVector.back()->GetPositionVector());
        clusterEndPointsMap.emplace(pCluster, clusterEndPoints);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SVMVertexSelectionAlgorithm::AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, ClusterList &availableShowerLikeClusters,
    const Cluster *const pCluster, ClusterList &showerCluster) const
{
    const auto existingEndPointsIter(clusterEndPointsMap.find(pCluster));
    if (existingEndPointsIter == clusterEndPointsMap.end())
        return false;
        
    const ClusterEndPoints &existingClusterEndPoints(existingEndPointsIter->second);
    
    for (auto iter = availableShowerLikeClusters.begin(); iter != availableShowerLikeClusters.end(); /* no increment */)
    {
        const auto &newEndPointsIter(clusterEndPointsMap.find(*iter));
        if (newEndPointsIter == clusterEndPointsMap.end())
            continue;
            
        const ClusterEndPoints &newClusterEndPoints(newEndPointsIter->second);
        bool satisfiesClosestDistance(true);
        
        if (m_useShowerClusteringApproximation)
        {
            const float startStartDistance((newClusterEndPoints.first - existingClusterEndPoints.first).GetMagnitude());
            const float startEndDistance((newClusterEndPoints.first - existingClusterEndPoints.second).GetMagnitude());
            const float endStartDistance((newClusterEndPoints.second - existingClusterEndPoints.first).GetMagnitude());
            const float endEndDistance((newClusterEndPoints.second - existingClusterEndPoints.second).GetMagnitude());

            const float smallestDistance(std::min(std::min(startStartDistance, startEndDistance), std::min(endStartDistance, endEndDistance)));
            satisfiesClosestDistance = (smallestDistance < m_showerClusteringDistance);
        }
        
        else
            satisfiesClosestDistance = (LArClusterHelper::GetClosestDistance(pCluster, *iter) < m_showerClusteringDistance);
        
        if (satisfiesClosestDistance)
        {
            showerCluster.push_back(*iter);
            availableShowerLikeClusters.erase(iter);
            return true;
        }
        
        ++iter;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

SVMVertexSelectionAlgorithm::EventFeatureInfo SVMVertexSelectionAlgorithm::CalculateEventFeatures(const ClusterList &clusterListU, 
    const ClusterList &clusterListV, const ClusterList &clusterListW, const VertexVector &vertexVector) const
{
    float eventEnergy(0.f);
    unsigned int nShoweryHits(0), nHits(0);
    
    this->IncrementShoweryParameters(clusterListU, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListV, nShoweryHits, nHits, eventEnergy);
    this->IncrementShoweryParameters(clusterListW, nShoweryHits, nHits, eventEnergy);
    
    const unsigned int nClusters(clusterListU.size() + clusterListV.size() + clusterListW.size());
    const float eventShoweryness = (nHits > 0) ? static_cast<float>(nShoweryHits) / static_cast<float>(nHits) : 0.f;
    
    const float xSpan = this->GetCandidateSpan(vertexVector, [](const Vertex *const pVertex){ return pVertex->GetPosition().GetX(); });
    const float ySpan = this->GetCandidateSpan(vertexVector, [](const Vertex *const pVertex){ return pVertex->GetPosition().GetY(); });
    const float zSpan = this->GetCandidateSpan(vertexVector, [](const Vertex *const pVertex){ return pVertex->GetPosition().GetZ(); });
    
    float eventVolume(0.f), longitudinality(0.f);
    if ((xSpan != 0.f) && (ySpan != 0.f)) // ySpan often 0 - to be investigated
    {
        eventVolume     = xSpan * ySpan * zSpan;
        longitudinality = zSpan / (xSpan + ySpan);
    }
    
    else if (xSpan == 0.f)
    {
        eventVolume     = ySpan * ySpan * zSpan;
        longitudinality = zSpan / (ySpan + ySpan);
    }
    
    else
    {
        eventVolume     = xSpan * xSpan * zSpan;
        longitudinality = zSpan / (xSpan + xSpan);
    }
    
    return EventFeatureInfo(eventShoweryness, eventEnergy, eventVolume, longitudinality, nHits, nClusters, vertexVector.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::IncrementShoweryParameters(const ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits,
    float &eventEnergy) const
{
    for (const Cluster * const pCluster : clusterList)
    {
        if (this->IsClusterShowerLike(pCluster))
            nShoweryHits += pCluster->GetNCaloHits();
        
        eventEnergy += pCluster->GetElectromagneticEnergy();
        nHits += pCluster->GetNCaloHits();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SVMVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    bool isClusterShowerLike(false);
    
    if (m_cheatTrackShowerId) // ATTN temporary
    {          
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));                                                                           
        isClusterShowerLike = (PHOTON == pMCParticle->GetParticleId()) || (E_MINUS == std::abs(pMCParticle->GetParticleId())); 
    }                                                                                                                                                           

    else
        isClusterShowerLike = (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
        
    return isClusterShowerLike;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float SVMVertexSelectionAlgorithm::GetCandidateSpan(const VertexVector &vertexVector, const std::function<float(const Vertex *const)> &getCoord) const
{
    float coordMin(std::numeric_limits<float>::max());
    float coordMax(std::numeric_limits<float>::min());
    
    for (const Vertex *const pVertex : vertexVector)
    {
        const float coord(getCoord(pVertex));
        if (coord < coordMin)
            coordMin = coord;
            
        if (coord > coordMax)
            coordMax = coord;
    }
    
    if ((coordMax > std::numeric_limits<float>::min()) && (coordMin < std::numeric_limits<float>::max()))
        return std::fabs(coordMax - coordMin);
        
    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, const ClusterListMap &clusterListMap,
    const SlidingFitDataListMap &slidingFitDataListMap, const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap, 
    const Vertex * const pVertex, VertexFeatureInfoMap &vertexFeatureInfoMap) const
{
    float bestFastScore(0.f); // not actually used - artefact of toolizing RPhi score and still using performance trick
    
    const float beamDeweighting(this->GetBeamDeweightingScore(beamConstants, pVertex));
    
    const float energyKick(SVMHelper::CalculateFeature<EnergyKickFeatureTool>(m_featureToolMap, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore));
        
    const float localAsymmetry(SVMHelper::CalculateFeature<LocalAsymmetryFeatureTool>(m_featureToolMap, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore));
        
    const float globalAsymmetry(SVMHelper::CalculateFeature<GlobalAsymmetryFeatureTool>(m_featureToolMap, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore));
        
    const float showerAsymmetry(SVMHelper::CalculateFeature<ShowerAsymmetryFeatureTool>(m_featureToolMap, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore));
        
    const float rPhiFeature(SVMHelper::CalculateFeature<RPhiFeatureTool>(m_featureToolMap, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore));
    
    VertexFeatureInfo vertexFeatureInfo(beamDeweighting, rPhiFeature, energyKick, localAsymmetry, globalAsymmetry, showerAsymmetry);
    vertexFeatureInfoMap.emplace(pVertex, vertexFeatureInfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex * const pVertex,
                                                           VertexScoreList &initialScoreList) const
{
    VertexFeatureInfo vertexFeatureInfo = vertexFeatureInfoMap.at(pVertex);
    
    const float beamDeweightingScore(vertexFeatureInfo.m_beamDeweighting / m_beamDeweightingConstant);
    const float energyKickScore(-vertexFeatureInfo.m_energyKick / m_energyKickConstant);
    const float localAsymmetryScore(vertexFeatureInfo.m_localAsymmetry / m_localAsymmetryConstant);
    const float globalAsymmetryScore(vertexFeatureInfo.m_globalAsymmetry / m_globalAsymmetryConstant);
    const float showerAsymmetryScore(vertexFeatureInfo.m_showerAsymmetry / m_showerAsymmetryConstant);

    initialScoreList.emplace_back(pVertex, beamDeweightingScore + energyKickScore + localAsymmetryScore + globalAsymmetryScore + showerAsymmetryScore);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::GetTopNVertices(VertexScoreList &initialScoreList, VertexList &topNVertices) const
{
    std::sort(initialScoreList.begin(), initialScoreList.end());
    
    for (const VertexScore &vertexScore : initialScoreList)
    {
        const Vertex * const pVertex = vertexScore.GetVertex();        
        bool farEnoughAway(true);
        
        for (const Vertex * const pTopNVertex: topNVertices)
        {
            if (pTopNVertex == pVertex)
            {
                farEnoughAway = false;
                break;
            }
            
            const float distance = (pTopNVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude();
            
            if (distance <= m_minTopNSeparation)
            {
                farEnoughAway = false;
                break;
            }
        }
        
        if (farEnoughAway)
            topNVertices.push_back(pVertex);
            
        if (topNVertices.size() >= m_topNSize)
            return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::GetBestVertex(const VertexList &topNVertices, const Vertex *&pBestVertex, float &bestVertexDr) const
{
    // Extract input collections
    const MCParticleList *pMCParticleList = nullptr;
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);
    
    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
    
    for (const Vertex * const pVertex : topNVertices)
    {
        float mcVertexDr(std::numeric_limits<float>::max());
        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint().GetX() + m_mcVertexXCorrection, pMCNeutrino->GetEndpoint().GetY(), 
                pMCNeutrino->GetEndpoint().GetZ());
                
            const float dr = (mcNeutrinoPosition - pVertex->GetPosition()).GetMagnitude();
            if (dr < mcVertexDr)
                mcVertexDr = dr;
        }
        
        if (mcVertexDr < bestVertexDr)
        {
            bestVertexDr = mcVertexDr;
            pBestVertex = pVertex;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::GetCheatedVertex(const VertexVector &vertexVector, const Vertex *&pCheatedVertex) const // ATTN temporary
{
    // Extract input collections
    const MCParticleList *pMCParticleList = nullptr;
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);
    
    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
    
    float bestVertexDr(std::numeric_limits<float>::max());
    for (const Vertex * const pVertex : vertexVector)
    {
        float mcVertexDr(std::numeric_limits<float>::max());
        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint().GetX() + m_mcVertexXCorrection, pMCNeutrino->GetEndpoint().GetY(), 
                pMCNeutrino->GetEndpoint().GetZ());
                
            const float dr = (mcNeutrinoPosition - pVertex->GetPosition()).GetMagnitude();
            if (dr < mcVertexDr)
                mcVertexDr = dr;
        }
        
        if (mcVertexDr < bestVertexDr)
        {
            bestVertexDr = mcVertexDr;
            pCheatedVertex = pVertex;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

FloatVector SVMVertexSelectionAlgorithm::GenerateFeatureList(const EventFeatureInfo &eventFeatureInfo, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const Vertex * const pVertex, const VertexList &topNVertices) const
{
    FloatVector featureList;
    this->AddEventFeaturesToVector(eventFeatureInfo, featureList);
    
    VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
    this->AddVertexFeaturesToVector(vertexFeatureInfo, featureList);

    bool encounteredChosenVertex(false);
    for (const Vertex * const pTopNVertex : topNVertices)
    {
        if (pVertex == pTopNVertex && !encounteredChosenVertex)
        {
            encounteredChosenVertex = true;
            continue;
        }
           
        VertexFeatureInfo otherVertexFeatureInfo(vertexFeatureInfoMap.at(pTopNVertex));
        this->AddVertexFeaturesToVector(otherVertexFeatureInfo, featureList);
    }
    
    
    if (topNVertices.size() < m_topNSize)
    {
        for (std::size_t i = topNVertices.size(); i < m_topNSize; ++i)
        {
            for (const Vertex * const pTopNVertex : topNVertices)
            {
                if (pVertex == pTopNVertex)
                    continue;
                   
                VertexFeatureInfo otherVertexFeatureInfo(vertexFeatureInfoMap.at(pTopNVertex));
                this->AddVertexFeaturesToVector(otherVertexFeatureInfo, featureList);
                
                break;
            }
        }
    }
    
    return featureList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo, FloatVector &featureVector) const
{
    featureVector.push_back(eventFeatureInfo.m_eventShoweryness);
    featureVector.push_back(eventFeatureInfo.m_eventEnergy);
    featureVector.push_back(eventFeatureInfo.m_eventVolume);
    featureVector.push_back(eventFeatureInfo.m_longitudinality);
    featureVector.push_back(static_cast<float>(eventFeatureInfo.m_nHits));
    featureVector.push_back(static_cast<float>(eventFeatureInfo.m_nClusters));
    featureVector.push_back(static_cast<float>(eventFeatureInfo.m_nCandidates));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo, FloatVector &featureVector) const
{
    featureVector.push_back(vertexFeatureInfo.m_beamDeweighting);
    featureVector.push_back(vertexFeatureInfo.m_energyKick);
    featureVector.push_back(vertexFeatureInfo.m_globalAsymmetry);
    featureVector.push_back(vertexFeatureInfo.m_localAsymmetry);;
    featureVector.push_back(vertexFeatureInfo.m_showerAsymmetry);
    featureVector.push_back(vertexFeatureInfo.m_rPhiFeature);
}

//------------------------------------------------------------------------------------------------------------------------------------------

SVMVertexSelectionAlgorithm::FeatureListVector SVMVertexSelectionAlgorithm::GeneratePermutedFeatureLists(const EventFeatureInfo &eventFeatureInfo, 
    const VertexFeatureInfoMap &vertexFeatureInfoMap, const Vertex * const pVertex, const VertexList &topNVertices) const
{
    if (topNVertices.size() != m_topNSize)
    {
        std::cout << "SVMVertexSelectionAlgorithm: Can only generate permuted feature lists for full sets of top-N vertices" << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
    }
        
    FloatVector staticFeatureList;
    this->AddEventFeaturesToVector(eventFeatureInfo, staticFeatureList);
    
    VertexFeatureInfo vertexFeatureInfo(vertexFeatureInfoMap.at(pVertex));
    this->AddVertexFeaturesToVector(vertexFeatureInfo, staticFeatureList);

    // The rest of the vector can be permuted
    VertexVector vertexVector;
    for (const Vertex * const pTopNVertex : topNVertices)
    {
        if (pVertex != pTopNVertex)
            vertexVector.push_back(pTopNVertex);
    }
    
    std::sort(vertexVector.begin(), vertexVector.end(), std::less<const Vertex *>());
    
    FeatureListVector featureListVector;    
    do
    {
        FloatVector featureList(staticFeatureList);
        for (const Vertex * const pTopNVertex : vertexVector)
        {
            const auto featureInfoIter = vertexFeatureInfoMap.find(pTopNVertex);
            if (featureInfoIter == vertexFeatureInfoMap.end())
            {
                std::cout << "SVMVertexSelectionAlgorithm: Could not find vertex in feature info map" << std::endl;
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
                
            this->AddVertexFeaturesToVector(featureInfoIter->second, featureList);
        }
        
        featureListVector.push_back(featureList);
        
    } while (std::next_permutation(vertexVector.begin(), vertexVector.end(), std::less<const Vertex *>()));
    
    return featureListVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, VertexScoreList &bestVertexScoreList,
    const VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const
{   
    // Now some fine-tuning using only the RPhi score.
    std::sort(bestVertexScoreList.begin(), bestVertexScoreList.end());
    
    if (!bestVertexScoreList.empty())
    {
        const Vertex * const pFavouriteVertex(bestVertexScoreList.front().GetVertex());
        const CartesianVector faveVertexPosition(pFavouriteVertex->GetPosition());
        
        for (const Vertex * const pVertex : vertexVector)
        { 
            if ((pVertex->GetPosition() - faveVertexPosition).GetMagnitude() < m_minTopNSeparation)
            {
                const float rPhiScore(vertexFeatureInfoMap.at(pVertex).m_rPhiFeature);
                finalVertexScoreList.emplace_back(pVertex, rPhiScore);
            }
        }
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SVMVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));
    
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));
    
    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SVMHelper::AddFeatureToolToMap(pAlgorithmTool, m_featureToolMap));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClassifyUsingPermutations", m_classifyUsingPermutations));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProduceAllTrainingPermutations", m_produceAllTrainingPermutations));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));
        
    if (m_trainingSetMode && m_mcParticleListName.empty())
    {
        std::cout << "SVMVertexSelectionAlgorithm: MCParticleListName is required for training set mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingOutputFile", m_trainingOutputFile));
        
    if (m_trainingSetMode && m_trainingOutputFile.empty())
    {
        std::cout << "SVMVertexSelectionAlgorithm: TrainingOutputFile is required for training set mode" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AllowClassifyDuringTraining", m_allowClassifyDuringTraining));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ParameterInputFile", m_parameterInputFile));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SVMName", m_svmName));
        
    if ((!m_trainingSetMode || m_allowClassifyDuringTraining))
    {
        if (m_parameterInputFile.empty() || m_svmName.empty())
        {
            std::cout << "SVMVertexSelectionAlgorithm: ParameterInputFile and SVMName must be set if training set mode is off or we allow " <<
                         "classification during training" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
        
        m_svMachine.Initialize(m_parameterInputFile, m_svmName);
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTopNSeparation", m_minTopNSeparation));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TopNSize", m_topNSize));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamDeweightingConstant", m_beamDeweightingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LocalAsymmetryConstant", m_localAsymmetryConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GlobalAsymmetryConstant", m_globalAsymmetryConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerAsymmetryConstant", m_showerAsymmetryConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyKickConstant", m_energyKickConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerClusterHits", m_minShowerClusterHits));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusteringApproximation", m_useShowerClusteringApproximation));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCVertexXCorrection", m_mcVertexXCorrection));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DrawThings", m_drawThings)); // ATTN temporary
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CheatTrackShowerId", m_cheatTrackShowerId)); // ATTN temporary
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CheatTheVertex", m_cheatTheVertex)); // ATTN temporary

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}
                                                            
} // namespace lar_content
