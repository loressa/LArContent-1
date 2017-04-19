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


#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h" 
#include "larpandoracontent/LArObjects/LArMCParticle.h"

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
    m_cheatTheVertex(false)       // ATTN temporary
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
    
    //--------------------------------------------------------------------------------------------
    
    // Extract input collections
    const MCParticleList *pMCParticleList = nullptr;
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);
    
    const CaloHitList *pCaloHitList = nullptr;
    PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList);
    
    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    this->SelectTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
    
    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);
    
    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList selectedCaloHitList;
    this->SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList);
    
    // Remove shared hits where target particle deposits below threshold energy fraction
    CaloHitList goodCaloHitList;
    this->SelectGoodCaloHits(&selectedCaloHitList, mcToPrimaryMCMap, goodCaloHitList);

    // Obtain maps: [good hit -> primary mc particle], [primary mc particle -> list of good hits]
    LArMonitoringHelper::CaloHitToMCMap goodHitToPrimaryMCMap;
    LArMonitoringHelper::MCContributionMap mcToGoodTrueHitListMap;
    LArMonitoringHelper::GetMCParticleToCaloHitMatches(&goodCaloHitList, mcToPrimaryMCMap, goodHitToPrimaryMCMap, mcToGoodTrueHitListMap);
    
    //--------------------------------------------------------------------------------------------
    
    const Vertex *pCheatedVertex(NULL);
    std::string interactionType = this->GetCheatedVertex(vertexVector, pCheatedVertex, mcToGoodTrueHitListMap);
    
    if (m_cheatTheVertex) // ATTN temporary
    {    
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
        
        if (pCheatedVertex)
            topNVertices.push_back(pCheatedVertex);
        
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
            
            SupportVectorMachine::DoubleVector featureList = this->GenerateFeatureList(eventFeatureInfo, vertexFeatureInfoMap, pVertex, topNVertices);
            if (allowedToClassify && featureList.size() != m_svMachine.GetNFeatures())
            {
                std::cout << "SVMVertexSelectionAlgorithm: the number of features (" << featureList.size() << ") did not match the expected " << 
                             "number of features (" << m_svMachine.GetNFeatures() << ")" << std::endl;
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            }
            
            if (allowedToClassify)
            {
                // ATTN temporary
                if (true)
                {
                    std::cout << "Features: ";
                    for (const double feature : featureList)
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
                    for (const SupportVectorMachine::DoubleVector &permutedFeatureList : this->GeneratePermutedFeatureLists(eventFeatureInfo, vertexFeatureInfoMap, pVertex, topNVertices))
                        SVMHelper::ProduceTrainingExample(m_trainingOutputFile + "_" + interactionType + ".txt", (pVertex == pBestVertex), permutedFeatureList);
                }
                
                else
                    SVMHelper::ProduceTrainingExample(m_trainingOutputFile + "_" + interactionType + ".txt", (pVertex == pBestVertex), featureList);
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
                    
                    SupportVectorMachine::DoubleVector featureList(this->GenerateFeatureList(eventFeatureInfo, vertexFeatureInfoMap, pVertex, vertexList));
                    if (featureList.size() != m_svMachine.GetNFeatures())
                    {
                        std::cout << "SVMVertexSelectionAlgorithm: the number of features (" << featureList.size() << ") did not match the expected " << 
                                     "number of features (" << m_svMachine.GetNFeatures() << ")" << std::endl;
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
    
    const double beamDeweighting(this->GetBeamDeweightingScore(beamConstants, pVertex));
    
    const double energyKick(SVMHelper::CalculateFeaturesOfType<EnergyKickFeatureTool>(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
        
    const double localAsymmetry(SVMHelper::CalculateFeaturesOfType<LocalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
        
    const double globalAsymmetry(SVMHelper::CalculateFeaturesOfType<GlobalAsymmetryFeatureTool>(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
        
    const double showerAsymmetry(SVMHelper::CalculateFeaturesOfType<ShowerAsymmetryFeatureTool>(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
        
    const double rPhiFeature(SVMHelper::CalculateFeaturesOfType<RPhiFeatureTool>(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
        
    const double rPhiFeaturse(SVMHelper::CalculateFeatures(m_featureToolVector, this, pVertex, slidingFitDataListMap, 
        clusterListMap, kdTreeMap, showerClusterListMap, beamDeweighting, bestFastScore).at(0));
    
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

std::string SVMVertexSelectionAlgorithm::GetCheatedVertex(const VertexVector &vertexVector, const Vertex *&pCheatedVertex, 
    const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap) const // ATTN temporary
{
    // Extract input collections
    const MCParticleList *pMCParticleList = nullptr;
    PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList);
    
    // Obtain vector: true neutrinos
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);
    
    std::string interactionType("UNKNOWN");
    
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
            {
                const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

                if (pLArMCNeutrino)
                    interactionType = this->ToString(this->GetInteractionType(pLArMCNeutrino, pMCParticleList, mcToGoodTrueHitListMap));
                    
                mcVertexDr = dr;
                
            }
        }

        if (mcVertexDr < bestVertexDr)
        {
            bestVertexDr = mcVertexDr;
            pCheatedVertex = pVertex;
        }
    }
    
    return interactionType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

SupportVectorMachine::DoubleVector SVMVertexSelectionAlgorithm::GenerateFeatureList(const EventFeatureInfo &eventFeatureInfo, const VertexFeatureInfoMap &vertexFeatureInfoMap,
    const Vertex * const pVertex, const VertexList &topNVertices) const
{
    SupportVectorMachine::DoubleVector featureList;
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

void SVMVertexSelectionAlgorithm::AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo, SupportVectorMachine::DoubleVector &featureVector) const
{
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventShoweryness));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventEnergy));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_eventVolume));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_longitudinality));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nHits));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nClusters));
    featureVector.push_back(static_cast<double>(eventFeatureInfo.m_nCandidates));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo, SupportVectorMachine::DoubleVector &featureVector) const
{
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_beamDeweighting));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_energyKick));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_globalAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_localAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_showerAsymmetry));
    featureVector.push_back(static_cast<double>(vertexFeatureInfo.m_rPhiFeature));
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
        
    SupportVectorMachine::DoubleVector staticFeatureList;
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
        SupportVectorMachine::DoubleVector featureList(staticFeatureList);
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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SVMHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));
        
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
        
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}




//------------------------------------------------------------------------------------------------------------------------------------------


void SVMVertexSelectionAlgorithm::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedCaloHitList) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;

            if (this->PassMCParticleChecks(pPrimaryParticle, pPrimaryParticle, pHitParticle))
                selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SVMVertexSelectionAlgorithm::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle) const
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) && (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > 2.5f)
            return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (this->PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedGoodCaloHitList) const
{
    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap primaryWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pMCParticle);

            if (mcToPrimaryMCMap.end() != mcIter)
                primaryWeightMap[mcIter->second] += weight;
        }

        MCParticleVector mcPrimaryVector;
        for (const auto &mapEntry : primaryWeightMap) mcPrimaryVector.push_back(mapEntry.first);
        std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestPrimaryParticle(nullptr);
        float bestPrimaryWeight(0.f), primaryWeightSum(0.f);

        for (const MCParticle *const pPrimaryMCParticle : mcPrimaryVector)
        {
            const float primaryWeight(primaryWeightMap.at(pPrimaryMCParticle));
            primaryWeightSum += primaryWeight;

            if (primaryWeight > bestPrimaryWeight)
            {
                bestPrimaryWeight = primaryWeight;
                pBestPrimaryParticle = pPrimaryMCParticle;
            }
        }

        if (!pBestPrimaryParticle || (primaryWeightSum < std::numeric_limits<float>::epsilon()) || ((bestPrimaryWeight / primaryWeightSum) < 0.9f))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

SVMVertexSelectionAlgorithm::InteractionType SVMVertexSelectionAlgorithm::GetInteractionType(const LArMCParticle *const pLArMCNeutrino, const MCParticleList *pMCParticleList, const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap) const
{
    // Obtain vector: primary mc particles
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);
    
    unsigned int nNonNeutrons(0), nMuons(0), nElectrons(0), nProtons(0), nPiPlus(0), nPiMinus(0), nNeutrons(0), nPhotons(0);
    
    for (const auto pMCPrimary : mcPrimaryVector)
    {
        LArMonitoringHelper::MCContributionMap::const_iterator goodTrueHitsIter = mcToGoodTrueHitListMap.find(pMCPrimary);

        if (mcToGoodTrueHitListMap.end() != goodTrueHitsIter)
        {
            const CaloHitList &caloHitList(goodTrueHitsIter->second);
            if (caloHitList.size() < 15)
                continue;
                
            int nGoodViews(0);
            if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList) >= 5)
                ++nGoodViews;
            
            if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList) >= 5)
                ++nGoodViews; 
                
            if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList) >= 5)
                ++nGoodViews;
                
            if (nGoodViews < 2)
                continue;
                
        }
        
        //if (!this->IsGoodMCPrimary(pMCPrimary))
        //    continue;
        
        if (2112 != pMCPrimary->GetParticleId()) ++nNonNeutrons;

        if (13 == pMCPrimary->GetParticleId()) ++nMuons;
        if (11 == pMCPrimary->GetParticleId()) ++nElectrons;
        else if (2212 == pMCPrimary->GetParticleId()) ++nProtons;
        else if (22 == pMCPrimary->GetParticleId()) ++nPhotons;
        else if (211 == pMCPrimary->GetParticleId()) ++nPiPlus;
        else if (-211 == pMCPrimary->GetParticleId()) ++nPiMinus;
        else if (2112 == pMCPrimary->GetParticleId()) ++nNeutrons;
    }

    if (1098 == pLArMCNeutrino->GetNuanceCode()) return NU_E_SCATTERING;

    if (1001 == pLArMCNeutrino->GetNuanceCode())
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCQEL_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCQEL_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCQEL_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCQEL_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCQEL_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCQEL_MU_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCQEL_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCQEL_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCQEL_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCQEL_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCQEL_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCQEL_E_P_P_P_P_P;
    }

    if (1002 == pLArMCNeutrino->GetNuanceCode())
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCQEL_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCQEL_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCQEL_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCQEL_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCQEL_P_P_P_P_P;
    }

    if ((pLArMCNeutrino->GetNuanceCode() >= 1003) && (pLArMCNeutrino->GetNuanceCode() <= 1005))
    {
        if ((1 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons)) return CCRES_MU;
        if ((2 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons)) return CCRES_MU_P;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons)) return CCRES_MU_P_P;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons)) return CCRES_MU_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons)) return CCRES_MU_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons)) return CCRES_MU_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_MU_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_MU_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_MU_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_MU_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nMuons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_MU_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nMuons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nMuons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nMuons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nMuons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nMuons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_MU_P_P_P_P_P_PIZERO;

        if ((1 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons)) return CCRES_E;
        if ((2 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons)) return CCRES_E_P;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons)) return CCRES_E_P_P;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons)) return CCRES_E_P_P_P;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons)) return CCRES_E_P_P_P_P;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons)) return CCRES_E_P_P_P_P_P;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPiPlus)) return CCRES_E_PIPLUS;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_PIPLUS;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_PIPLUS;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPiPlus)) return CCRES_E_P_P_P_P_P_PIPLUS;

        if ((2 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (1 == nPhotons)) return CCRES_E_PHOTON;
        if ((3 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (1 == nPhotons)) return CCRES_E_P_PHOTON;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_PHOTON;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (1 == nPhotons)) return CCRES_E_P_P_P_P_P_PHOTON;

        if ((3 == nNonNeutrons) && (1 == nElectrons) && (0 == nProtons) && (2 == nPhotons)) return CCRES_E_PIZERO;
        if ((4 == nNonNeutrons) && (1 == nElectrons) && (1 == nProtons) && (2 == nPhotons)) return CCRES_E_P_PIZERO;
        if ((5 == nNonNeutrons) && (1 == nElectrons) && (2 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (1 == nElectrons) && (3 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (1 == nElectrons) && (4 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_PIZERO;
        if ((8 == nNonNeutrons) && (1 == nElectrons) && (5 == nProtons) && (2 == nPhotons)) return CCRES_E_P_P_P_P_P_PIZERO;
    }

    if ((pLArMCNeutrino->GetNuanceCode() >= 1006) && (pLArMCNeutrino->GetNuanceCode() <= 1009))
    {
        if ((1 == nNonNeutrons) && (1 == nProtons)) return NCRES_P;
        if ((2 == nNonNeutrons) && (2 == nProtons)) return NCRES_P_P;
        if ((3 == nNonNeutrons) && (3 == nProtons)) return NCRES_P_P_P;
        if ((4 == nNonNeutrons) && (4 == nProtons)) return NCRES_P_P_P_P;
        if ((5 == nNonNeutrons) && (5 == nProtons)) return NCRES_P_P_P_P_P;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiPlus)) return NCRES_PIPLUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiPlus)) return NCRES_P_PIPLUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_PIPLUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_PIPLUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_PIPLUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiPlus)) return NCRES_P_P_P_P_P_PIPLUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPiMinus)) return NCRES_PIMINUS;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPiMinus)) return NCRES_P_PIMINUS;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_PIMINUS;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_PIMINUS;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_PIMINUS;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPiMinus)) return NCRES_P_P_P_P_P_PIMINUS;

        if ((1 == nNonNeutrons) && (0 == nProtons) && (1 == nPhotons)) return NCRES_PHOTON;
        if ((2 == nNonNeutrons) && (1 == nProtons) && (1 == nPhotons)) return NCRES_P_PHOTON;
        if ((3 == nNonNeutrons) && (2 == nProtons) && (1 == nPhotons)) return NCRES_P_P_PHOTON;
        if ((4 == nNonNeutrons) && (3 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_PHOTON;
        if ((5 == nNonNeutrons) && (4 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_PHOTON;
        if ((6 == nNonNeutrons) && (5 == nProtons) && (1 == nPhotons)) return NCRES_P_P_P_P_P_PHOTON;

        if ((2 == nNonNeutrons) && (0 == nProtons) && (2 == nPhotons)) return NCRES_PIZERO;
        if ((3 == nNonNeutrons) && (1 == nProtons) && (2 == nPhotons)) return NCRES_P_PIZERO;
        if ((4 == nNonNeutrons) && (2 == nProtons) && (2 == nPhotons)) return NCRES_P_P_PIZERO;
        if ((5 == nNonNeutrons) && (3 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_PIZERO;
        if ((6 == nNonNeutrons) && (4 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_PIZERO;
        if ((7 == nNonNeutrons) && (5 == nProtons) && (2 == nPhotons)) return NCRES_P_P_P_P_P_PIZERO;
    }

    if (pLArMCNeutrino->GetNuanceCode() == 1091) return CCDIS;
    if (pLArMCNeutrino->GetNuanceCode() == 1092) return NCDIS;
    if (pLArMCNeutrino->GetNuanceCode() == 1096) return NCCOH;
    if (pLArMCNeutrino->GetNuanceCode() == 1097) return CCCOH;

    return OTHER_INTERACTION;
}

/**
 *  @brief  Get string representing interaction type
 * 
 *  @param  interactionType
 * 
 *  @return the interaction type string
 */
std::string SVMVertexSelectionAlgorithm::ToString(const InteractionType interactionType) const
{
    switch (interactionType)
    {
    case CCQEL_MU: return "CCQEL_MU";
    case CCQEL_MU_P: return "CCQEL_MU_P";
    case CCQEL_MU_P_P: return "CCQEL_MU_P_P";
    case CCQEL_MU_P_P_P: return "CCQEL_MU_P_P_P";
    case CCQEL_MU_P_P_P_P: return "CCQEL_MU_P_P_P_P";
    case CCQEL_MU_P_P_P_P_P: return "CCQEL_MU_P_P_P_P_P";
    case CCQEL_E: return "CCQEL_E";
    case CCQEL_E_P: return "CCQEL_E_P";
    case CCQEL_E_P_P: return "CCQEL_E_P_P";
    case CCQEL_E_P_P_P: return "CCQEL_E_P_P_P";
    case CCQEL_E_P_P_P_P: return "CCQEL_E_P_P_P_P";
    case CCQEL_E_P_P_P_P_P: return "CCQEL_E_P_P_P_P_P";
    case NCQEL_P: return "NCQEL_P";
    case NCQEL_P_P: return "NCQEL_P_P";
    case NCQEL_P_P_P: return "NCQEL_P_P_P";
    case NCQEL_P_P_P_P: return "NCQEL_P_P_P_P";
    case NCQEL_P_P_P_P_P: return "NCQEL_P_P_P_P_P";
    case CCRES_MU: return "CCRES_MU";
    case CCRES_MU_P: return "CCRES_MU_P";
    case CCRES_MU_P_P: return "CCRES_MU_P_P";
    case CCRES_MU_P_P_P: return "CCRES_MU_P_P_P";
    case CCRES_MU_P_P_P_P: return "CCRES_MU_P_P_P_P";
    case CCRES_MU_P_P_P_P_P: return "CCRES_MU_P_P_P_P_P";
    case CCRES_MU_PIPLUS: return "CCRES_MU_PIPLUS";
    case CCRES_MU_P_PIPLUS: return "CCRES_MU_P_PIPLUS";
    case CCRES_MU_P_P_PIPLUS: return "CCRES_MU_P_P_PIPLUS";
    case CCRES_MU_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_PIPLUS";
    case CCRES_MU_P_P_P_P_P_PIPLUS: return "CCRES_MU_P_P_P_P_P_PIPLUS";
    case CCRES_MU_PHOTON: return "CCRES_MU_PHOTON";
    case CCRES_MU_P_PHOTON: return "CCRES_MU_P_PHOTON";
    case CCRES_MU_P_P_PHOTON: return "CCRES_MU_P_P_PHOTON";
    case CCRES_MU_P_P_P_PHOTON: return "CCRES_MU_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_PHOTON";
    case CCRES_MU_P_P_P_P_P_PHOTON: return "CCRES_MU_P_P_P_P_P_PHOTON";
    case CCRES_MU_PIZERO: return "CCRES_MU_PIZERO";
    case CCRES_MU_P_PIZERO: return "CCRES_MU_P_PIZERO";
    case CCRES_MU_P_P_PIZERO: return "CCRES_MU_P_P_PIZERO";
    case CCRES_MU_P_P_P_PIZERO: return "CCRES_MU_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_PIZERO";
    case CCRES_MU_P_P_P_P_P_PIZERO: return "CCRES_MU_P_P_P_P_P_PIZERO";
    case CCRES_E: return "CCRES_E";
    case CCRES_E_P: return "CCRES_E_P";
    case CCRES_E_P_P: return "CCRES_E_P_P";
    case CCRES_E_P_P_P: return "CCRES_E_P_P_P";
    case CCRES_E_P_P_P_P: return "CCRES_E_P_P_P_P";
    case CCRES_E_P_P_P_P_P: return "CCRES_E_P_P_P_P_P";
    case CCRES_E_PIPLUS: return "CCRES_E_PIPLUS";
    case CCRES_E_P_PIPLUS: return "CCRES_E_P_PIPLUS";
    case CCRES_E_P_P_PIPLUS: return "CCRES_E_P_P_PIPLUS";
    case CCRES_E_P_P_P_PIPLUS: return "CCRES_E_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_PIPLUS";
    case CCRES_E_P_P_P_P_P_PIPLUS: return "CCRES_E_P_P_P_P_P_PIPLUS";
    case CCRES_E_PHOTON: return "CCRES_E_PHOTON";
    case CCRES_E_P_PHOTON: return "CCRES_E_P_PHOTON";
    case CCRES_E_P_P_PHOTON: return "CCRES_E_P_P_PHOTON";
    case CCRES_E_P_P_P_PHOTON: return "CCRES_E_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_PHOTON";
    case CCRES_E_P_P_P_P_P_PHOTON: return "CCRES_E_P_P_P_P_P_PHOTON";
    case CCRES_E_PIZERO: return "CCRES_E_PIZERO";
    case CCRES_E_P_PIZERO: return "CCRES_E_P_PIZERO";
    case CCRES_E_P_P_PIZERO: return "CCRES_E_P_P_PIZERO";
    case CCRES_E_P_P_P_PIZERO: return "CCRES_E_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_PIZERO";
    case CCRES_E_P_P_P_P_P_PIZERO: return "CCRES_E_P_P_P_P_P_PIZERO";
    case NCRES_P: return "NCRES_P";
    case NCRES_P_P: return "NCRES_P_P";
    case NCRES_P_P_P: return "NCRES_P_P_P";
    case NCRES_P_P_P_P: return "NCRES_P_P_P_P";
    case NCRES_P_P_P_P_P: return "NCRES_P_P_P_P_P";
    case NCRES_PIPLUS: return "NCRES_PIPLUS";
    case NCRES_P_PIPLUS: return "NCRES_P_PIPLUS";
    case NCRES_P_P_PIPLUS: return "NCRES_P_P_PIPLUS";
    case NCRES_P_P_P_PIPLUS: return "NCRES_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_PIPLUS";
    case NCRES_P_P_P_P_P_PIPLUS: return "NCRES_P_P_P_P_P_PIPLUS";
    case NCRES_PIMINUS: return "NCRES_PIMINUS";
    case NCRES_P_PIMINUS: return "NCRES_P_PIMINUS";
    case NCRES_P_P_PIMINUS: return "NCRES_P_P_PIMINUS";
    case NCRES_P_P_P_PIMINUS: return "NCRES_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_PIMINUS";
    case NCRES_P_P_P_P_P_PIMINUS: return "NCRES_P_P_P_P_P_PIMINUS";
    case NCRES_PHOTON: return "NCRES_PHOTON";
    case NCRES_P_PHOTON: return "NCRES_P_PHOTON";
    case NCRES_P_P_PHOTON: return "NCRES_P_P_PHOTON";
    case NCRES_P_P_P_PHOTON: return "NCRES_P_P_P_PHOTON";
    case NCRES_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_PHOTON";
    case NCRES_P_P_P_P_P_PHOTON: return "NCRES_P_P_P_P_P_PHOTON";
    case NCRES_PIZERO: return "NCRES_PIZERO";
    case NCRES_P_PIZERO: return "NCRES_P_PIZERO";
    case NCRES_P_P_PIZERO: return "NCRES_P_P_PIZERO";
    case NCRES_P_P_P_PIZERO: return "NCRES_P_P_P_PIZERO";
    case NCRES_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_PIZERO";
    case NCRES_P_P_P_P_P_PIZERO: return "NCRES_P_P_P_P_P_PIZERO";
    case CCDIS: return "CCDIS";
    case NCDIS: return "NCDIS";
    case CCCOH: return "CCCOH";
    case NCCOH: return "NCCOH";
    case OTHER_INTERACTION: return "OTHER_INTERACTION";
    case ALL_INTERACTIONS: return "ALL_INTERACTIONS";
    case NU_E_SCATTERING: return "NU_E_SCATTERING";
    default: return "UNKNOWN";
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void SVMVertexSelectionAlgorithm::SelectTrueNeutrinos(const MCParticleList *const pAllMCParticleList, MCParticleVector &selectedMCNeutrinoVector) const
{
    MCParticleVector allMCNeutrinoVector;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pAllMCParticleList, allMCNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : allMCNeutrinoVector)
    {
        // ATTN Now demand that input mc neutrinos LArMCParticles, with addition of interaction type
        const LArMCParticle *const pLArMCNeutrino(dynamic_cast<const LArMCParticle*>(pMCNeutrino));

        if (pLArMCNeutrino && (0 != pLArMCNeutrino->GetNuanceCode()))
            selectedMCNeutrinoVector.push_back(pMCNeutrino);
    }
}
                                                            
} // namespace lar_content
