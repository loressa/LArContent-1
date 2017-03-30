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

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h" 
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include <fstream>
#include <tuple>

using namespace pandora;

namespace lar_content
{

EnergyKickVertexSelectionAlgorithm::EnergyKickVertexSelectionAlgorithm() :
    m_minClusterCaloHits(12),
    m_slidingFitWindow(100),
    m_rOffset(10.f),
    m_xOffset(0.06),
    m_epsilon(0.06),
    m_maxAsymmetryDistance(5.f),
    m_minAsymmetryCosAngle(0.9962),
    m_maxAsymmetryNClusters(2),
    m_showerDeweightingConstant(1.f),
    m_showerCollapsingConstant(1.f),
    m_minShowerSpineLength(15.f),
    m_showerClusteringDistance(3.f),
    m_vertexClusterDistance(4.f),
    m_minShowerClusterHits(1),
    m_useShowerClusteringApproximation(false),    
    m_cheatTrackShowerId(false),
    
    m_fastScoreCheck(true),
    m_fastScoreOnly(false),
    m_fullScore(true),
    m_kernelEstimateSigma(0.048f),
    m_kappa(0.42f),
    m_maxHitVertexDisplacement1D(100.f),
    m_minFastScoreFraction(0.8f),
    m_fastHistogramNPhiBins(200),
    m_fastHistogramPhiMin(-1.1f * M_PI),
    m_fastHistogramPhiMax(+1.1f * M_PI),
    m_enableFolding(true),
    
    m_beamDeweightingConstant(0.4),
    m_globalAsymmetryConstant(1.f),
    m_showerAsymmetryConstant(1.f),
    m_asymmetryConstant(3.f),
    
    m_trainTypeOne(false),
    m_trainTypeTwo(false)
{
    m_supportVectors = this->GetSupportVectors();
    m_modifiedAlphas = this->GetModifiedAlphas();
    m_svMu = this->GetSVMu();
    m_svSigma = this->GetSVSigma();
    m_svBias = this->GetSVBias();
    m_svScale = this->GetSVScale();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants,
    HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    
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
    
    VertexInfoVector vertexInfoVector;
    VertexInfoVector allVerticesInfoVector;
    
    
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
 
    this->CalculateShowerClusterMap(clustersU, showerClusterMapU);
    this->CalculateShowerClusterMap(clustersV, showerClusterMapV);
    this->CalculateShowerClusterMap(clustersW, showerClusterMapW);
    
    SlidingFitDataList singleClusterSlidingFitDataListU, singleClusterSlidingFitDataListV, singleClusterSlidingFitDataListW;
    
    this->CalculateClusterSlidingFits(clustersU, singleClusterSlidingFitDataListU);
    this->CalculateClusterSlidingFits(clustersV, singleClusterSlidingFitDataListV);
    this->CalculateClusterSlidingFits(clustersW, singleClusterSlidingFitDataListW);
    
    const float eventHitShoweryness = this->GetEventHitShoweryness(clustersU, clustersV, clustersW);
    const float eventClusterShoweryness = this->GetEventClusterShoweryness(clustersU, clustersV, clustersW);
    
    EventInputs eventInputs(eventHitShoweryness, eventClusterShoweryness);
    
    float bestFastScore(0.f);
    float rPhiScore(0.f);
    
    for (const Vertex *const pVertex : vertexVector)
    {
        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float beamDeweightingScore(this->IsBeamModeOn() ? std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant()) : 1.f);
        const float beamDeweighting(this->IsBeamModeOn() ? -(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant() : 0.f);

        float energyKick(0.f), energyAsymmetry(0.f), globalEnergyAsymmetry(0.f), showerEnergyAsymmetry(0.f);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U), energyKick, energyAsymmetry, globalEnergyAsymmetry, showerClusterMapU, showerEnergyAsymmetry, singleClusterSlidingFitDataListU);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V), energyKick, energyAsymmetry, globalEnergyAsymmetry, showerClusterMapV, showerEnergyAsymmetry, singleClusterSlidingFitDataListV);
        this->IncrementEnergyScoresForView(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W), energyKick, energyAsymmetry, globalEnergyAsymmetry, showerClusterMapW, showerEnergyAsymmetry, singleClusterSlidingFitDataListW);

        const float energyKickScore(-energyKick / m_epsilon);
        const float energyAsymmetryScore(energyAsymmetry / m_asymmetryConstant);
        const float globalEnergyAsymmetryScore(globalEnergyAsymmetry / m_globalAsymmetryConstant);
        const float showerEnergyAsymmetryScore(showerEnergyAsymmetry / m_showerAsymmetryConstant);

        const float oldVertexScore = beamDeweighting + energyKickScore + energyAsymmetryScore + globalEnergyAsymmetryScore + showerEnergyAsymmetryScore;

        

    
        energyAsymmetry /= 3.f;
        globalEnergyAsymmetry /= 3.f;
        showerEnergyAsymmetry /= 3.f;
        energyKick /= 3.f;

        (void) vertexScoreList;

        //-------------------------------------------------------------------------------------------------------------------------------
        
        KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

        this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
        this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
        this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);
        
        bool done(false);
        
        if (m_fastScoreCheck || m_fastScoreOnly)
        {
            const float fastScore(this->GetFastScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));

            if (false)//m_fastScoreOnly)
            {
                rPhiScore = fastScore;
                done = true;
            }

            //if (fastScore < m_minFastScoreFraction * bestFastScore)
           // {
            //    rPhiScore = 0.f;
            //    done = true;
           // }

            if (fastScore > bestFastScore)
                bestFastScore = fastScore;
        }
        
        if (!m_fastScoreOnly && !done)
        {
            rPhiScore = m_fullScore ? this->GetFullScore(kernelEstimateU, kernelEstimateV, kernelEstimateW) :
                                      this->GetMidwayScore(kernelEstimateU, kernelEstimateV, kernelEstimateW);
        }
        
        TrainingSetInputs trainingSetInputs(beamDeweightingScore, rPhiScore, energyKick, energyAsymmetry, globalEnergyAsymmetry, showerEnergyAsymmetry);

        allVerticesInfoVector.emplace_back(pVertex, trainingSetInputs, 0.f, "", 0.f);
        
        
        //-------------------------------------------------------------------------------------------------------------------------------
        std::string interactionType;
        float mcVertexDr(std::numeric_limits<float>::max());
        for (const MCParticle *const pMCNeutrino : mcNeutrinoVector)
        {
            const CartesianVector mcNeutrinoPosition(pMCNeutrino->GetEndpoint());
            const float dr = (mcNeutrinoPosition - pVertex->GetPosition()).GetMagnitude();
            
            if (dr < mcVertexDr)
            {
                mcVertexDr = dr;
                
                const LArMCParticle *const pLArMCNeutrino = dynamic_cast<const LArMCParticle*>(pMCNeutrino);

                if (pLArMCNeutrino)
                    interactionType = this->ToString(this->GetInteractionType(pLArMCNeutrino, pMCParticleList, mcToGoodTrueHitListMap));
            }
        }

        //-------------------------------------------------------------------------------------------------------------------------------
        
        if (mcVertexDr != std::numeric_limits<float>::max())
        {
            vertexInfoVector.emplace_back(pVertex, trainingSetInputs, oldVertexScore, interactionType, mcVertexDr);
            if (m_trainTypeOne)
            {
                /*
                std::cout << "Beam deweighting:          " << beamDeweightingScore << std::endl;
                std::cout << "r/phi:                     " << rPhiScore << std::endl;
                std::cout << "Energy kick:               " << energyKick << std::endl;
                std::cout << "Energy asymmetry:          " << energyAsymmetry << std::endl;
                std::cout << "Global energy asym:        " << globalEnergyAsymmetry << std::endl;
                std::cout << "Shower energy asym:        " << showerEnergyAsymmetry << std::endl;
                std::cout << "Event hit showeryness:     " << eventHitShoweryness << std::endl;
                std::cout << "Event cluster showeryness: " << eventClusterShoweryness << std::endl;
                std::cout << "    => Vertex dr:          " << mcVertexDr << std::endl;
                std::cout << "Interaction type:          " << interactionType << std::endl;
                std::cout << std::endl;
                
                std::cout << m_trainingSetPrefix + interactionType + "_outputs.txt" << std::endl;
                */
            
                std::ofstream inputsFile;
                inputsFile.open(m_trainingSetPrefix + interactionType + "_inputs.txt", std::ios_base::app);
                inputsFile << beamDeweightingScore << " " << rPhiScore << " " << energyKick << " " << energyAsymmetry << " " << globalEnergyAsymmetry <<
                            " " << showerEnergyAsymmetry << " " << eventHitShoweryness << " " << eventClusterShoweryness << "\n"; 
                
                std::ofstream outputsFile;
                outputsFile.open(m_trainingSetPrefix + interactionType + "_outputs.txt", std::ios_base::app);
                outputsFile << mcVertexDr << "\n";
            }
            
        }
    }
    
    if (m_trainTypeOne)
        return;
    
    VertexInfoVector favouriteVertices;
    
    std::sort(vertexInfoVector.begin(), vertexInfoVector.end(), 
              [](const VertexInfo &lhs, const VertexInfo &rhs)
              {
                  return std::get<2>(lhs) > std::get<2>(rhs);
              });

    for (const VertexInfo &vertexInfo : vertexInfoVector)
    {
        const Vertex * const pVertex = std::get<0>(vertexInfo);
        
        bool farEnoughAway(true);
        
        for (const VertexInfo &item : favouriteVertices)
        {
            const Vertex * const pListedVertex = std::get<0>(item);
            
            if (pListedVertex == pVertex)
            {
                farEnoughAway = false;
                break;
            }
            
            const float distance = (pListedVertex->GetPosition() - pVertex->GetPosition()).GetMagnitude();
            
            if (distance <= 5.f)
            {
                farEnoughAway = false;
                break;
            }
        }
        
        if (farEnoughAway)
            favouriteVertices.push_back(vertexInfo);
            
        if (favouriteVertices.size() >= 5)
            break;
    }
        
    if (m_trainTypeTwo)
    {
        if (favouriteVertices.size() == 5)
        {
            float closestDr = std::numeric_limits<float>::max();
            const Vertex *pBestVertex(NULL);
            
            for (const VertexInfo &item : favouriteVertices)
            {
                const float dr = std::get<4>(item);
                
                if (dr < closestDr)
                {
                    closestDr = dr;
                    pBestVertex = std::get<0>(item);
                }
            }
            
            if (pBestVertex)
            {
                std::string interactionType = std::get<3>(favouriteVertices.front());
                
                std::ofstream inputsFile;
                inputsFile.open(m_trainingSetPrefix + interactionType + "_inputs.txt", std::ios_base::app);
                
                std::ofstream outputsFile;
                outputsFile.open(m_trainingSetPrefix + interactionType + "_outputs.txt", std::ios_base::app);
                
                const auto favouriteVerticesCopy = favouriteVertices;
                
                int vtxNumber(1);
                for (const VertexInfo &item : favouriteVertices)
                {
                    std::cout << "Vertex number " << std::to_string(vtxNumber++) << " : " << std::get<0>(item) << std::endl;
                    
                    /*
                    const CartesianVector position = std::get<0>(item)->GetPosition();
                    const CartesianVector positionU = LArGeometryHelper::ProjectPosition(this->GetPandora(), position, TPC_VIEW_U);
                    const CartesianVector positionV = LArGeometryHelper::ProjectPosition(this->GetPandora(), position, TPC_VIEW_V);
                    const CartesianVector positionW = LArGeometryHelper::ProjectPosition(this->GetPandora(), position, TPC_VIEW_W);
                    
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionU, "VTX " + std::to_string(vtxNumber) + " U", RED, 1);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionV, "VTX " + std::to_string(vtxNumber) + " V", RED, 1);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionW, "VTX " + std::to_string(vtxNumber) + " W", RED, 1);
                    
                    ++vtxNumber;
                    */
                    
                    int isWinningVertex = (std::get<0>(item) == pBestVertex) ? 1 : 0;
                    TrainingSetInputs trainingSetInputs = std::get<1>(item);
                    
                    std::cout << "        -> Number " << std::get<0>(trainingSetInputs) << std::endl;
                    
                   // const float mcVertexDr = std::get<4>(item);
                   
                   inputsFile << std::get<0>(eventInputs) << " " <<
                                 std::get<1>(eventInputs) << " " <<
                                 std::get<0>(trainingSetInputs) << " " <<
                                 std::get<1>(trainingSetInputs) << " " <<
                                 std::get<2>(trainingSetInputs) << " " <<
                                 std::get<3>(trainingSetInputs) << " " <<
                                 std::get<4>(trainingSetInputs) << " " <<
                                 std::get<5>(trainingSetInputs);
                            
                    
                    /*
                    std::cout << " *** THIS VERTEX *** " << std::endl;
                    std::cout << "Beam deweighting:          " << std::get<0>(trainingSetInputs) << std::endl;
                    std::cout << "r/phi:                     " << std::get<1>(trainingSetInputs) << std::endl;
                    std::cout << "Energy kick:               " << std::get<2>(trainingSetInputs) << std::endl;
                    std::cout << "Energy asymmetry:          " << std::get<3>(trainingSetInputs) << std::endl;
                    std::cout << "Global energy asym:        " << std::get<4>(trainingSetInputs) << std::endl;
                    std::cout << "Shower energy asym:        " << std::get<5>(trainingSetInputs) << std::endl;
                    std::cout << "Event hit showeryness:     " << std::get<6>(trainingSetInputs) << std::endl;
                    std::cout << "Event cluster showeryness: " << std::get<7>(trainingSetInputs) << std::endl;
                    std::cout << "    => Vertex dr:          " << mcVertexDr << std::endl;
                    std::cout << "    => Is winning vertex:  " << isWinningVertex << std::endl;
                    std::cout << "Interaction type:          " << interactionType << std::endl;
                    std::cout << std::endl;
                    
                    std::cout << m_trainingSetPrefix + interactionType + "_outputs.txt" << std::endl;
                    */
                    
                    int innerVertexNumber(1);
                    for (const VertexInfo &otherItem : favouriteVerticesCopy)
                    {
                        std::cout << "    -> Vertex number " << std::to_string(innerVertexNumber) << " : " << std::get<0>(otherItem) << std::endl;
                        
                        if (std::get<0>(otherItem) != std::get<0>(item))
                        {
                            
                        std::cout << "    -> Vertex number " << std::to_string(innerVertexNumber++) << " : " << std::get<0>(otherItem) << std::endl;
                        TrainingSetInputs otherTrainingSetInputs = std::get<1>(otherItem);
                        std::cout << "        -> Number " << std::get<0>(otherTrainingSetInputs) << std::endl;
                        
                        
                        //int otherIsWinningVertex = (std::get<0>(otherItem) == pBestVertex) ? 1 : 0;
                        //std::string otherInteractionType = std::get<3>(otherItem);
                        //const float otherMcVertexDr = std::get<4>(otherItem);
                        
                        inputsFile << " " << std::get<0>(otherTrainingSetInputs)
                                   << " " << std::get<1>(otherTrainingSetInputs)
                                   << " " << std::get<2>(otherTrainingSetInputs)
                                   << " " << std::get<3>(otherTrainingSetInputs)
                                   << " " << std::get<4>(otherTrainingSetInputs)
                                   << " " << std::get<5>(otherTrainingSetInputs);
                        
                        /*
                        std::cout << "     *** OTHER VERTEX *** " << std::endl;
                        std::cout << "    Beam deweighting:          " << std::get<0>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    r/phi:                     " << std::get<1>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Energy kick:               " << std::get<2>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Energy asymmetry:          " << std::get<3>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Global energy asym:        " << std::get<4>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Shower energy asym:        " << std::get<5>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Event hit showeryness:     " << std::get<6>(otherTrainingSetInputs) << std::endl;
                        std::cout << "    Event cluster showeryness: " << std::get<7>(otherTrainingSetInputs) << std::endl;
                        std::cout << "        => Vertex dr:          " << otherMcVertexDr << std::endl;
                        std::cout << "        => Is winning vertex:  " << otherIsWinningVertex << std::endl;
                        std::cout << "    Interaction type:          " << otherInteractionType << std::endl;
                        std::cout << std::endl;
                        
                        std::cout << "    " << m_trainingSetPrefix + otherInteractionType + "_outputs.txt" << std::endl;
                        */
                        }
                    }
                    
                    inputsFile << "\n";
                    outputsFile << isWinningVertex << "\n";
                    
                    /*
                    
                    inputsFile << beamDeweightingScore << " " << rPhiScore << " " << energyKick << " " << energyAsymmetry << " " << globalEnergyAsymmetry <<
                                " " << showerEnergyAsymmetry << " " << eventHitShoweryness << " " << eventClusterShoweryness << "\n"; 
                    
                    
                    outputsFile << mcVertexDr << "\n";
                    */
                    
                }
            }
        }
    }
    
    if (m_trainTypeTwo)
        return;
    
    if (favouriteVertices.empty())
        throw;
        
        for (const auto &item : favouriteVertices)
    {
        const Vertex * const pVertex(std::get<0>(item));
        std::cout << "Considering vertiuces " << pVertex << std::endl;
    }

    for (const auto &item : favouriteVertices)
    {
        const Vertex * const pVertex(std::get<0>(item));
        TrainingSetInputsList trainingSetInputsList{std::get<1>(item)};
        
        for (const auto &otherItem : favouriteVertices)
        {
            const Vertex * const pOtherVertex(std::get<0>(otherItem));
            
            if (pOtherVertex == pVertex)
                continue;
            
            trainingSetInputsList.push_back(std::get<1>(otherItem));
        }
        
        const double vertexScore = this->GetSVMScore(eventInputs, trainingSetInputsList);
        
        std::cout << "*** SCORE FOR THIS VERTEX : " << vertexScore << std::endl;
        
        vertexScoreList.emplace_back(pVertex, vertexScore);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

double EnergyKickVertexSelectionAlgorithm::GetSVMScore(const EventInputs &eventInputs, TrainingSetInputsList trainingSetInputsList) const
{
    if (trainingSetInputsList.empty())
        throw;
    
    std::cout << trainingSetInputsList.size() << std::endl;
    
    while (trainingSetInputsList.size() < 5)
    {
        std::cout << "Making it bigger" << std::endl;
        trainingSetInputsList.push_back(trainingSetInputsList.front());
        
    }
    
    double totalScore(0.f);
    
    for (int i = 0; i < m_supportVectors.size(); ++i)
    {
        const FeatureVector supportVector(m_supportVectors.at(i));
        
        // Produce the both normalized and scaled features.
        std::vector<double> scaledFeatures;
        
        scaledFeatures.push_back((std::get<0>(eventInputs) - std::get<0>(m_svMu)) / (std::get<0>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(eventInputs) - std::get<1>(m_svMu)) / (std::get<1>(m_svSigma) * m_svScale));
        
        const TrainingSetInputs &trainingSetInputs1 = trainingSetInputsList.at(0);
        const TrainingSetInputs &trainingSetInputs2 = trainingSetInputsList.at(1);
        const TrainingSetInputs &trainingSetInputs3 = trainingSetInputsList.at(2);
        const TrainingSetInputs &trainingSetInputs4 = trainingSetInputsList.at(3);
        const TrainingSetInputs &trainingSetInputs5 = trainingSetInputsList.at(4);
        
        scaledFeatures.push_back((std::get<0>(trainingSetInputs1) - std::get<2>(m_svMu)) / (std::get<2>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(trainingSetInputs1) - std::get<3>(m_svMu)) / (std::get<3>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<2>(trainingSetInputs1) - std::get<4>(m_svMu)) / (std::get<4>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<3>(trainingSetInputs1) - std::get<5>(m_svMu)) / (std::get<5>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<4>(trainingSetInputs1) - std::get<6>(m_svMu)) / (std::get<6>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<5>(trainingSetInputs1) - std::get<7>(m_svMu)) / (std::get<7>(m_svSigma) * m_svScale));
        
        scaledFeatures.push_back((std::get<0>(trainingSetInputs2) - std::get<8>(m_svMu)) / (std::get<8>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(trainingSetInputs2) - std::get<9>(m_svMu)) / (std::get<9>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<2>(trainingSetInputs2) - std::get<10>(m_svMu)) / (std::get<10>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<3>(trainingSetInputs2) - std::get<11>(m_svMu)) / (std::get<11>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<4>(trainingSetInputs2) - std::get<12>(m_svMu)) / (std::get<12>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<5>(trainingSetInputs2) - std::get<13>(m_svMu)) / (std::get<13>(m_svSigma) * m_svScale));
        
        scaledFeatures.push_back((std::get<0>(trainingSetInputs3) - std::get<14>(m_svMu)) / (std::get<14>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(trainingSetInputs3) - std::get<15>(m_svMu)) / (std::get<15>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<2>(trainingSetInputs3) - std::get<16>(m_svMu)) / (std::get<16>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<3>(trainingSetInputs3) - std::get<17>(m_svMu)) / (std::get<17>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<4>(trainingSetInputs3) - std::get<18>(m_svMu)) / (std::get<18>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<5>(trainingSetInputs3) - std::get<19>(m_svMu)) / (std::get<19>(m_svSigma) * m_svScale));
        
        scaledFeatures.push_back((std::get<0>(trainingSetInputs4) - std::get<20>(m_svMu)) / (std::get<20>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(trainingSetInputs4) - std::get<21>(m_svMu)) / (std::get<21>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<2>(trainingSetInputs4) - std::get<22>(m_svMu)) / (std::get<22>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<3>(trainingSetInputs4) - std::get<23>(m_svMu)) / (std::get<23>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<4>(trainingSetInputs4) - std::get<24>(m_svMu)) / (std::get<24>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<5>(trainingSetInputs4) - std::get<25>(m_svMu)) / (std::get<25>(m_svSigma) * m_svScale));
        
        scaledFeatures.push_back((std::get<0>(trainingSetInputs5) - std::get<26>(m_svMu)) / (std::get<26>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<1>(trainingSetInputs5) - std::get<27>(m_svMu)) / (std::get<27>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<2>(trainingSetInputs5) - std::get<28>(m_svMu)) / (std::get<28>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<3>(trainingSetInputs5) - std::get<29>(m_svMu)) / (std::get<29>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<4>(trainingSetInputs5) - std::get<30>(m_svMu)) / (std::get<30>(m_svSigma) * m_svScale));
        scaledFeatures.push_back((std::get<5>(trainingSetInputs5) - std::get<31>(m_svMu)) / (std::get<31>(m_svSigma) * m_svScale));
        
        // Get the bits of the support vector.
        // Produce the both normalized and scaled features.
        std::vector<double> supportVectorElements;
        
        supportVectorElements.push_back(std::get<0>(supportVector));
        supportVectorElements.push_back(std::get<1>(supportVector));
        supportVectorElements.push_back(std::get<2>(supportVector));
        supportVectorElements.push_back(std::get<3>(supportVector));
        supportVectorElements.push_back(std::get<4>(supportVector));
        supportVectorElements.push_back(std::get<5>(supportVector));
        supportVectorElements.push_back(std::get<6>(supportVector));
        supportVectorElements.push_back(std::get<7>(supportVector));
        supportVectorElements.push_back(std::get<8>(supportVector));
        supportVectorElements.push_back(std::get<9>(supportVector));
        supportVectorElements.push_back(std::get<10>(supportVector));
        supportVectorElements.push_back(std::get<11>(supportVector));
        supportVectorElements.push_back(std::get<12>(supportVector));
        supportVectorElements.push_back(std::get<13>(supportVector));
        supportVectorElements.push_back(std::get<14>(supportVector));
        supportVectorElements.push_back(std::get<15>(supportVector));
        supportVectorElements.push_back(std::get<16>(supportVector));
        supportVectorElements.push_back(std::get<17>(supportVector));
        supportVectorElements.push_back(std::get<18>(supportVector));
        supportVectorElements.push_back(std::get<19>(supportVector));
        supportVectorElements.push_back(std::get<20>(supportVector));
        supportVectorElements.push_back(std::get<21>(supportVector));
        supportVectorElements.push_back(std::get<22>(supportVector));
        supportVectorElements.push_back(std::get<23>(supportVector));
        supportVectorElements.push_back(std::get<24>(supportVector));
        supportVectorElements.push_back(std::get<25>(supportVector));
        supportVectorElements.push_back(std::get<26>(supportVector));
        supportVectorElements.push_back(std::get<27>(supportVector));
        supportVectorElements.push_back(std::get<28>(supportVector));
        supportVectorElements.push_back(std::get<29>(supportVector));
        supportVectorElements.push_back(std::get<30>(supportVector));
        supportVectorElements.push_back(std::get<31>(supportVector));
        
        totalScore += this->EvaluatePolynomialKernel(scaledFeatures, supportVectorElements) * m_modifiedAlphas.at(i);
    }
    
    return totalScore + m_svBias;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double EnergyKickVertexSelectionAlgorithm::EvaluatePolynomialKernel(std::vector<double> x, std::vector<double> y) const
{
    double total(0.f);
    
    if (x.size() != y.size())
        throw;
        
    for (int i = 0; i < x.size(); ++i)
        total += x.at(i) * y.at(i);
    
    total += 1.f;
    
    return total*total;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void EnergyKickVertexSelectionAlgorithm::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
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

bool EnergyKickVertexSelectionAlgorithm::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
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

void EnergyKickVertexSelectionAlgorithm::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
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

EnergyKickVertexSelectionAlgorithm::InteractionType EnergyKickVertexSelectionAlgorithm::GetInteractionType(const LArMCParticle *const pLArMCNeutrino, const MCParticleList *pMCParticleList, const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap) const
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
std::string EnergyKickVertexSelectionAlgorithm::ToString(const InteractionType interactionType) const
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

void EnergyKickVertexSelectionAlgorithm::SelectTrueNeutrinos(const MCParticleList *const pAllMCParticleList, MCParticleVector &selectedMCNeutrinoVector) const
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

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetEventHitShoweryness(const ClusterList &inputClusterListU, const ClusterList &inputClusterListV, const ClusterList &inputClusterListW) const
{
    unsigned int totNHits(0U);
    unsigned int totNShoweryHits(0U);
    
    for (const Cluster * const pCluster : inputClusterListU)
    {
        totNHits += pCluster->GetNCaloHits();
        
        if (this->IsClusterShowerLike(pCluster))
            totNShoweryHits += pCluster->GetNCaloHits();
    }
    
    for (const Cluster * const pCluster : inputClusterListV)
    {
        totNHits += pCluster->GetNCaloHits();
        
        if (this->IsClusterShowerLike(pCluster))
            totNShoweryHits += pCluster->GetNCaloHits();
    }
    
    for (const Cluster * const pCluster : inputClusterListW)
    {
        totNHits += pCluster->GetNCaloHits();
        
        if (this->IsClusterShowerLike(pCluster))
            totNShoweryHits += pCluster->GetNCaloHits();
    }
    
    return static_cast<float>(totNShoweryHits) / static_cast<float>(totNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetEventClusterShoweryness(const ClusterList &inputClusterListU, const ClusterList &inputClusterListV, const ClusterList &inputClusterListW) const
{
    unsigned int totNClusters(inputClusterListU.size() + inputClusterListV.size() + inputClusterListW.size());
    unsigned int totNShoweryClusters(0U);
    
    for (const Cluster * const pCluster : inputClusterListU)
    {        
        if (this->IsClusterShowerLike(pCluster))
            ++totNShoweryClusters;
    }
    
    for (const Cluster * const pCluster : inputClusterListV)
    {
        if (this->IsClusterShowerLike(pCluster))
            ++totNShoweryClusters;
    }
    
    for (const Cluster * const pCluster : inputClusterListW)
    {
        if (this->IsClusterShowerLike(pCluster))
            ++totNShoweryClusters;
    }
    
    return static_cast<float>(totNShoweryClusters) / static_cast<float>(totNClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::CalculateShowerClusterMap(const ClusterList &inputClusterList, ShowerClusterMap &showerClusterMap) const
{
    ShowerClusterList showerClusterList;
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

void EnergyKickVertexSelectionAlgorithm::CalculateClusterSlidingFits(const ClusterList &inputClusterList, SlidingFitDataList &singleClusterSlidingFitDataList) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    ClusterVector sortedClusters(inputClusterList.begin(), inputClusterList.end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster * const pCluster : sortedClusters)
    {
        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;
            
        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const int slidingFitWindow(std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * m_slidingFitWindow)));
        singleClusterSlidingFitDataList.emplace_back(pCluster, slidingFitWindow, slidingFitPitch);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::IncrementEnergyScoresForView(const CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, float &globalEnergyAsymmetry, 
    const ShowerClusterMap &showerClusterMap, float &showerEnergyAsymmetry, const SlidingFitDataList &singleClusterSlidingFitDataList) const
{    
    unsigned int totHits(0);
    bool useEnergy(true), useAsymmetry(true);
    float totEnergy(0.f), totEnergyKick(0.f), totHitKick(0.f);
    CartesianVector energyWeightedDirectionSum(0.f, 0.f, 0.f), hitWeightedDirectionSum(0.f, 0.f, 0.f);
    ClusterVector asymmetryClusters;

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

            this->IncrementEnergyKickParameters(pCluster, clusterDisplacement, clusterDirection, totEnergyKick, totEnergy, totHitKick, totHits);

            if (LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster) < m_maxAsymmetryDistance)
            {
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
    const CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits) const
{
    float impactParameter(clusterDisplacement.GetCrossProduct(clusterDirection).GetMagnitude());   
    const float displacement(clusterDisplacement.GetMagnitude());
        
    if (this->IsClusterShowerLike(pCluster))
    {       
        totEnergyKick += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy() * 
                        (m_showerCollapsingConstant * impactParameter + m_xOffset) / (displacement + m_rOffset);
                        
        totEnergy += m_showerDeweightingConstant * pCluster->GetElectromagneticEnergy();
        
        totHitKick += static_cast<float>(m_showerDeweightingConstant) * static_cast<float>(pCluster->GetNCaloHits()) * 
                      (m_showerCollapsingConstant * impactParameter + m_xOffset) / (displacement + m_rOffset);
                      
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

bool EnergyKickVertexSelectionAlgorithm::IsClusterShowerLike(const Cluster *const pCluster) const
{
    auto findIter = m_showerLikeClusterMap.find(pCluster);
    if (findIter != m_showerLikeClusterMap.end())
        return findIter->second;
    
    bool isClusterShowerLike(false);
    
    if (m_cheatTrackShowerId)  
    {          
        // V1
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));                                                                           
        isClusterShowerLike = (PHOTON == pMCParticle->GetParticleId()) || (E_MINUS == std::abs(pMCParticle->GetParticleId())); 

        /* V2
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));                                                                     
        const MCParticle *const pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle)); 
        isClusterShowerLike = (PHOTON == pPrimaryMCParticle->GetParticleId()) || (E_MINUS == std::abs(pPrimaryMCParticle->GetParticleId())); 
        */
    }                                                                                                                                                           

    else
        isClusterShowerLike = (pCluster->GetParticleId() == E_MINUS && LArClusterHelper::GetLength(pCluster) < m_minShowerSpineLength);
    
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
    const float totHitEnergy(beforeVtxHitEnergy + afterVtxHitEnergy);
    const unsigned int totHitCount(beforeVtxHitCount + afterVtxHitCount);
    
    if (useEnergyMetrics && (totHitEnergy > std::numeric_limits<float>::epsilon()))
        return std::fabs((afterVtxHitEnergy - beforeVtxHitEnergy)) / totHitEnergy;

    if (0 == totHitCount)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return std::fabs((static_cast<float>(afterVtxHitCount) - static_cast<float>(beforeVtxHitCount))) / static_cast<float>(totHitCount);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    // Is the below correct?
    histogramU.Scale(1.f/histogramU.GetCumulativeSum());
    histogramV.Scale(1.f/histogramV.GetCumulativeSum());
    histogramW.Scale(1.f/histogramW.GetCumulativeSum());

    // ATTN Need to renormalise histograms if ever want to directly compare fast and full scores
    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binContentU(histogramU.GetBinContent(xBin));
        const float binContentV(histogramV.GetBinContent(xBin));
        const float binContentW(histogramW.GetBinContent(xBin));
        figureOfMerit += binContentU * binContentU + binContentV * binContentV + binContentW * binContentW;
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    // Is the below correct?
    histogramU.Scale(1.f/histogramU.GetCumulativeSum());
    histogramV.Scale(1.f/histogramV.GetCumulativeSum());
    histogramW.Scale(1.f/histogramW.GetCumulativeSum());

    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binCenter(histogramU.GetXLow() + (static_cast<float>(xBin) + 0.5f) * histogramU.GetXBinWidth());        
        figureOfMerit += histogramU.GetBinContent(xBin) * kernelEstimateU.Sample(binCenter);
        figureOfMerit += histogramV.GetBinContent(xBin) * kernelEstimateV.Sample(binCenter);
        figureOfMerit += histogramW.GetBinContent(xBin) * kernelEstimateW.Sample(binCenter);
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    float figureOfMeritU(0.f);
    float figureOfMeritV(0.f);
    float figureOfMeritW(0.f);
    
    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
    {
        const auto sample = kernelEstimateU.Sample(contribution.first);
        figureOfMeritU += contribution.second * sample;
    }

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
    {
        const auto sample = kernelEstimateV.Sample(contribution.first);
        figureOfMeritV += contribution.second * sample;
    }

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
    {
        const auto sample = kernelEstimateW.Sample(contribution.first);
        figureOfMeritW += contribution.second * sample;
    }

    return (figureOfMeritU + figureOfMeritV + figureOfMeritW)/50000.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    for (const auto &hit : found)
    {
        const CartesianVector displacement(hit.data->GetPositionVector() - vertexPosition2D);
        const float magnitude(displacement.GetMagnitude());

        if (magnitude < std::numeric_limits<float>::epsilon())
            continue;

        float phi(this->atan2Fast(displacement.GetZ(), displacement.GetX()));
        float weight(1.f / (std::sqrt(magnitude + std::fabs(m_kappa))));

        if (m_enableFolding && (phi < 0.f))
        {
            phi += M_PI;
            weight *= -1.f;
        }

        kernelEstimate.AddContribution(phi, weight);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::atan2Fast(const float y, const float x) const
{
    const float ONE_QTR_PI(0.25f * M_PI);
    const float THR_QTR_PI(0.75f * M_PI);

    const float abs_y(std::max(std::fabs(y), std::numeric_limits<float>::epsilon()));
    const float abs_x(std::fabs(x));

    const float r((x < 0.f) ? (x + abs_y) / (abs_y + abs_x) : (abs_x - abs_y) / (abs_x + abs_y));
    const float angle(((x < 0.f) ? THR_QTR_PI : ONE_QTR_PI) + (0.1963f * r * r - 0.9817f) * r);

    return ((y < 0.f) ? -angle : angle); // negate if in quad III or IV
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyKickVertexSelectionAlgorithm::KernelEstimate::Sample(const float x) const
{
    const ContributionList &contributionList(this->GetContributionList());
    ContributionList::const_iterator lowerIter(contributionList.lower_bound(x - 3.f * m_sigma));
    ContributionList::const_iterator upperIter(contributionList.upper_bound(x + 3.f * m_sigma));

    float sample(0.f);
    const float gaussConstant(1.f / std::sqrt(2.f * M_PI * m_sigma * m_sigma));

    for (ContributionList::const_iterator iter = lowerIter; iter != upperIter; ++iter)
    {
        const float deltaSigma((x - iter->first) / m_sigma);
        const float gaussian(gaussConstant * std::exp(-0.5f * deltaSigma * deltaSigma));
        sample += iter->second * gaussian;
    }

    return sample;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EnergyKickVertexSelectionAlgorithm::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    
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

    if ((m_rOffset < std::numeric_limits<float>::epsilon()) || (m_epsilon < std::numeric_limits<float>::epsilon()))
    {
        std::cout << "EnergyKickVertexSelection: Invalid parameter(s), ROffset " << m_rOffset << ", Epsilon " << m_epsilon << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryDistance", m_maxAsymmetryDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAsymmetryCosAngle", m_minAsymmetryCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAsymmetryNClusters", m_maxAsymmetryNClusters));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerDeweightingConstant", m_showerDeweightingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerCollapsingConstant", m_showerCollapsingConstant));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerSpineLength", m_minShowerSpineLength));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerClusteringDistance", m_showerClusteringDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexClusterDistance", m_vertexClusterDistance));
      
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerClusterHits", m_minShowerClusterHits));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseShowerClusteringApproximation", m_useShowerClusteringApproximation));   
  
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CheatTrackShowerId", m_cheatTrackShowerId));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreCheck", m_fastScoreCheck));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreOnly", m_fastScoreOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullScore", m_fullScore));
 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Kappa", m_kappa));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinFastScoreFraction", m_minFastScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramNPhiBins", m_fastHistogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMin", m_fastHistogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMax", m_fastHistogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFolding", m_enableFolding));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetPrefix", m_trainingSetPrefix));
        
        
        
    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
