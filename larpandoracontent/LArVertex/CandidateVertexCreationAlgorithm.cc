/**
 *  @file   larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h"

#include <utility>

using namespace pandora;

namespace lar_content
{

CandidateVertexCreationAlgorithm::CandidateVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_slidingFitWindow(20),
    m_minClusterCaloHits(0),
    m_minClusterLengthSquared(0.f * 0.f),
    m_maxClusterXDiscrepancy(4.f),
    m_chiSquaredCut(2.f),
    m_enableCrossingCandidates(false),
    m_enableEnergyCandidates(false),
    m_minCrossingClusterSize(10),
    m_extrapolationLength(20.0f),
    m_extrapolationStepSize(0.1f),
    m_minClusterCrossingApproach(0.25f),
    m_postCrossingSkipDistance(5.0f),
    m_minEnergyVertexClusterSize(60),
    m_oneBinDistanceFractionalDeviationThreshold(0.15),
    m_twoBinDistanceFractionalDeviationThreshold(0.5),
    m_threeBinDistanceFractionalDeviationThreshold(0.7),
    m_minAverageBinEnergy(0.25),
    m_maxAverageBinEnergy(4.3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::Run()
{
    try
    {
        ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
        this->SelectClusters(clusterVectorU, clusterVectorV, clusterVectorW);

        const VertexList *pVertexList(NULL); std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        this->CreateClusterEndPointComparisonVertices(clusterVectorU, clusterVectorV);
        this->CreateClusterEndPointComparisonVertices(clusterVectorU, clusterVectorW);
        this->CreateClusterEndPointComparisonVertices(clusterVectorV, clusterVectorW);
        
        if (m_enableCrossingCandidates)
            this->CreateCrossingVertices(clusterVectorU, clusterVectorV, clusterVectorW);
            
        if (m_enableEnergyCandidates)
        {
            this->CreateEnergyVertices(clusterVectorU, clusterVectorV, clusterVectorW);
            this->CreateEnergyVertices(clusterVectorV, clusterVectorW, clusterVectorU);
            this->CreateEnergyVertices(clusterVectorW, clusterVectorU, clusterVectorV);
        }

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::SelectClusters(ClusterVector &clusterVectorU, ClusterVector &clusterVectorV, ClusterVector &clusterVectorW)
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterVector &selectedClusterVector((TPC_VIEW_U == hitType) ? clusterVectorU : (TPC_VIEW_V == hitType) ? clusterVectorV : clusterVectorW);

        if (!selectedClusterVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pCluster : sortedClusters)
        {
            if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
                continue;

            try
            {
                this->AddToSlidingFitCache(pCluster);
                selectedClusterVector.push_back(pCluster);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateClusterEndPointComparisonVertices(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2) const
{
    for (const Cluster *const pCluster1 : clusterVector1)
    {
        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));

        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster1));
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());

        for (const Cluster *const pCluster2 : clusterVector2)
        {
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const TwoDSlidingFitResult &fitResult2(this->GetCachedSlidingFitResult(pCluster2));
            const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
            const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());

            this->CreateVertex(maxLayerPosition1, hitType1, fitResult2);
            this->CreateVertex(minLayerPosition1, hitType1, fitResult2);
            this->CreateVertex(maxLayerPosition2, hitType2, fitResult1);
            this->CreateVertex(minLayerPosition2, hitType2, fitResult1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateCrossingVertices(const ClusterVector &clusterVectorU, const ClusterVector &clusterVectorV,
    const ClusterVector &clusterVectorW) const 
{
    CartesianPointList crossingsU;
    CartesianPointList crossingsV;
    CartesianPointList crossingsW;

    this->Find2DClusterCrossings(clusterVectorU, crossingsU);
    this->Find2DClusterCrossings(clusterVectorV, crossingsV);
    this->Find2DClusterCrossings(clusterVectorW, crossingsW);

    this->CreateMatchedVertices(crossingsU, crossingsV, TPC_VIEW_U, TPC_VIEW_V);
    this->CreateMatchedVertices(crossingsU, crossingsW, TPC_VIEW_U, TPC_VIEW_W);
    this->CreateMatchedVertices(crossingsV, crossingsW, TPC_VIEW_V, TPC_VIEW_W);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEnergyVertices(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2, const ClusterVector &clusterVector3) const
{
    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iter1End = clusterVector1.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster = *iter1;
        
        if (pCluster->GetNCaloHits() < m_minEnergyVertexClusterSize)
            continue;

        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);

        CartesianPointList energyAlongRLvector;
        this->CreateEnergyAlongRLVector(pCluster, energyAlongRLvector);

        CartesianPointList filteredEnergyAlongRLvector;
        this->FilterEnergyVector(energyAlongRLvector, filteredEnergyAlongRLvector);

        FloatVector energySpikeRLvector;
        this->FindEnergySpike(filteredEnergyAlongRLvector, energySpikeRLvector);

        this->CreateVerticesFromSpikes(energySpikeRLvector, filteredEnergyAlongRLvector, caloHitList, clusterVector1, clusterVector2, clusterVector3);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateVertex(const CartesianVector &position1, const HitType hitType1, const TwoDSlidingFitResult &fitResult2) const
{
    const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());
    
    if ((((position1.GetX() < minLayerPosition2.GetX()) && (position1.GetX() < maxLayerPosition2.GetX())) ||
        ((position1.GetX() > minLayerPosition2.GetX()) && (position1.GetX() > maxLayerPosition2.GetX()))) &&
        (std::fabs(position1.GetX() - minLayerPosition2.GetX()) > m_maxClusterXDiscrepancy) &&
        (std::fabs(position1.GetX() - maxLayerPosition2.GetX()) > m_maxClusterXDiscrepancy))
    {
        return;
    }

    CartesianVector position2(0.f, 0.f, 0.f);
    if (STATUS_CODE_SUCCESS != fitResult2.GetExtrapolatedPositionAtX(position1.GetX(), position2))
        return;

    const HitType hitType2(LArClusterHelper::GetClusterHitType(fitResult2.GetCluster()));

    float chiSquared(0.f);
    CartesianVector position3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);
    const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_W));

    if (chiSquared > m_chiSquaredCut)
        return;

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = position3D;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::Find2DClusterCrossings(const ClusterVector &clusterVector, CartesianPointList &crossingsVector) const
{
    for (ClusterVector::const_iterator iter1 = clusterVector.begin(), iter1End = clusterVector.end(); iter1 != iter1End; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;

        if (pCluster1->GetNCaloHits() < m_minCrossingClusterSize)
            continue;

        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(*iter1));

        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());

        CartesianVector minZHitPosition(0.f, 0.f, 0.f), maxZHitPosition(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(pCluster1, minZHitPosition, maxZHitPosition);

        crossingsVector.push_back(minLayerPosition1);
        crossingsVector.push_back(maxLayerPosition1);

        crossingsVector.push_back(minZHitPosition);
        crossingsVector.push_back(maxZHitPosition);

        CartesianPointList spacePointVector1;
        this->GetExtrapolatedClusterSpacepoints(pCluster1, spacePointVector1);

        for (ClusterVector::const_iterator iter2 = clusterVector.begin(), iter2End = clusterVector.end(); iter2 != iter2End; ++iter2)
        {
            const Cluster *const pCluster2 = *iter2;
            
            if ((pCluster1 == pCluster2) || pCluster2->GetNCaloHits() < m_minCrossingClusterSize)
                continue;

            CartesianPointList spacePointVector2;
            this->GetExtrapolatedClusterSpacepoints(pCluster2, spacePointVector2);

            FloatVector distancesBetweenClusters;

            for (CartesianVector &position1: spacePointVector1)
            {
                for (CartesianVector &position2: spacePointVector2)
                    distancesBetweenClusters.push_back((position1-position2).GetMagnitude());
            }

            std::sort(spacePointVector1.begin(), spacePointVector1.end(), SortSpacePointsByZ);
            std::sort(spacePointVector2.begin(), spacePointVector2.end(), SortSpacePointsByZ);

            this->FindCrossingsFromSpacepoints(spacePointVector1, spacePointVector2, crossingsVector);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::GetExtrapolatedClusterSpacepoints(const Cluster *const pCluster, CartesianPointList &spacePointVector) const
{
    CaloHitList caloHitList1;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList1);

    for (const CaloHit *const pCaloHit : caloHitList1)
    {
        CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        spacePointVector.push_back(caloHitPosition);
    }

    const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster));

    float minLayerRL(fitResult1.GetL(fitResult1.GetMinLayer()));
    float maxLayerRL(fitResult1.GetL(fitResult1.GetMaxLayer()));

    for (float i = m_extrapolationStepSize; i < m_extrapolationLength; i += m_extrapolationStepSize)
    {
        CartesianVector tempExtrapolatedPositionUnder(0.f, 0.f, 0.f);
        CartesianVector tempExtrapolatedPositionOver(0.f, 0.f, 0.f);

        fitResult1.GetExtrapolatedPosition(minLayerRL - i, tempExtrapolatedPositionUnder);
        fitResult1.GetExtrapolatedPosition(maxLayerRL + i, tempExtrapolatedPositionOver);

        spacePointVector.push_back(tempExtrapolatedPositionUnder);
        spacePointVector.push_back(tempExtrapolatedPositionOver);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindCrossingsFromSpacepoints(const CartesianPointList &spacePointVector1, const CartesianPointList &spacePointVector2,
    CartesianPointList &crossingsVector) const
{
    int skipCounter(0);
    bool shouldSkip(false);

    for (const CartesianVector &position1: spacePointVector1)
    {
        if (shouldSkip)
        {
            if (skipCounter <= (m_postCrossingSkipDistance/m_extrapolationStepSize))
            {
                skipCounter++;
                continue;
            }
        }

        skipCounter = 0;
        shouldSkip = false;

        for (const CartesianVector &position2: spacePointVector2)
        {
            if ((position1-position2).GetMagnitude() < m_minClusterCrossingApproach)
            {
                crossingsVector.push_back(position1);
                crossingsVector.push_back(position2);
                
                shouldSkip = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateMatchedVertices(CartesianPointList &crossingsVector1, CartesianPointList &crossingsVector2, HitType hitType1, HitType hitType2) const
{
    for (CartesianVector &position1: crossingsVector1)
    {
        std::vector<std::pair<CartesianVector*,float>> matched3DPositions;
        FloatVector chiSquareds;
        
        for (CartesianVector &position2: crossingsVector2)
        {
            if (std::fabs(position1.GetX() - position2.GetX()) > 0.1)
                continue;

            float chiSquared(0.f);
            CartesianVector position3D(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);
            const CartesianVector vertexProjectionU(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_U));
            const CartesianVector vertexProjectionV(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_V));
            const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, TPC_VIEW_W));

            if (chiSquared > 2.0)
                return;

            std::pair<CartesianVector*,float> positionChiSquaredPair;
            positionChiSquaredPair = std::make_pair(&position3D, chiSquared);
            
            matched3DPositions.push_back(positionChiSquaredPair);
            chiSquareds.push_back(chiSquared);
        }

        for (std::pair<CartesianVector*,float> &pair : matched3DPositions)
        {
            if (pair.second == (*std::min_element(chiSquareds.begin(), chiSquareds.end())))
            {
                PandoraContentApi::Vertex::Parameters parameters;
                parameters.m_position = *(pair.first);
                parameters.m_vertexLabel = VERTEX_INTERACTION;
                parameters.m_vertexType = VERTEX_3D;

                const Vertex *pVertex(NULL);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEnergyAlongRLVector(const Cluster *const pCluster, CartesianPointList &energyAlongRLvector) const
{
    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    CaloHitList caloHitList;
    orderedCaloHitList.GetCaloHitList(caloHitList);

    const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(pCluster));

    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());

        float caloHitEnergy(pCaloHit->GetElectromagneticEnergy());
        float rL(0.f), rT(0.f);

        slidingFitResult.GetLocalPosition(caloHitPosition, rL, rT);
        CartesianVector energyAlongRL(rL, 0.f, 1000*caloHitEnergy); //units of MeV

        energyAlongRLvector.push_back(energyAlongRL);
    }

    std::sort(energyAlongRLvector.begin(), energyAlongRLvector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FilterEnergyVector(const CartesianPointList &unfilteredEnergyVector, CartesianPointList &filteredEnergyVector) const
{
    for (CartesianPointList::const_iterator pairIter = std::next(unfilteredEnergyVector.begin(), 1), pairIterEnd = std::prev(unfilteredEnergyVector.end(), 1); pairIter != pairIterEnd; ++pairIter)
    {
        float chargeRatioNext(((*(std::next(pairIter, 1))).GetZ())/((*pairIter).GetZ()));
        float chargeRatioPrevious(((*(std::prev(pairIter, 1))).GetZ())/((*pairIter).GetZ()));

        if (!(chargeRatioPrevious < 0.8 && chargeRatioNext < 0.8))
            filteredEnergyVector.push_back(*pairIter);
    }

    std::sort(filteredEnergyVector.begin(), filteredEnergyVector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::BinEnergyRLVector(const CartesianPointList &energyAlongRLvector, CartesianPointList &binnedEnergyAlongRLvector) const
{
    int clusterLength(100);
    for (float i = 0.5; i != clusterLength; i += 0.5)
    {
        int nHitsBin(0);
        float averageBinPosition(0.f);
        float averageBinEnergy(0.f);

        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (energyRL.GetX() > i)
                break;

            if (energyRL.GetX() < i && energyRL.GetX() > (i - 1))
            {
                nHitsBin++;
                averageBinPosition += energyRL.GetX();
                averageBinEnergy += energyRL.GetZ();
            }
        }

        if (nHitsBin == 0)
            continue;

        averageBinPosition /= nHitsBin;
        averageBinEnergy /= nHitsBin;

        CartesianVector bin(averageBinPosition, 0.f, averageBinEnergy);
        binnedEnergyAlongRLvector.push_back(bin);
    }

    std::sort(binnedEnergyAlongRLvector.begin(), binnedEnergyAlongRLvector.end(), SortEnergyVectorByRL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindEnergySpike(const CartesianPointList &energyAlongRLvector, FloatVector &spikeRLvector) const
{
    CartesianPointList binnedEnergyAlongRLvector;
    this->BinEnergyRLVector(energyAlongRLvector, binnedEnergyAlongRLvector);

    FloatVector binRLvector;
    this->FindBinWithSpike(binnedEnergyAlongRLvector, binRLvector);

    this->ConvertBinRLToSpikeRL(binRLvector, spikeRLvector, energyAlongRLvector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindBinWithSpike(const CartesianPointList &binnedEnergyAlongRLvector, FloatVector &binRLvector) const
{
    float averageEnergyDifference(0.f);

    for (CartesianPointList::const_iterator pairIter = std::next(binnedEnergyAlongRLvector.begin(), 1), pairIterEnd = std::prev(binnedEnergyAlongRLvector.end(), 1); pairIter != pairIterEnd; ++pairIter)
    {
        float thisBinAverageEnergy((*pairIter).GetZ());
        float nextBinAverageEnergy((*std::next(pairIter, 1)).GetZ());
        float energyDifference(std::fabs(nextBinAverageEnergy / thisBinAverageEnergy));
        averageEnergyDifference += (energyDifference);
    }

    averageEnergyDifference /= binnedEnergyAlongRLvector.size();

    if (averageEnergyDifference > 1.5)
        return;

    static float bestScore(1.5f);

    for (CartesianPointList::const_iterator pairIter = std::next(binnedEnergyAlongRLvector.begin(), 1), pairIterEnd = std::prev(binnedEnergyAlongRLvector.end(), 2); pairIter != pairIterEnd; ++pairIter)
    {
        float thisBinAveragePosition((*pairIter).GetX());
        float thisBinAverageEnergy((*pairIter).GetZ());

        float previousBinAverageEnergy((*std::prev(pairIter, 1)).GetZ());
        float nextBinAverageEnergy((*std::next(pairIter, 1)).GetZ());

        float previousPreviousBinAverageEnergy((*std::prev(pairIter, 2)).GetZ());
        float nextNextBinAverageEnergy((*std::next(pairIter, 2)).GetZ());

        float previousPreviousPreviousBinAverageEnergy((*std::prev(pairIter, 3)).GetZ());
        float nextNextNextBinAverageEnergy((*std::next(pairIter, 3)).GetZ());

        float previousPreviousPreviousPreviousBinAverageEnergy((*std::prev(pairIter, 4)).GetZ());
        float nextNextNextNextBinAverageEnergy((*std::next(pairIter, 4)).GetZ());

        float nextBinFractionalDeviation(std::fabs(1 - std::fabs(nextBinAverageEnergy / thisBinAverageEnergy)));
        float previousBinFractionalDeviation(std::fabs(1 - std::fabs(previousBinAverageEnergy / thisBinAverageEnergy)));

        float nextNextBinFractionalDeviation(std::fabs(1 - (std::fabs(nextNextBinAverageEnergy / thisBinAverageEnergy))));
        float previousPreviousBinFractionalDeviation(std::fabs(1 - std::fabs(previousPreviousBinAverageEnergy / thisBinAverageEnergy)));

        float nextNextNextBinFractionalDeviation(std::fabs(1 - (std::fabs(nextNextNextBinAverageEnergy / thisBinAverageEnergy))));
        float previousPreviousPreviousBinFractionalDeviation(std::fabs(1 - std::fabs(previousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy)));

        if (((nextBinFractionalDeviation > m_oneBinDistanceFractionalDeviationThreshold) && (nextNextBinFractionalDeviation > m_twoBinDistanceFractionalDeviationThreshold) && (nextNextNextBinFractionalDeviation > m_threeBinDistanceFractionalDeviationThreshold)
        && (previousBinFractionalDeviation < m_oneBinDistanceFractionalDeviationThreshold) && (previousPreviousBinFractionalDeviation < m_twoBinDistanceFractionalDeviationThreshold))
        || ((previousBinFractionalDeviation > m_oneBinDistanceFractionalDeviationThreshold) && (previousPreviousBinFractionalDeviation > m_twoBinDistanceFractionalDeviationThreshold) && (previousPreviousPreviousBinFractionalDeviation > m_threeBinDistanceFractionalDeviationThreshold)
        && (nextBinFractionalDeviation < m_oneBinDistanceFractionalDeviationThreshold) && (nextNextBinFractionalDeviation < m_twoBinDistanceFractionalDeviationThreshold)))
        {
            if (thisBinAverageEnergy < m_minAverageBinEnergy || thisBinAverageEnergy > m_maxAverageBinEnergy) //this is because currently the energy is capped at a certain level
                continue;

            float score((nextBinAverageEnergy / thisBinAverageEnergy + nextNextBinAverageEnergy / thisBinAverageEnergy + nextNextNextBinAverageEnergy / thisBinAverageEnergy + nextNextNextNextBinAverageEnergy / thisBinAverageEnergy)
                / (previousBinAverageEnergy / thisBinAverageEnergy + previousPreviousBinAverageEnergy / thisBinAverageEnergy));

            float scoreTwo((previousBinAverageEnergy / thisBinAverageEnergy + previousPreviousBinAverageEnergy / thisBinAverageEnergy + previousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy + previousPreviousPreviousPreviousBinAverageEnergy / thisBinAverageEnergy)
                / (nextBinAverageEnergy / thisBinAverageEnergy + nextNextBinAverageEnergy / thisBinAverageEnergy));

            float workingScore(0.f);

            if (score > scoreTwo)
                workingScore = score;
            else
                workingScore = scoreTwo;

            if (workingScore > bestScore)
            {
                binRLvector.clear();
                binRLvector.push_back(thisBinAveragePosition);
                bestScore = workingScore;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ConvertBinRLToSpikeRL(const FloatVector &binRLvector, FloatVector &spikeRLvector, const CartesianPointList &energyAlongRLvector) const
{
    for (const float &binRL : binRLvector)
    {
        float closestApproach(100.f);
        float closestMatchingRL(0.f);

        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (std::fabs(energyRL.GetX() - binRL) < closestApproach)
            {
                closestApproach = std::fabs(energyRL.GetX() - binRL);
                closestMatchingRL = energyRL.GetX();
            }
        }

        for (const CartesianVector &energyRL : energyAlongRLvector)
        {
            if (energyRL.GetX() == closestMatchingRL)
            {
                spikeRLvector.push_back(energyRL.GetX());
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::ConvertRLtoCaloHit(const float spikeRL, const CartesianPointList &filteredEnergyAlongRLvector,
    const CaloHitList &caloHitList, CartesianVector &hitPosition) const
{
    std::vector<std::pair<float, float>> XtoRLvector;

    float targetEnergy(0.f);

    for (const CartesianVector &energyRL : filteredEnergyAlongRLvector)
    {
        if (energyRL.GetX() == spikeRL) // TODO ... ???
            targetEnergy = energyRL.GetZ();
    }

    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit(*hitIter);
        const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());
        const float caloHitEnergy(1000* pCaloHit->GetElectromagneticEnergy()); //MeV
        
        if (caloHitEnergy == targetEnergy) // TODO ... ???
            hitPosition = caloHitPosition;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindMatchingHitsInDifferentView(const ClusterVector &clusterVector, const CartesianVector &energySpikePosition,
    CartesianPointList &matchedHits) const
{
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iter1End = clusterVector.end(); iter != iter1End; ++iter)
    {
        const Cluster *const pCluster = *iter;

        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList caloHitList;
        orderedCaloHitList.GetCaloHitList(caloHitList);

        for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit(*hitIter);
            const CartesianVector caloHitPosition(pCaloHit->GetPositionVector());

            if (caloHitPosition.GetX() + 0.25 > energySpikePosition.GetX() && caloHitPosition.GetX() - 0.25 < energySpikePosition.GetX())
                matchedHits.push_back(caloHitPosition);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateVerticesFromSpikes(const FloatVector &energySpikeRLvector, const CartesianPointList &filteredEnergyAlongRLvector,
    const CaloHitList &caloHitList, const ClusterVector &clusterVector1, const ClusterVector &clusterVector2, const ClusterVector &clusterVector3) const
{
    for (const float energySpikeRL : energySpikeRLvector)
    {
        CartesianVector energySpikePosition(0.f, 0.f, 0.f);
        this->ConvertRLtoCaloHit(energySpikeRL, filteredEnergyAlongRLvector, caloHitList, energySpikePosition);

        CartesianPointList energySpikes1, matchedHits2, matchedHits3;

        energySpikes1.push_back(energySpikePosition);

        this->FindMatchingHitsInDifferentView(clusterVector2, energySpikePosition, matchedHits2);
        this->FindMatchingHitsInDifferentView(clusterVector3, energySpikePosition, matchedHits3);

        if (!clusterVector1.empty() && !clusterVector2.empty())
        {
            const Cluster *const pCluster1(clusterVector1.at(0));
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));

            const Cluster *const pCluster2(clusterVector2.at(0));
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            this->CreateMatchedVertices(energySpikes1, matchedHits2, hitType1, hitType2);
        }

        if (!clusterVector1.empty() && !clusterVector3.empty())
        {
            const Cluster *const pCluster1(clusterVector1.at(0));
            const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
            
            const Cluster *const pCluster3(clusterVector3.at(0));
            const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));
            
            this->CreateMatchedVertices(energySpikes1, matchedHits3, hitType1, hitType3);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &CandidateVertexCreationAlgorithm::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CandidateVertexCreationAlgorithm::SortSpacePointsByZ(CartesianVector &vector1, CartesianVector &vector2)
{
    return vector1.GetZ() < vector2.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CandidateVertexCreationAlgorithm::SortEnergyVectorByRL(CartesianVector &vector1, CartesianVector &vector2)
{
    return vector1.GetX() < vector2.GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterXDiscrepancy", m_maxClusterXDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableCrossingCandidates", m_enableCrossingCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableEnergyCandidates", m_enableEnergyCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCrossingClusterSize", m_minCrossingClusterSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationLength", m_extrapolationLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExtrapolationStepSize", m_extrapolationStepSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCrossingApproach", m_minClusterCrossingApproach));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostCrossingSkipDistance", m_postCrossingSkipDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinEnergyVertexClusterSize", m_minEnergyVertexClusterSize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content