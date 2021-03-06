/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LongitudinalExtensionAlgorithm::LongitudinalExtensionAlgorithm() :
    m_clusterMinLength(5.f),
    m_clusterMinLayerOccupancy(0.75f),
    m_nodeMaxDisplacement(1.5f),
    m_nodeMaxCosRelativeAngle(0.906f),
    m_emissionMaxLongitudinalDisplacement(15.f),
    m_emissionMaxTransverseDisplacement(2.5f),
    m_emissionMaxCosRelativeAngle(0.985f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        if (LArClusterHelper::GetLayerOccupancy(pCluster) < m_clusterMinLayerOccupancy)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Convert each input cluster into a pointing cluster
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        try
        {
            pointingClusterList.push_back(LArPointingCluster(*iter));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Form associations between pairs of pointing clusters
    for (LArPointingClusterList::const_iterator iterI = pointingClusterList.begin(), iterEndI = pointingClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster &clusterI = *iterI;

        for (LArPointingClusterList::const_iterator iterJ = iterI, iterEndJ = pointingClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster &clusterJ = *iterJ;

            if (clusterI.GetCluster() == clusterJ.GetCluster())
                continue;

            this->FillClusterAssociationMatrix(clusterI, clusterJ, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillClusterAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pClusterI(clusterI.GetCluster());
    const Cluster *const pClusterJ(clusterJ.GetCluster());

    if (pClusterI == pClusterJ)
        return;

    // Check that new layer occupancy would be reasonable
    if (LArClusterHelper::GetLayerOccupancy(pClusterI, pClusterJ) < m_clusterMinLayerOccupancy)
        return;

    // Identify closest pair of vertices 
    LArPointingCluster::Vertex targetVertexI, targetVertexJ;

    try
    {
        LArPointingClusterHelper::GetClosestVertices(clusterI, clusterJ, targetVertexI, targetVertexJ);
    }
    catch (StatusCodeException &)
    {
        return;
    }

    // (Just in case...)
    if (!(targetVertexI.IsInitialized() && targetVertexJ.IsInitialized()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector &vertexPositionI(targetVertexI.GetPosition());
    const CartesianVector &vertexPositionJ(targetVertexJ.GetPosition());
    const CartesianVector &vertexDirectionI(targetVertexI.GetDirection());
    const CartesianVector &vertexDirectionJ(targetVertexJ.GetDirection());

    // Check for reasonable proximity between vertices
    const float distanceSquared((vertexPositionI - vertexPositionJ).GetMagnitudeSquared());

    if (distanceSquared > m_emissionMaxLongitudinalDisplacement * m_emissionMaxLongitudinalDisplacement)
        return;

    // Check that vertices have a reasonable linear fit
    if (targetVertexI.GetRms() > 1.f || targetVertexJ.GetRms() > 1.f)
        return;

    // Association type
    ClusterAssociation::AssociationType associationType(ClusterAssociation::NONE);

    // Requirements for Nodes
    if (distanceSquared < 2.f * m_nodeMaxDisplacement * m_nodeMaxDisplacement)
    {
        associationType = ClusterAssociation::WEAK;

        if (distanceSquared < m_nodeMaxDisplacement * m_nodeMaxDisplacement)
        {
            const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));

            if (cosTheta > m_nodeMaxCosRelativeAngle)
            {
                associationType = ClusterAssociation::STRONG;
            }
        }
    }

    // Requirements for Emissions
    const float clusterLengthI(LArPointingClusterHelper::GetLength(clusterI));
    const float clusterLengthJ(LArPointingClusterHelper::GetLength(clusterJ));

    if (associationType < ClusterAssociation::STRONG)
    {
        const float cosTheta(-vertexDirectionI.GetDotProduct(vertexDirectionJ));
        const float cosThetaI((vertexPositionI - vertexPositionJ).GetUnitVector().GetDotProduct(vertexDirectionI));
        const float cosThetaJ((vertexPositionJ - vertexPositionI).GetUnitVector().GetDotProduct(vertexDirectionJ));

        float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
        LArPointingClusterHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL1, rT1);
        LArPointingClusterHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, vertexPositionI, rL2, rT2);

        if ((rL1 > -2.5f && rL1 < std::min(0.66f * clusterLengthJ, m_emissionMaxLongitudinalDisplacement)) &&
            (rL2 > -2.5f && rL2 < std::min(0.66f * clusterLengthI, m_emissionMaxLongitudinalDisplacement)) &&
            (rT1 < m_emissionMaxTransverseDisplacement) && (rT2 < m_emissionMaxTransverseDisplacement) &&
            (targetVertexI.GetRms() < 0.5f && targetVertexJ.GetRms() < 0.5f) &&
            (cosTheta > m_emissionMaxCosRelativeAngle) && (std::fabs(cosThetaI) > 0.25f) && (std::fabs(cosThetaJ) > 0.25f))
        {
            associationType = ClusterAssociation::STRONG;
        }
    }

    // Record the association
    if (ClusterAssociation::NONE != associationType)
    {
        const ClusterAssociation::VertexType vertexTypeI(targetVertexI.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        const ClusterAssociation::VertexType vertexTypeJ(targetVertexJ.IsInnerVertex() ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        (void) clusterAssociationMatrix[pClusterI].insert(ClusterAssociationMap::value_type(pClusterJ, ClusterAssociation(vertexTypeI, vertexTypeJ, associationType, clusterLengthJ)));
        (void) clusterAssociationMatrix[pClusterJ].insert(ClusterAssociationMap::value_type(pClusterI, ClusterAssociation(vertexTypeJ, vertexTypeI, associationType, clusterLengthI)));
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &inputAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{
    // Decide which associations will become merges
    // To make the merge A <-> B, both A -> B and B -> A must be strong associations
    // with the largest figures of merit of all the A -> X and B -> Y associations

    // First step: remove double-counting from the map of associations
    // i.e. if the map has A <-> B, B <-> C, A <-> C, then remove A <-> C
    ClusterAssociationMatrix clusterAssociationMatrix;

    ClusterVector sortedInputClusters;
    for (const auto &mapEntry : inputAssociationMatrix) sortedInputClusters.push_back(mapEntry.first);
    std::sort(sortedInputClusters.begin(), sortedInputClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : sortedInputClusters)
    {
        const ClusterAssociationMap &associationMap1(inputAssociationMatrix.at(pCluster1));

        for (const Cluster *const pCluster2 : sortedInputClusters)
        {
            if (pCluster1 == pCluster2)
                continue;

            const ClusterAssociationMap &associationMap2(inputAssociationMatrix.at(pCluster2));

            ClusterAssociationMap::const_iterator iter12 = associationMap1.find(pCluster2);
            if (associationMap1.end() == iter12)
                continue;

            ClusterAssociationMap::const_iterator iter21 = associationMap2.find(pCluster1);
            if (associationMap2.end() == iter21)
                continue;

            const ClusterAssociation &association12(iter12->second);
            const ClusterAssociation &association21(iter21->second);

            bool isAssociated(true);

            ClusterVector sortedAssociationClusters;
            for (const auto &mapEntry : associationMap1) sortedAssociationClusters.push_back(mapEntry.first);
            std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

            for (const Cluster *const pCluster3 : sortedAssociationClusters)
            {
                const ClusterAssociation &association13(associationMap1.at(pCluster3));

                ClusterAssociationMap::const_iterator iter23 = associationMap2.find(pCluster3);
                if (associationMap2.end() == iter23)
                    continue;

                const ClusterAssociation &association23(iter23->second);

                if (association12.GetParent() == association13.GetParent() &&
                    association23.GetParent() == association21.GetParent() &&
                    association13.GetDaughter() != association23.GetDaughter())
                {
                    isAssociated = false;
                    break;
                }
            }

            if (isAssociated)
            {
                (void) clusterAssociationMatrix[pCluster1].insert(ClusterAssociationMap::value_type(pCluster2, association12));
                (void) clusterAssociationMatrix[pCluster2].insert(ClusterAssociationMap::value_type(pCluster1, association21));
            }
        }
    }


    // Second step: find the best associations A -> X and B -> Y
    ClusterAssociationMatrix intermediateAssociationMatrix;

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clusterAssociationMatrix) sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedClusters)
    {
        const ClusterAssociationMap &clusterAssociationMap(clusterAssociationMatrix.at(pParentCluster));

        const Cluster *pBestClusterInner(nullptr);
        ClusterAssociation bestAssociationInner(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        const Cluster *pBestClusterOuter(nullptr);
        ClusterAssociation bestAssociationOuter(ClusterAssociation::UNDEFINED, ClusterAssociation::UNDEFINED, ClusterAssociation::NONE, 0.f);

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : clusterAssociationMap) sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &clusterAssociation(clusterAssociationMap.at(pDaughterCluster));

            // Inner associations
            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationInner.GetFigureOfMerit())
                {
                    bestAssociationInner = clusterAssociation;

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                        pBestClusterInner = pDaughterCluster;
                    else
                        pBestClusterInner = nullptr;
                }
            }

            // Outer associations
            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() > bestAssociationOuter.GetFigureOfMerit())
                {
                    bestAssociationOuter = clusterAssociation;

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                        pBestClusterOuter = pDaughterCluster;
                    else
                        pBestClusterOuter = nullptr;
                }
            }
        }

        if (pBestClusterInner)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterInner, bestAssociationInner));

        if (pBestClusterOuter)
            (void) intermediateAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pBestClusterOuter, bestAssociationOuter));
    }


    // Third step: make the merge if A -> X and B -> Y is in fact A -> B and B -> A
    ClusterVector intermediateSortedClusters;
    for (const auto &mapEntry : intermediateAssociationMatrix) intermediateSortedClusters.push_back(mapEntry.first);
    std::sort(intermediateSortedClusters.begin(), intermediateSortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : intermediateSortedClusters)
    {
        const ClusterAssociationMap &parentAssociationMap(intermediateAssociationMatrix.at(pParentCluster));

        ClusterVector sortedAssociationClusters;
        for (const auto &mapEntry : parentAssociationMap) sortedAssociationClusters.push_back(mapEntry.first);
        std::sort(sortedAssociationClusters.begin(), sortedAssociationClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedAssociationClusters)
        {
            const ClusterAssociation &parentToDaughterAssociation(parentAssociationMap.at(pDaughterCluster));

            ClusterAssociationMatrix::const_iterator iter5 = intermediateAssociationMatrix.find(pDaughterCluster);

            if (intermediateAssociationMatrix.end() == iter5)
                continue;

            const ClusterAssociationMap &daughterAssociationMap(iter5->second);

            ClusterAssociationMap::const_iterator iter6 = daughterAssociationMap.find(pParentCluster);

            if (daughterAssociationMap.end() == iter6)
                continue;

            const ClusterAssociation &daughterToParentAssociation(iter6->second);

            if (parentToDaughterAssociation.GetParent() == daughterToParentAssociation.GetDaughter() &&
                parentToDaughterAssociation.GetDaughter() == daughterToParentAssociation.GetParent())
            {
                ClusterList &parentList(clusterMergeMap[pParentCluster]);

                if (parentList.end() == std::find(parentList.begin(), parentList.end(), pDaughterCluster))
                    parentList.push_back(pDaughterCluster);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLayerOccupancy", m_clusterMinLayerOccupancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NodeMaxDisplacement", m_nodeMaxDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NodeMaxCosRelativeAngle", m_nodeMaxCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EmissionMaxLongitudinalDisplacement", m_emissionMaxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EmissionMaxTransverseDisplacement", m_emissionMaxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EmissionMaxCosRelativeAngle", m_emissionMaxCosRelativeAngle));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
