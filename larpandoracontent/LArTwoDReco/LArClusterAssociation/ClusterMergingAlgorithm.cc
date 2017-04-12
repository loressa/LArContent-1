/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ClusterMergingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;

    if (m_inputClusterListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));
    }

    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << m_inputClusterListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    while (true)
    {
        ClusterVector unsortedVector, clusterVector;
        this->GetListOfCleanClusters(pClusterList, unsortedVector);
        this->GetSortedListOfCleanClusters(unsortedVector, clusterVector);

        ClusterMergeMap clusterMergeMap;
        this->PopulateClusterMergeMap(clusterVector, clusterMergeMap);

        if (clusterMergeMap.empty())
            break;

        this->MergeClusters(clusterVector, clusterMergeMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::MergeClusters(ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterSet clusterVetoList;

    for (const Cluster *const pSeedCluster : clusterVector)
    {
        ClusterList mergeList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterMergeMap, clusterVetoList, mergeList);

        ClusterVector mergeVector(mergeList.begin(), mergeList.end());
        std::sort(mergeVector.begin(), mergeVector.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pAssociatedCluster : mergeVector)
        {
            if (clusterVetoList.count(pAssociatedCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (!pAssociatedCluster->IsAvailable())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            (void) clusterVetoList.insert(pAssociatedCluster);

            if (m_inputClusterListName.empty())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster,
                    m_inputClusterListName, m_inputClusterListName));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(const Cluster *const pSeedCluster, const ClusterMergeMap &clusterMergeMap, ClusterList& associatedClusterList) const
{
    ClusterSet clusterVetoList;
    this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterMergeMap, clusterVetoList, associatedClusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(const Cluster *const pSeedCluster, const Cluster *const pCurrentCluster, const ClusterMergeMap &clusterMergeMap,
    const ClusterSet &clusterVetoList, ClusterList &associatedClusterList) const
{
    if (clusterVetoList.count(pCurrentCluster))
        return;

    ClusterMergeMap::const_iterator iter1 = clusterMergeMap.find(pCurrentCluster);

    if (iter1 == clusterMergeMap.end())
        return;

    ClusterVector associatedClusters(iter1->second.begin(), iter1->second.end());
    std::sort(associatedClusters.begin(), associatedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pAssociatedCluster : associatedClusters)
    {
        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (associatedClusterList.end() != std::find(associatedClusterList.begin(), associatedClusterList.end(), pAssociatedCluster))
            continue;

        associatedClusterList.push_back(pAssociatedCluster);
        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, clusterMergeMap, clusterVetoList, associatedClusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::GetSortedListOfCleanClusters(const ClusterVector &inputClusters, ClusterVector &outputClusters) const
{
    ClusterVector pfoClusters, availableClusters;

    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->IsAvailable())
        {
            pfoClusters.push_back(pCluster);
        }
        else
        {
            availableClusters.push_back(pCluster);
        }
    }

    std::sort(pfoClusters.begin(), pfoClusters.end(), LArClusterHelper::SortByNHits);
    std::sort(availableClusters.begin(), availableClusters.end(), LArClusterHelper::SortByNHits);

    outputClusters.insert(outputClusters.end(), pfoClusters.begin(), pfoClusters.end());
    outputClusters.insert(outputClusters.end(), availableClusters.begin(), availableClusters.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
