/**
 *  @file   larpandoracontent/LArStitching/StitchingAlgorithm.h
 * 
 *  @brief  Header file for the Stitching algorithm class.
 * 
 *  $Log: $
 */
#ifndef PANDORA_STITCHING_ALGORITHM_H
#define PANDORA_STITCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMopUpBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

class StitchingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingAlgorithm class
 */
class StitchingAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    typedef std::unordered_map<const pandora::ParticleFlowObject*, int> PfoToVolumeIdMap;

    /**
     *  @brief  StitchingInfo class
     */
    class StitchingInfo
    {
    public:
        // Placeholder
        PfoToVolumeIdMap    m_pfoToVolumeIdMap;         ///< The pfo to volume id map
    };

    /**
     *  @brief  Get the new/recreated cluster list name
     * 
     *  @return the new/recreated cluster list name
     */
    const std::string &GetNewClusterListName() const;

    /**
     *  @brief  Get the new/recreated vertex list name
     * 
     *  @return the new/recreated vertex list name
     */
    const std::string &GetNewVertexListName() const;

    /**
     *  @brief  Get the new/recreated pfo list name
     * 
     *  @return the new/recreated pfo list name
     */
    const std::string &GetNewPfoListName() const;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<StitchingTool*> StitchingToolVector;
    StitchingToolVector     m_algorithmToolVector;      ///< The algorithm tool vector

    std::string             m_newClusterListName;       ///< The new/recreated cluster list name
    std::string             m_newVertexListName;        ///< The new/recreated vertex list name
    std::string             m_newPfoListName;           ///< The new/recreated pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingTool class
 */
class StitchingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    virtual void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *StitchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new StitchingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewClusterListName() const
{
    return m_newClusterListName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewVertexListName() const
{
    return m_newVertexListName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewPfoListName() const
{
    return m_newPfoListName;
}

} // namespace lar_content

#endif // #ifndef PANDORA_STITCHING_ALGORITHM_H
