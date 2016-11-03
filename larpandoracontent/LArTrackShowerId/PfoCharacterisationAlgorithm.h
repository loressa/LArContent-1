/**
 *  @file   larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoCharacterisationAlgorithm class
 */
class PfoCharacterisationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    PfoCharacterisationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Whether pfo is identified as a clear track
     *
     *  @param  pPfo address of the relevant pfo
     * 
     *  @return boolean
     */
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_trackPfoListName;         ///< The track pfo list name
    std::string             m_showerPfoListName;        ///< The shower pfo list name
    pandora::StringVector   m_inputPfoListNames;        ///< The names of the input pfo lists

    bool                    m_updateClusterIds;         ///< Whether to update daughter cluster particle id labels to match pfo id
    bool                    m_postBranchAddition;       ///< Whether to use configuration for shower clusters post branch addition
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PfoCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new PfoCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
