/**
 *  @file   larpandoracontent/LArVertex/EnergyKickFeatureTool.h
 * 
 *  @brief  Header file for the energy kick feature tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_ENERGY_KICK_FEATURE_TOOL_H
#define LAR_ENERGY_KICK_FEATURE_TOOL_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnergyKickFeatureTool class
 */
class EnergyKickFeatureTool : public VertexSelectionBaseAlgorithm::VertexFeatureTool<float>
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default constructor
     */
    EnergyKickFeatureTool();

    /**
     *  @brief  Run the tool
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pVertex address of the vertex
     *  @param  slidingFitDataListMap map of the sliding fit data lists
     * 
     *  @return the energy kick feature
     */
    float Run(const VertexSelectionBaseAlgorithm * const pAlgorithm, const pandora::Vertex * const pVertex, 
        const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &slidingFitDataListMap, const VertexSelectionBaseAlgorithm::ClusterListMap &, 
        const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    /**
     *  @brief  Get the energy kick feature for a given view
     * 
     *  @param  vertexPosition2D the projection of the vertex position in this view
     *  @param  slidingFitDataList the list of sliding fit data objects in this view
     * 
     *  @return the energy kick feature
     */
    float GetEnergyKickForView(const pandora::CartesianVector &vertexPosition2D, const VertexSelectionBaseAlgorithm::SlidingFitDataList &slidingFitDataList) const;
                 
    /**
     *  @brief  Increment the energy kick parameters for a given cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  clusterDisplacement the cluster displacement
     *  @param  clusterDirection the cluster direction
     *  @param  totEnergyKick the total energy kick
     *  @param  totEnergy the total energy
     *  @param  totHitKick the total hit kick
     *  @param  totHits the total number of hits
     */
    void IncrementEnergyKickParameters(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterDisplacement, 
        const pandora::CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits) const;

    float    m_rOffset;    ///< The r offset parameter in the energy score
    float    m_xOffset;    ///< The x offset parameter in the energy score
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *EnergyKickFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new EnergyKickFeatureTool();
}

} // namespace lar_content

#endif // #ifndef LAR_ENERGY_KICK_FEATURE_TOOL_H
