/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the energy kick vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H
#define LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnergyKickVertexSelectionAlgorithm class
 */
class EnergyKickVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
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
    EnergyKickVertexSelectionAlgorithm();

private:
     /**
     *  @brief Sliding fit data class.
     */
    class SlidingFitData
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster pointer to the cluster
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, const float slidingFitPitch);
        
        SlidingFitData() : m_minLayerDirection(0.f, 0.f, 0.f), m_maxLayerDirection(0.f, 0.f, 0.f), m_minLayerPosition(0.f, 0.f, 0.f), m_maxLayerPosition(0.f, 0.f, 0.f), m_pCluster(NULL) {};
        
        /**
         *  @brief  Get the min layer direction
         * 
         *  @return the min layer direction
         */
        const pandora::CartesianVector &GetMinLayerDirection() const;

        /**
         *  @brief  Get the max layer direction
         * 
         *  @return the max layer direction
         */
        const pandora::CartesianVector &GetMaxLayerDirection() const;

        /**
         *  @brief  Get the min layer position
         * 
         *  @return the min layer position
         */
        const pandora::CartesianVector &GetMinLayerPosition() const;

        /**
         *  @brief  Get the max layer position
         * 
         *  @return the max layer position
         */
        const pandora::CartesianVector &GetMaxLayerPosition() const;
        
        /**
         *  @brief  Get a pointer to the corresponding cluster
         * 
         *  @return pointer to the corresponding cluster
         */
        const pandora::Cluster *GetCluster() const;

    private:
        pandora::CartesianVector    m_minLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_maxLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_minLayerPosition;     ///< The position of the fit at the max layer
        pandora::CartesianVector    m_maxLayerPosition;     ///< The position of the fit at the max layer
        const pandora::Cluster     *m_pCluster;             ///< Pointer to the corresponding cluster
    };

    typedef std::vector<SlidingFitData> SlidingFitDataList;


     /**
     *  @brief Sliding fit data class.
     */
    class ShowerCluster
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster pointer to the cluster
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        ShowerCluster(const pandora::ClusterList &clusterList);
        
        /**
         *  @brief  Get a pointer to the corresponding cluster
         * 
         *  @return pointer to the corresponding cluster
         */
        pandora::ClusterList GetClusters() const;

    private:
        pandora::ClusterList        m_clusterList;             ///< Pointer to the corresponding cluster
    };

    typedef std::map<const pandora::Cluster *, ShowerCluster> ShowerClusterMap;
    typedef std::vector<ShowerCluster> ShowerClusterList;

    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;
        
    ShowerClusterMap CalculateShowerClusterMap(const SlidingFitDataList &slidingFitDataList) const;

    /**
     *  @brief  Calculate the sliding fits data objects for the clusters in a given view
     * 
     *  @param  slidingFitDataListU to receive the list of sliding fit data objects for u clusters
     *  @param  slidingFitDataListV to receive the list of sliding fit data objects for v clusters
     *  @param  slidingFitDataListW to receive the list of sliding fit data objects for w clusters
     */
    void CalculateClusterSlidingFits(SlidingFitDataList &slidingFitDataListU, SlidingFitDataList &slidingFitDataListV, SlidingFitDataList &slidingFitDataListW) const;

    /**
     *  @brief  Increment the energy kick and energy asymmetry for a given vertex in a given view
     * 
     *  @param  vertexPosition2D the projection of the vertex's position into this view
     *  @param  energyKick the energy kick to increment
     *  @param  energyAsymmetry the energy asymmetry to increment
     *  @param  slidingFitDataList the list of sliding fit data objects in this view
     *  
     *  @return the energy score
     */
    void IncrementEnergyScoresForView(const pandora::CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, float &localEnergyAsymmetry,
        const SlidingFitDataList &slidingFitDataList, const ShowerClusterMap &showerClusterMap, pandora::HitType hitType) const;

    /**
     *  @brief  Increment the parameters used to calculate the energy kick for a given cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  clusterDisplacement distance from the vertex to the closest cluster endpoint (max or min layer of sliding fit)
     *  @param  clusterDirection the cluster direction at closest endpoint (max or min layer of sliding fit)
     *  @param  totEnergyKick the total energy kick to increment
     *  @param  totEnergy the total energy to increment
     *  @param  totHitKick the total hit kick to increment
     *  @param  totHits the total number of hits to increment
     */
    void IncrementEnergyKickParameters(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterDisplacement,
        const pandora::CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits, const SlidingFitDataList &closeSlidingFitData, const ShowerClusterMap &showerClusterMap, pandora::HitType hitType) const;

    float IsClusterShowerLike(const pandora::Cluster *const pCluster, const SlidingFitDataList &closeSlidingFitData, const ShowerClusterMap &showerClusterMap, pandora::HitType hitType) const;
    bool IsClusterShowerLike(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Increment the parameters used to calculate the energy asymmetry for a given cluster and a given vertex
     * 
     *  @param  weight the weight to apply to the local direction measure (cluster energy or number of hits)
     *  @param  clusterDirection the cluster direction at closest endpoint (max or min layer of sliding fit)
     *  @param  localWeightedDirectionSum current local event axis sum using energy or hit weighting
     * 
     */
    void IncrementEnergyAsymmetryParameters(const float weight, const pandora::CartesianVector &clusterDirection,
         pandora::CartesianVector &localWeightedDirectionSum) const;

    /**
     *  @brief  Calculate the energy asymmetry for a vertex in a given view using the calculated parameters
     * 
     *  @param  useEnergyMetrics whether to use the energy metrics, or to revert to hit-based metrics
     *  @param  vertexPosition2D the projection of the vertex's position into this view
     *  @param  slidingFitDataList the list of sliding fits for all the 2D clusters
     *  @param  localWeightedDirection current local event axis using energy or hit weighting
     */
    void IncrementEnergyAsymmetry(float &energyAsymmetry, float &localEnergyAsymmetry, pandora::ClusterList &localAsymmetryClusterList, const bool useEnergyMetrics, const pandora::CartesianVector &vertexPosition2D,
        const SlidingFitDataList &slidingFitDataList, const pandora::CartesianVector &localWeightedDirection) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;        ///< The list of cluster list names
    unsigned int            m_minClusterCaloHits;           ///< The min number of hits parameter in the energy score
    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits

    float                   m_rOffset;                      ///< The r offset parameter in the energy score
    float                   m_xOffset;                      ///< The x offset parameter in the energy score
    float                   m_epsilon;                      ///< The epsilon parameter in the energy score

    float                   m_asymmetryConstant;            ///< The asymmetry constant parameter in the energy score
    float                   m_maxAsymmetryDistance;         ///< The max distance between cluster (any hit) and vertex to calculate asymmetry score
    
    float                   m_localAsymmetryConstant;       ///< ...
    unsigned int            m_minLocalAsymmetryClusterHits; ///< ...
    bool                    m_useLocalAsymmetry;            ///< ...

    float                   m_beamDeweightingConstant;      ///< ...
    float                   m_showerDeweightingConstant;    ///< ...
    float                   m_showerCollapsingConstant;     ///< ...
    float                   m_minShowerSpineLength;         ///< ...
    
    bool                    m_useClusterCharacterisation;   ///< ...
    float                   m_maxTrackLikeDeviation;        ///< ...
    
    float m_showerClusteringDistance;
    
    float m_showerAngleConstant;
    float m_showerDistanceConstant;
    float m_vertexClusterDistance;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EnergyKickVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new EnergyKickVertexSelectionAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMinLayerDirection() const
{
    return m_minLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMaxLayerDirection() const
{
    return m_maxLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMinLayerPosition() const
{
    return m_minLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMaxLayerPosition() const
{
    return m_maxLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::ClusterList EnergyKickVertexSelectionAlgorithm::ShowerCluster::GetClusters() const
{
    return m_clusterList;
}

} // namespace lar_content

#endif // #ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H
