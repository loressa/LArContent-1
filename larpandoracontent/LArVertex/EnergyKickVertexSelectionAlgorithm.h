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
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

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
     *  @brief Kernel estimate class
     */
    class KernelEstimate
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  sigma the width associated with the kernel estimate
         */
        KernelEstimate(const float sigma);

        /**
         *  @brief  Sample the parameterised distribution at a specified x coordinate
         * 
         *  @param  x the position at which to sample
         * 
         *  @return the sample value
         */
        float Sample(const float x) const;

        typedef std::multimap<float, float> ContributionList;   ///< Map from x coord to weight, ATTN avoid map.find, etc. with float key

        /**
         *  @brief  Get the contribution list
         * 
         *  @return the contribution list
         */
        const ContributionList &GetContributionList() const;

        /**
         *  @brief  Get the assigned width
         * 
         *  @return the assigned width
         */
        float GetSigma() const;

        /**
         *  @brief  Add a contribution to the distribution
         * 
         *  @param  x the position
         *  @param  weight the weight
         */
        void AddContribution(const float x, const float weight);

    private:
        ContributionList            m_contributionList;         ///< The contribution list
        const float                 m_sigma;                    ///< The assigned width
    };

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
        SlidingFitData(const pandora::ClusterVector &clusterVector, const int slidingFitWindow, const float slidingFitPitch, int minClusterCaloHits);
        
        SlidingFitData() : m_minLayerDirection(0.f, 0.f, 0.f), m_maxLayerDirection(0.f, 0.f, 0.f), m_minLayerPosition(0.f, 0.f, 0.f), m_maxLayerPosition(0.f, 0.f, 0.f) {};
        
        /**
         *  @brief  Get the min layer direction
         * 
         *  @return the min layer directionslidingFitData
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
        const pandora::ClusterVector &GetClusterVector() const;

    private:
        pandora::CartesianVector    m_minLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_maxLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_minLayerPosition;     ///< The position of the fit at the max layer
        pandora::CartesianVector    m_maxLayerPosition;     ///< The position of the fit at the max layer
        pandora::ClusterVector  m_clusterVector;             ///< Pointer to the corresponding cluster
    };

    typedef std::vector<SlidingFitData> SlidingFitDataList;
    
    
     /**
     *  @brief Sliding fit data class.
     */
    class ShowerCluster
    {/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the energy kick vertex selection algorithm class.
 * 
 *  $Log: $
 */

    public:

        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster pointer to the cluster
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        ShowerCluster(const pandora::ClusterVector &clusterVector, const int slidingFitWindow, 
        const float slidingFitPitch);
        
        /**
         *  @brief  Get a pointer to the corresponding cluster
         * 
         *  @return pointer to the corresponding cluster
         */
        const pandora::ClusterVector &GetClusters() const;
        
        const TwoDSlidingFitResult &GetFit() const { return m_twoDSlidingFitResult; }
        
    private:
        pandora::ClusterVector        m_clusterVector;             ///< Pointer to the corresponding cluster
        TwoDSlidingFitResult m_twoDSlidingFitResult;        
    };

    typedef std::map<const pandora::Cluster *, ShowerCluster> ShowerClusterMap;
    typedef std::vector<ShowerCluster> ShowerClusterList;


   /**
 *  @brief   InteractionType enum
 */
enum InteractionType
{
    // TODO Move to dynamic interaction type identification and labelling
    CCQEL_MU,
    CCQEL_MU_P,
    CCQEL_MU_P_P,
    CCQEL_MU_P_P_P,
    CCQEL_MU_P_P_P_P,
    CCQEL_MU_P_P_P_P_P,
    CCQEL_E,
    CCQEL_E_P,
    CCQEL_E_P_P,
    CCQEL_E_P_P_P,
    CCQEL_E_P_P_P_P,
    CCQEL_E_P_P_P_P_P,
    NCQEL_P,
    NCQEL_P_P,
    NCQEL_P_P_P,
    NCQEL_P_P_P_P,
    NCQEL_P_P_P_P_P,
    CCRES_MU,
    CCRES_MU_P,
    CCRES_MU_P_P,
    CCRES_MU_P_P_P,
    CCRES_MU_P_P_P_P,
    CCRES_MU_P_P_P_P_P,
    CCRES_MU_PIPLUS,
    CCRES_MU_P_PIPLUS,
    CCRES_MU_P_P_PIPLUS,
    CCRES_MU_P_P_P_PIPLUS,
    CCRES_MU_P_P_P_P_PIPLUS,
    CCRES_MU_P_P_P_P_P_PIPLUS,
    CCRES_MU_PHOTON,
    CCRES_MU_P_PHOTON,
    CCRES_MU_P_P_PHOTON,
    CCRES_MU_P_P_P_PHOTON,
    CCRES_MU_P_P_P_P_PHOTON,
    CCRES_MU_P_P_P_P_P_PHOTON,
    CCRES_MU_PIZERO,
    CCRES_MU_P_PIZERO,
    CCRES_MU_P_P_PIZERO,
    CCRES_MU_P_P_P_PIZERO,
    CCRES_MU_P_P_P_P_PIZERO,
    CCRES_MU_P_P_P_P_P_PIZERO,
    CCRES_E,
    CCRES_E_P,
    CCRES_E_P_P,
    CCRES_E_P_P_P,
    CCRES_E_P_P_P_P,
    CCRES_E_P_P_P_P_P,
    CCRES_E_PIPLUS,
    CCRES_E_P_PIPLUS,
    CCRES_E_P_P_PIPLUS,
    CCRES_E_P_P_P_PIPLUS,
    CCRES_E_P_P_P_P_PIPLUS,
    CCRES_E_P_P_P_P_P_PIPLUS,
    CCRES_E_PHOTON,
    CCRES_E_P_PHOTON,
    CCRES_E_P_P_PHOTON,
    CCRES_E_P_P_P_PHOTON,
    CCRES_E_P_P_P_P_PHOTON,
    CCRES_E_P_P_P_P_P_PHOTON,
    CCRES_E_PIZERO,
    CCRES_E_P_PIZERO,
    CCRES_E_P_P_PIZERO,
    CCRES_E_P_P_P_PIZERO,
    CCRES_E_P_P_P_P_PIZERO,
    CCRES_E_P_P_P_P_P_PIZERO,
    NCRES_P,
    NCRES_P_P,
    NCRES_P_P_P,
    NCRES_P_P_P_P,
    NCRES_P_P_P_P_P,
    NCRES_PIPLUS,
    NCRES_P_PIPLUS,
    NCRES_P_P_PIPLUS,
    NCRES_P_P_P_PIPLUS,
    NCRES_P_P_P_P_PIPLUS,
    NCRES_P_P_P_P_P_PIPLUS,
    NCRES_PIMINUS,
    NCRES_P_PIMINUS,
    NCRES_P_P_PIMINUS,
    NCRES_P_P_P_PIMINUS,
    NCRES_P_P_P_P_PIMINUS,
    NCRES_P_P_P_P_P_PIMINUS,
    NCRES_PHOTON,
    NCRES_P_PHOTON,
    NCRES_P_P_PHOTON,
    NCRES_P_P_P_PHOTON,
    NCRES_P_P_P_P_PHOTON,
    NCRES_P_P_P_P_P_PHOTON,
    NCRES_PIZERO,
    NCRES_P_PIZERO,
    NCRES_P_P_PIZERO,
    NCRES_P_P_P_PIZERO,
    NCRES_P_P_P_P_PIZERO,
    NCRES_P_P_P_P_P_PIZERO,
    CCDIS,
    NCDIS,
    CCCOH,
    NCCOH,
    OTHER_INTERACTION,
    NU_E_SCATTERING,
    ALL_INTERACTIONS // ATTN use carefully!
};

    void SelectCaloHits(const pandora::CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    pandora::CaloHitList &selectedCaloHitList) const;
    
    bool PassMCParticleChecks(const pandora::MCParticle *const pOriginalPrimary, const pandora::MCParticle *const pThisMCParticle,
    const pandora::MCParticle *const pHitMCParticle) const;
    
    void SelectGoodCaloHits(const pandora::CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    pandora::CaloHitList &selectedGoodCaloHitList) const;

    InteractionType GetInteractionType(const LArMCParticle *const pLArMCNeutrino, const pandora::MCParticleList *pMCParticleList, const LArMonitoringHelper::MCContributionMap &mcToGoodTrueHitListMap) const;
    std::string ToString(const InteractionType interactionType) const;

    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    void SelectTrueNeutrinos(const pandora::MCParticleList *const pAllMCParticleList, pandora::MCParticleVector &selectedMCNeutrinoVector) const;

    float GetEventHitShoweryness(const pandora::ClusterList &inputClusterListU, const pandora::ClusterList &inputClusterListV, const pandora::ClusterList &inputClusterListW) const;
    float GetEventClusterShoweryness(const pandora::ClusterList &inputClusterListU, const pandora::ClusterList &inputClusterListV, const pandora::ClusterList &inputClusterListW) const;

    void CalculateShowerClusterMap(const pandora::ClusterList &clusterList, ShowerClusterMap &showerClusterMap) const;

    /**
     *  @brief  Calculate the sliding fits data objects for the clusters in a given view
     * 
     *  @param  slidingFitDataListU to receive the list of sliding fit data objects for u clusters
     *  @param  slidingFitDataListV to receive the list of sliding fit data objects for v clusters
     *  @param  slidingFitDataListW to receive the list of sliding fit data objects for w clusters
     */
    void CalculateClusterSlidingFits(const pandora::ClusterList &inputClusterList, SlidingFitDataList &singleClusterSlidingFitDataList) const;

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
    void IncrementEnergyScoresForView(const pandora::CartesianVector &vertexPosition2D, float &energyKick, float &energyAsymmetry, float &globalEnergyAsymmetry,
        const ShowerClusterMap &showerClusterMap, float &showerEnergyAsymmetry, const SlidingFitDataList &singleClusterSlidingFitDataList) const;

    void IncrementGlobalEnergyAsymmetry(float &globalEnergyAsymmetry, const bool useEnergyMetrics, const pandora::CartesianVector &vertexPosition2D,
                    const SlidingFitDataList &slidingFitDataList, const pandora::CartesianVector &localWeightedDirectionSum) const;

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
        const pandora::CartesianVector &clusterDirection, float &totEnergyKick, float &totEnergy, float &totHitKick, unsigned int &totHits) const;

    bool IsClusterShowerLike(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Increment the parameters used to calculate the energy asymmetry for a given cluster and a given vertex
     * 
     *  @param  weight the weight to apply to the local direction measure (cluster energy or number of hits)
     *  @param  clusterDirection the cluster direction at closest endpoint (max or min layer of sliding fit)
     *  @param  localWeightedDirectionSum current local event axis sum using energy or hit weighting
     * 
     *  @return whether the asymmetry calculation is viable
     */
    bool IncrementEnergyAsymmetryParameters(const float weight, const pandora::CartesianVector &clusterDirection,
        pandora::CartesianVector &localWeightedDirectionSum) const;

    /**
     *  @brief  Calculate the energy asymmetry for a vertex in a given view using the calculated parameters
     * 
     *  @param  useEnergyMetrics whether to use the energy metrics, or to revert to hit-based metrics
     *  @param  vertexPosition2D the projection of the vertex's position into this view
     *  @param  asymmetryClusters the list of clusters considered for the energy asymmetry calculation
     *  @param  localWeightedDirection current local event axis using energy or hit weighting
     */
    float CalculateEnergyAsymmetry(const bool useEnergyMetrics, const pandora::CartesianVector &vertexPosition2D,
        const pandora::ClusterVector &asymmetryClusters, const pandora::CartesianVector &localWeightedDirection) const;
        
    /**
     *  @brief  Get the score for a trio of kernel estimations, using fast histogram approach
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the fast score
     */
    float GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Get the score for a trio of kernel estimations, using kernel density estimation but with reduced (binned) sampling
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the midway score
     */
    float GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Get the score for a trio of kernel estimations, using kernel density estimation and full hit-by-hit sampling
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the full score
     */
    float GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Use hits in clusters (in the provided kd tree) to fill a provided kernel estimate with hit-vertex relationship information
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *  @param  kernelEstimate to receive the populated kernel estimate
     */
    void FillKernelEstimate(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its spatial position in relation to other selected candidates
     * 
     *  @param  pVertex the address of the vertex
     *  @param  selectedVertexList the selected vertex list
     * 
     *  @return boolean
     */
    bool AcceptVertexLocation(const pandora::Vertex *const pVertex, const pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Fast estimate of std::atan2 function. Rather coarse (max |error| > 0.01) but should suffice for this use-case.
     * 
     *  @param  y the y coordinate
     *  @param  x the x coordinate
     * 
     *  @return estimate of std::atan2
     */
    float atan2Fast(const float y, const float x) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;        ///< The list of cluster list names
    unsigned int            m_minClusterCaloHits;           ///< The min number of hits parameter in the energy score
    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits

    float                   m_rOffset;                      ///< The r offset parameter in the energy score
    float                   m_xOffset;                      ///< The x offset parameter in the float &globalEnergyAsymmetryenergy score
    float                   m_epsilon;                      ///< The epsilon parameter in the energy score

    float                   m_maxAsymmetryDistance;         ///< The max distance between cluster (any hit) and vertex to calculate asymmetry score
    float                   m_minAsymmetryCosAngle;         ///< The min opening angle cosine used to determine viability of asymmetry score
    unsigned int            m_maxAsymmetryNClusters;        ///< The max number of associated clusters to calculate the asymmetry
    
    
    float m_showerDeweightingConstant;
    float m_showerCollapsingConstant;
    float m_minShowerSpineLength;
    float m_showerClusteringDistance;
    float m_vertexClusterDistance;
    
    typedef std::pair<const ShowerCluster * const, const SlidingFitData * const> ShowerVertexPair;
    typedef std::map<ShowerVertexPair, float> ShowerDataMap;
    typedef std::map<const pandora::Cluster *const, bool> ShowerLikeClusterMap;
    
    mutable ShowerDataMap m_showerDataMap;
    mutable ShowerLikeClusterMap m_showerLikeClusterMap;
    
    unsigned int m_minShowerClusterHits;
    bool m_useShowerClusteringApproximation;
    bool m_cheatTrackShowerId;
    
    
    
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_caloHitListName;           ///< Name of input MC particle list
    
    
    bool            m_fastScoreCheck;               ///< Whether to use the fast histogram based score to selectively avoid calling full or midway scores
    bool            m_fastScoreOnly;                ///< Whether to use the fast histogram based score only
    bool            m_fullScore;                    ///< Whether to use the full kernel density estimation score, as opposed to the midway score
 
    float           m_kernelEstimateSigma;          ///< The Gaussian width to use for kernel estimation
    float           m_kappa;                        ///< Hit-deweighting offset, of form: weight = 1 / sqrt(distance + kappa), units cm
    float           m_maxHitVertexDisplacement1D;   ///< Max hit-vertex displacement in *any one dimension* for contribution to kernel estimation

    float           m_minFastScoreFraction;         ///< Fast score must be at least this fraction of best fast score to calculate full score
    unsigned int    m_fastHistogramNPhiBins;        ///< Number of bins to use for fast score histograms
    float           m_fastHistogramPhiMin;          ///< Min value for fast score histograms
    float           m_fastHistogramPhiMax;          ///< Max value for fast score histograms

    bool            m_enableFolding;                ///< Whether to enable folding of -pi -> +pi phi distribution into 0 -> +pi region only
    
    std::string m_trainingSetPrefix;
    
 
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

inline const pandora::ClusterVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetClusterVector() const
{
    return m_clusterVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterVector &EnergyKickVertexSelectionAlgorithm::ShowerCluster::GetClusters() const
{
    return m_clusterVector;
}


inline EnergyKickVertexSelectionAlgorithm::KernelEstimate::KernelEstimate(const float sigma) :
    m_sigma(sigma)
{
    if (m_sigma < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const EnergyKickVertexSelectionAlgorithm::KernelEstimate::ContributionList &EnergyKickVertexSelectionAlgorithm::KernelEstimate::GetContributionList() const
{
    return m_contributionList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float EnergyKickVertexSelectionAlgorithm::KernelEstimate::GetSigma() const
{
    return m_sigma;
}


} // namespace lar_content

#endif // #ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H