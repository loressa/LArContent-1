/**
 *  @file   larpandoracontent/LArVertex/SVMVertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the SVM vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SVM_VERTEX_SELECTION_ALGORITHM_H
#define LAR_SVM_VERTEX_SELECTION_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{    
/**
 *  @brief  SVMVertexSelectionAlgorithm class
 */
class SVMVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
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
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief Vertex feature info class
     */
    class VertexFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  beamDeweighting the beam deweighting feature
         *  @param  rPhiFeature the r/phi feature
         *  @param  energyKick the energy kick feature
         *  @param  localAsymmetry the local asymmetry feature
         *  @param  globalAsymmetry the global asymmetry feature
         *  @param  showerAsymmetry the shower asymmetry feature
         */
        VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, const float energyKick, const float localAsymmetry,
                          const float globalAsymmetry, const float showerAsymmetry);
                          
        float    m_beamDeweighting;    ///< The beam deweighting feature
        float    m_rPhiFeature;        ///< The r/phi feature
        float    m_energyKick;         ///< The energy kick feature
        float    m_localAsymmetry;     ///< The local asymmetry feature
        float    m_globalAsymmetry;    ///< The global asymmetry feature
        float    m_showerAsymmetry;    ///< The shower asymmetry feature
    };
    
    typedef std::map<const pandora::Vertex * const, VertexFeatureInfo> VertexFeatureInfoMap;
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief Event feature info class
     */
    class EventFeatureInfo
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  eventShoweryness the event showeryness feature
         *  @param  eventEnergy the energy of the event
         *  @param  eventVolume the volume of the event
         *  @param  longitudinality the longitudinality of the event
         *  @param  nHits the number of hits in the event
         *  @param  nClusters the number of clusters in the event
         *  @param  nCandidates the total number of vertex candidates
         */
        EventFeatureInfo(const float eventShoweryness, const float eventEnergy, const float eventVolume, const float longitudinality,
            const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates);
                          
        float           m_eventShoweryness;        ///< The event showeryness feature
        float           m_eventEnergy;             ///< The event energy
        float           m_eventVolume;             ///< The volume of the event
        float           m_longitudinality;         ///< The longitudinality of the event
        unsigned int    m_nHits;                   ///< The number of hits in the event
        unsigned int    m_nClusters;               ///< The number of clusters in the event
        unsigned int    m_nCandidates;             ///< The total number of vertex candidates
    };
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default constructor
     */
    SVMVertexSelectionAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::pair<pandora::CartesianVector, pandora::CartesianVector> ClusterEndPoints;
    typedef std::map<const pandora::Cluster *const, ClusterEndPoints> ClusterEndPointsMap;
    typedef std::vector<pandora::FloatVector> FeatureListVector;
    typedef std::vector<pandora::VertexVector> VectorOfVertexVectors; 
    
    /**
     *  @brief  Get the vertex score list
     * 
     *  @param  vertexVector the vector of vertices
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the hit kd tree for the U view
     *  @param  kdTreeV the hit kd tree for the V view
     *  @param  kdTreeW the hit kd tree for the W view
     *  @param  vertexScoreList the vertex score list to fill
     */
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, 
        HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;
    
    /**
     *  @brief  Calculate the shower cluster map for a cluster list
     * 
     *  @param  inputClusterList the input cluster list
     *  @param  showerClusterList the shower cluster list to populate
     */
    void CalculateShowerClusterList(const pandora::ClusterList &inputClusterList, ShowerClusterList &showerClusterList) const;
    
    /**
     *  @brief  Add the endpoints of any shower-like clusters to the map
     * 
     *  @param  clusterList the list of clusters
     *  @param  showerLikeClusters the list of shower-like clusters to populate
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints to populate
     */
    void GetShowerLikeClusterEndPoints(const pandora::ClusterList &clusterList, pandora::ClusterList &showerLikeClusters, 
        ClusterEndPointsMap &clusterEndPointsMap) const;
                                        
    /**
     *  @brief  Try to add an available cluster to a given shower cluster, considering distances from a given member of that shower cluster
     * 
     *  @param  clusterEndPointsMap the map of shower-like cluster endpoints 
     *  @param  availableShowerLikeClusters the list of shower-like clusters still available
     *  @param  pCluster the cluster in the shower cluster from which to consider distances
     *  @param  showerCluster the shower cluster
     */
    bool AddClusterToShower(const ClusterEndPointsMap &clusterEndPointsMap, pandora::ClusterList &availableShowerLikeClusters, 
        const pandora::Cluster *const pCluster, pandora::ClusterList &showerCluster) const;
    
    /**
     *  @brief  Calculate the event parameters
     * 
     *  @param  clusterListU the U-view cluster list
     *  @param  clusterListV the V-view cluster list
     *  @param  clusterListW the W-view cluster list
     *  @param  vertexVector the vector of vertex candidates
     * 
     *  @return the event feature info object
     */
    EventFeatureInfo CalculateEventFeatures(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV, 
        const pandora::ClusterList &clusterListW, const pandora::VertexVector &vertexVector) const;
                        
    /**
     *  @brief  Find whether a cluster is shower-like
     * 
     *  @param  pCluster the cluster
     * 
     *  @return whether the cluster is shower-like
     */
    bool IsClusterShowerLike(const pandora::Cluster *const pCluster) const;
               
    /**
     *  @brief  Increment the showery hit parameters for a cluster list
     * 
     *  @param  clusterList the cluster list
     *  @param  nShoweryHits the number of showery hits
     *  @param  nHits the number of hits
     *  @param  eventEnergy the event energy
     */
    void IncrementShoweryParameters(const pandora::ClusterList &clusterList, unsigned int &nShoweryHits, unsigned int &nHits, 
       float &eventEnergy) const;
    
    /**
     *  @brief  Get the span of a set of vertex candidates in a certain dimension
     * 
     *  @param  vertexVector the vector of vertex candidates
     *  @param  getCoord the function to extract the required coordinate from the vertex
     * 
     *  @return the span
     */
    float GetCandidateSpan(const pandora::VertexVector &vertexVector, const std::function<float(const pandora::Vertex *const)> &getCoord) const;
    
    /**
     *  @brief  Populate the vertex feature info map for a given vertex
     * 
     *  @param  beamConstants the beam constants
     *  @param  slidingFitDataListMap the sliding fit data list map
     *  @param  clusterListMap the cluster list map
     *  @param  kdTreeMap the kd tree map
     *  @param  pVertex the vertex
     *  @param  vertexFeatureInfoMap the map to populate
     */
    void PopulateVertexFeatureInfoMap(const BeamConstants &beamConstants, const ClusterListMap &clusterListMap, const SlidingFitDataListMap &slidingFitDataListMap,
        const ShowerClusterListMap &showerClusterListMap, const KDTreeMap &kdTreeMap, const pandora::Vertex *const pVertex, 
        VertexFeatureInfoMap &vertexFeatureInfoMap) const;
    
    /**
     *  @brief  Populate the initial vertex score list for a given vertex
     * 
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pVertex the vertex
     *  @param  initialScoreList the score list to populate
     */
    void PopulateInitialScoreList(VertexFeatureInfoMap &vertexFeatureInfoMap, const pandora::Vertex *const pVertex, VertexScoreList &initialScoreList) const;
    
    /**
     *  @brief  Get the list of top-N separated vertices
     * 
     *  @param  initialScoreList the initial score list
     *  @param  topNVertices the top-N vertex list to populate
     */
    void GetTopNVertices(VertexScoreList &initialScoreList, pandora::VertexList &topNVertices) const;
    
    /**
     *  @brief  Use the MC information to get the best vertex from the top-N list
     * 
     *  @param  topNVertices ahe top-N list
     *  @param  pBestVertex address of the best vertex
     *  @param  bestVertexDr dR of the best vertex
     */
    void GetBestVertex(const pandora::VertexList &topNVertices, const pandora::Vertex *&pBestVertex, float &bestVertexDr) const;
    
    void GetCheatedVertex(const pandora::VertexVector &vertexVector, const pandora::Vertex *&pCheatedVertex) const; // ATTN temporary

    /**
     *  @brief  Generate the feature list for a given vertex
     * 
     *  @param  eventFeatureInfo the event feature info
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pVertex the vertex
     *  @param  topNVertices the top-N vertex list
     * 
     *  @return the feature list
     */
    pandora::FloatVector GenerateFeatureList(const EventFeatureInfo &eventFeatureInfo, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        const pandora::Vertex *const pVertex, const pandora::VertexList &topNVertices) const;
    
    /**
     *  @brief  Add the event features to a vector in the correct order
     * 
     *  @param  eventFeatureInfo the event feature info
     *  @param  featureVector the vector of floats to append
     */
    void AddEventFeaturesToVector(const EventFeatureInfo &eventFeatureInfo, pandora::FloatVector &featureVector) const;
    
    /**
     *  @brief  Add the vertex features to a vector in the correct order
     * 
     *  @param  vertexFeatureInfo the vertex feature info
     *  @param  featureVector the vector of floats to append
     */
    void AddVertexFeaturesToVector(const VertexFeatureInfo &vertexFeatureInfo, pandora::FloatVector &featureVector) const;
    
    /**
     *  @brief  Generate all the permuted feature lists for a given vertex
     * 
     *  @param  eventFeatureInfo the event feature info
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  pVertex the vertex
     *  @param  topNVertices the top-N vertex list
     * 
     *  @return the vector of feature lists
     */
    FeatureListVector GeneratePermutedFeatureLists(const EventFeatureInfo &eventFeatureInfo, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        const pandora::Vertex *const pVertex, const pandora::VertexList &topNVertices) const;
    
    /**
     *  @brief  Populate the final vertex score list using the r/phi score to find the best vertex in the vicinity
     * 
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  bestVertexScoreList the list of best vertex scores
     *  @param  vertexVector the vector of all vertex candidates
     *  @param  finalVertexScoreList the final vertex score list to populate
     */
    void PopulateFinalVertexScoreList(const VertexFeatureInfoMap &vertexFeatureInfoMap, VertexScoreList &bestVertexScoreList,
        const pandora::VertexVector &vertexVector, VertexScoreList &finalVertexScoreList) const;
    
    VertexFeatureToolBase::FeatureToolMap m_featureToolMap;   ///< The feature tool map
    SupportVectorMachine<37>              m_svMachine;        ///< The support vector machine
                                                              
    std::string           m_trainingOutputFile;               ///< The training output file
    std::string           m_parameterInputFile;               ///< The parameter input file
    std::string           m_svmName;                          ///< The name of the SVM to find
    std::string           m_mcParticleListName;               ///< The MC particle list for creating training examples
                                                              
    pandora::StringVector m_inputClusterListNames;            ///< The list of cluster list names
    
    bool                  m_classifyUsingPermutations;        ///< Whether to classify using permutations instead of top-5
    bool                  m_trainingSetMode;                  ///< Whether to train
    bool                  m_allowClassifyDuringTraining;      ///< Whether classification is allowed during training
    bool                  m_produceAllTrainingPermutations;   ///< Whether to produce training sets with all top-N! permutations
    float                 m_mcVertexXCorrection;              ///< The correction to the x-coordinate of the MC vertex position
                                                              
    unsigned int          m_minClusterCaloHits;               ///< The min number of hits parameter in the energy score
    unsigned int          m_slidingFitWindow;                 ///< The layer window for the sliding linear fits
    float                 m_minShowerSpineLength;             ///< The minimum length at which all are considered to be tracks
                                                              
    float                 m_beamDeweightingConstant;          ///< The beam deweighting constant
    float                 m_localAsymmetryConstant;           ///< The local asymmetry constant
    float                 m_globalAsymmetryConstant;          ///< The global asymmetry constant
    float                 m_showerAsymmetryConstant;          ///< The shower asymmetry constant
    float                 m_energyKickConstant;               ///< The energy kick constant
                                                              
    float                 m_minTopNSeparation;                ///< The minimum separation of the best top-N vertices
    unsigned int          m_topNSize;                         ///< The number of best vertices to consider
    
    float                 m_showerClusteringDistance;         ///< The shower clustering distance
    unsigned int          m_minShowerClusterHits;             ///< The minimum number of shower cluster hits
    bool                  m_useShowerClusteringApproximation; ///< Whether to use the shower clustering distance approximation
    
    bool m_drawThings;         ///< ATTN temporary
    bool m_cheatTrackShowerId; ///< ATTN temporary
    bool m_cheatTheVertex;     ///< ATTN temporary
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SVMVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new SVMVertexSelectionAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SVMVertexSelectionAlgorithm::VertexFeatureInfo::VertexFeatureInfo(const float beamDeweighting, const float rPhiFeature, 
    const float energyKick, const float localAsymmetry, const float globalAsymmetry, const float showerAsymmetry) :
    m_beamDeweighting(beamDeweighting),
    m_rPhiFeature(rPhiFeature),
    m_energyKick(energyKick),
    m_localAsymmetry(localAsymmetry),
    m_globalAsymmetry(globalAsymmetry),
    m_showerAsymmetry(showerAsymmetry)
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
inline SVMVertexSelectionAlgorithm::EventFeatureInfo::EventFeatureInfo(const float eventShoweryness, const float eventEnergy, 
    const float eventVolume, const float longitudinality, const unsigned int nHits, const unsigned int nClusters, const unsigned int nCandidates) :
    m_eventShoweryness(eventShoweryness),
    m_eventEnergy(eventEnergy),
    m_eventVolume(eventVolume),
    m_longitudinality(longitudinality),
    m_nHits(nHits),
    m_nClusters(nClusters),
    m_nCandidates(nCandidates)
{
}

} // namespace lar_content

#endif // #ifndef LAR_SVM_VERTEX_SELECTION_ALGORITHM_H
