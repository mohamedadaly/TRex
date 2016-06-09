/*
-----------------------------------------------------------------------
Copyright: 2010-2015, iMinds-Vision Lab, University of Antwerp
           2014-2015, CWI, Amsterdam

Contact: astra@uantwerpen.be
Website: http://sf.net/projects/astra-toolbox

This file is part of the ASTRA Toolbox.


The ASTRA Toolbox is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The ASTRA Toolbox is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the ASTRA Toolbox. If not, see <http://www.gnu.org/licenses/>.

-----------------------------------------------------------------------
$Id$
*/

#ifndef _INC_ASTRA_TREXALGORITHM
#define _INC_ASTRA_TREXALGORITHM

#include "Globals.h"
#include "Config.h"

#include "Algorithm.h"
//#include "ReconstructionAlgorithm2D.h"
#include "SartAlgorithm.h"

#include "Projector2D.h"
#include "Float32ProjectionData2D.h"
#include "Float32VolumeData2D.h"
#include "Float32VolumeData3DMemory.h"
#include "AstraObjectManager.h"

#include "DataProjector.h"

namespace astra {

// ----------------------------------------------------------------------------
// Base class for priors
class _AstraExport CTRexPrior {
public:
	float32 m_fSigma;
	CTRexPrior() : m_fSigma(1.f) {}	
	CTRexPrior(float32 _fSigma) : m_fSigma(_fSigma) {}	
	virtual ~CTRexPrior() {}

	// Multiply the volume by the matrix K
	virtual void K(const CFloat32VolumeData2D* pX,  
		CFloat32VolumeData3DMemory* pKx) = 0;

	// Multiply the volume by the transpose of K
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX) = 0;

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV) = 0;

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm() = 0;

	// Return the depth of the object returned by applying K
	virtual int depth() = 0;
};

// Anisotropic TV
class _AstraExport CTRexPriorATV : public CTRexPrior {
public:
	CTRexPriorATV() : CTRexPrior() {}
	CTRexPriorATV(float _fSigma) : CTRexPrior(_fSigma) {}

	// Multiply the volume by the matrix K * sigma
	virtual void K(const CFloat32VolumeData2D* pX, 
		CFloat32VolumeData3DMemory* pKx);

	// Multiply the volume by the transpose of K * sigma
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX);

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV);

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm();

	virtual int depth() { return 2; }
};

// Isotropic TV
class _AstraExport CTRexPriorITV : public CTRexPriorATV {
public:
	CTRexPriorITV() : CTRexPriorATV() {}
	CTRexPriorITV(float _fSigma) : CTRexPriorATV(_fSigma) {}

	// Apply the proximal operator on the input u (of size Kx).
	virtual void prox(const CFloat32VolumeData3DMemory* pU, float32 rho,
		CFloat32VolumeData3DMemory* pV);

};

// Sum of Absolute Differences
class _AstraExport CTRexPriorSAD : public CTRexPriorATV {
public:
	CTRexPriorSAD() : CTRexPriorATV() {}
	CTRexPriorSAD(float _fSigma) : CTRexPriorATV(_fSigma) {}

	// Multiply the volume by the matrix K * sigma
	virtual void K(const CFloat32VolumeData2D* pX, 
		CFloat32VolumeData3DMemory* pKx);

	// Multiply the volume by the transpose of K * sigma
	virtual void Kt(const CFloat32VolumeData3DMemory* pKx, 
		CFloat32VolumeData2D* pX);

	// Return the squared norm of K ||K||^2_2
	virtual float32 norm();

	virtual int depth() { return 8; }

};
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Base class for data terms
class _AstraExport CTRexData{
public:
	CReconstructionAlgorithm2D* m_pAlg;
	Config* m_pCfg;
	// ProxInput for the prox operator
	CFloat32VolumeData2D* m_pProxInput;



	CTRexData() {
		m_pAlg = NULL;
		m_pCfg = NULL;
		m_pProxInput = NULL;
	}

	virtual ~CTRexData() {
		ASTRA_DELETE(m_pAlg);
		ASTRA_DELETE(m_pCfg);
		CData2DManager::getSingleton().remove(
			CData2DManager::getSingleton().getIndex(m_pProxInput));
		m_pProxInput = NULL;
	}

	virtual bool init(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda) = 0;
	virtual void prox(int32 iter) 
	{
		ASTRA_ASSERT(m_pAlg);
		// run for the required iterations
		m_pAlg->run(iter);
	}
};

// Least Squares data term
class _AstraExport CTRexDataLS : public CTRexData {
public:
	virtual bool init(const Config& cfg, CVolumeGeometry2D* geom,
					   float32 fLambda);
};

// Weighted least squares
class _AstraExport CTRexDataWLS : public CTRexData {
public:

};

// ----------------------------------------------------------------------------

/**
 * \brief
 * This class contains the implementation of the TRex algorithm for 2D reconstruction.
 *
 * The update step of pixel \f$v_j\f$ for projection \f$phi\f$ and iteration \f$k\f$ is given by:
 * \f[
 *	v_j^{(k+1)} = v_j^{(k)} + \frac{\sum_{p_i \in P_\phi} \left(  \lambda \frac{p_i - \sum_{r=1}^{N} w_{ir}v_r^{(k)}} {\sum_{r=1}^{N}w_{ir} }    \right)} {\sum_{p_i \in P_\phi}w_{ij}}
 * \f]
 *
 * \par XML Configuration
 * \astra_xml_item{ProjectorId, integer, Identifier of a projector as it is stored in the ProjectorManager.}
 * \astra_xml_item{ProjectionDataId, integer, Identifier of a projection data object as it is stored in the DataManager.}
 * \astra_xml_item{ReconstructionDataId, integer, Identifier of a volume data object as it is stored in the DataManager.}
 * \astra_xml_item_option{ReconstructionMaskId, integer, not used, Identifier of a volume data object that acts as a reconstruction mask. 1 = reconstruct on this pixel. 0 = don't reconstruct on this pixel.}
 * \astra_xml_item_option{SinogramMaskId, integer, not used, Identifier of a projection data object that acts as a projection mask. 1 = reconstruct using this ray. 0 = don't use this ray while reconstructing.}
 * \astra_xml_item_option{UseMinConstraint, bool, false, Use minimum value constraint.}
 * \astra_xml_item_option{MinConstraintValue, float, 0, Minimum constraint value.}
 * \astra_xml_item_option{UseMaxConstraint, bool, false, Use maximum value constraint.}
 * \astra_xml_item_option{MaxConstraintValue, float, 255, Maximum constraint value.}
 * \astra_xml_item_option{ProjectionOrder, string, "sequential", the order in which the projections are updated. 'sequential', 'random' or 'custom'}
 * \astra_xml_item_option{ProjectionOrderList, vector of float, not used, if ProjectionOrder='custom': use this order.}
 *
 * \par MATLAB example
 * \astra_code{
 *		cfg = astra_struct('SART');\n
 *		cfg.ProjectorId = proj_id;\n
 *		cfg.ProjectionDataId = sino_id;\n
 *		cfg.ReconstructionDataId = recon_id;\n
 *		cfg.option.MaskId = mask_id;\n
 *		cfg.option.UseMinConstraint = 'yes';\n 
 *		cfg.option.UseMaxConstraint = 'yes';\n
 *		cfg.option.MaxConstraintValue = 1024;\n
 *		cfg.option.ProjectionOrder = 'custom';\n
*		cfg.option.ProjectionOrderList = randperm(100);\n
 *		alg_id = astra_mex_algorithm('create'\, cfg);\n
 *		astra_mex_algorithm('iterate'\, alg_id\, 10);\n
 *		astra_mex_algorithm('delete'\, alg_id);\n
 * }
 */
class _AstraExport CTRexAlgorithm : public CReconstructionAlgorithm2D {
protected:

    /** Initial clearing. Only to be used by constructors.
     */
    virtual void _clear();

    /** Check the values of this object.  If everything is ok, the object can be set to the initialized state.
     * The following statements are then guaranteed to hold:
     * - valid projector
     * - valid data objects
     * - projection order all within range
     */
    virtual bool _check();

	// Data term: LS or WLS
	CTRexData* m_pData;

	// Prior term: ITV, ATV, SAD
	CTRexPrior* m_pPrior;

	// ADMM algorithm parameters: rho is the main parameter, 
	// and mu is the linearizaiton parameter.
	float32 m_fRho, m_fMu;

	// Regularization parameter that balances prior and data terms.
	float32 m_fSigma;

	// Weights for WLS that multiply the sinogram and the projection matrix.
	// Should be the sqrt of the WLS matrix.
	CFloat32ProjectionData2D* m_pW;

	// WLS nth root to apply before feeding to the prox operator for WLS data 
	// term
	float32 m_fWlsRoot;

	// Number of inner iterations
	int m_iInnterIter;

	// Inner iteration proximal operator metrics
	CFloat32VolumeData2D* m_pProxMetrics;

	//// The algorithm to solve the data term proximal operator.
	//CReconstructionAlgorithm2D* m_pDataProxOperator;

	//// Config that is passed to the tomography proximal operator.
	//Config* m_pTomoProxOperatorConfig;

public:
    
    // type of the algorithm, needed to register with CAlgorithmFactory
    static std::string type;	
    
    /** Default constructor, containing no code.
     */
    CTRexAlgorithm();
    
    /** Constructor.
     *
     * @param _pProjector		Projector Object.
     * @param _pSinogram		ProjectionData2D object containing the sinogram data.
     * @param _pReconstruction	VolumeData2D object for storing the reconstructed volume.
     */
    CTRexAlgorithm(CProjector2D* _pProjector, 
                   CFloat32ProjectionData2D* _pSinogram, 
                   CFloat32VolumeData2D* _pReconstruction);

    /** Destructor.
     */
    virtual ~CTRexAlgorithm();
    
    /** Clear this class.
     */
    virtual void clear();

    /** Initialize the algorithm with a config object.
     *
     * @param _cfg Configuration Object
     * @return initialization successful?
     */
    virtual bool initialize(const Config& _cfg);

    ///** Initialize class, no optionals, use sequential order.
    // *
    // * @param _pProjector		Projector Object.
    // * @param _pSinogram		ProjectionData2D object containing the sinogram data.
    // * @param _pReconstruction	VolumeData2D object for storing the reconstructed volume.
    // * @return initialization successful?	 
    // */
    //virtual bool initialize(CProjector2D* _pProjector, 
    //						CFloat32ProjectionData2D* _pSinogram,
    //						CFloat32VolumeData2D* _pReconstruction);

    ///** Initialize class, use custom order.
    // *
    // * @param _pProjector			Projector Object.
    // * @param _pSinogram			ProjectionData2D object containing the sinogram data.
    // * @param _pReconstruction		VolumeData2D object for storing the reconstructed volume.
    // * @param _piProjectionOrder	array containing a projection order.
    // * @param _iProjectionCount		number of elements in _piProjectionOrder.
    // * @return initialization successful?	 
    // */
    //virtual bool initialize(CProjector2D* _pProjector, 
    //						CFloat32ProjectionData2D* _pSinogram, 
    //						CFloat32VolumeData2D* _pReconstruction,
    //						int* _piProjectionOrder, 
    //						int _iProjectionCount);

    /** Get all information parameters
     *
     * @return map with all boost::any object
     */
    virtual map<string,boost::any> getInformation();

    /** Get a single piece of information represented as a boost::any
     *
     * @param _sIdentifier identifier string to specify which piece of information you want
     * @return boost::any object
     */
    virtual boost::any getInformation(std::string _sIdentifier);

    /** Perform a number of iterations.  Each iteration is a forward and backprojection of 
     * a single projection index.
     *
     * @param _iNrIterations amount of iterations to perform.
     */
    virtual void run(int _iNrIterations = 1);

    /** Get a description of the class.
     *
     * @return description string
     */
    virtual std::string description() const;

protected:


    ////< Order of the projections.
    //int* m_piProjectionOrder;
    ////< Number of projections specified in m_piProjectionOrder.
    //int m_iProjectionCount;
    ////< Current index in the projection order array.
    //int m_iCurrentProjection;

};

// inline functions
inline std::string CTRexAlgorithm::description() const { return CTRexAlgorithm::type; };


} // end namespace

#endif
