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

#include "astra/ArtProxOperatorAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"
#include "astra/Logging.h"


using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CArtProxOperatorAlgorithm::type = "ART-PROX";


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CArtProxOperatorAlgorithm::_clear()
{
	//CReconstructionAlgorithm2D::_clear();
	CArtAlgorithm::_clear();
	m_fLambda = 1.0f;
	m_pC = NULL;
	m_pY = NULL;

}

//---------------------------------------------------------------------------------------
// Clear - Public
void CArtProxOperatorAlgorithm::clear()
{
	//CReconstructionAlgorithm2D::clear();
	CArtAlgorithm::clear();
	//if (m_piProjectionOrder) {
	//	delete[] m_piProjectionOrder;
	//	m_piProjectionOrder = NULL;
	//}
	//m_iProjectionCount = 0;
	//m_iCurrentProjection = 0;
	//m_bIsInitialized = false;
	//m_iIterationCount = 0;
	ASTRA_DELETE(m_pY);
	ASTRA_DELETE(m_pC);
}

//----------------------------------------------------------------------------------------
// Constructor
CArtProxOperatorAlgorithm::CArtProxOperatorAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
CArtProxOperatorAlgorithm::CArtProxOperatorAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	_clear();
	//CArtAlgorithm::CArtAlgorithm(_pProjector, _pSinogram, _pReconstruction);
	//initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
//// Constructor
//CArtProxOperatorAlgorithm::CArtProxOperatorAlgorithm(CProjector2D* _pProjector, 
//							   CFloat32ProjectionData2D* _pSinogram, 
//							   CFloat32VolumeData2D* _pReconstruction,
//							   int* _piProjectionOrder, 
//							   int _iProjectionCount)
//{
//	_clear();
//	CArtAlgorithm::CArtAlgorithm(_pProjector, _pSinogram, _pReconstruction, 
//		_piProjectionOrder, _iProjectionCount);
//	//initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
//}

//----------------------------------------------------------------------------------------
// Destructor
CArtProxOperatorAlgorithm::~CArtProxOperatorAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CArtProxOperatorAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("BicavProxOperator", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class
	if (!CArtAlgorithm::initialize(_cfg)) {
		return false;
	}

	//// projection order
	//m_iCurrentProjection = 0;
	//m_iProjectionCount = m_pProjector->getProjectionGeometry()->getProjectionAngleCount();
	//string projOrder = _cfg.self.getOption("ProjectionOrder", "sequential");
	//CC.markOptionParsed("ProjectionOrder");
	//if (projOrder == "sequential") {
	//	m_piProjectionOrder = new int[m_iProjectionCount];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = i;
	//	}
	//} else if (projOrder == "random") {
	//	m_piProjectionOrder = new int[m_iProjectionCount];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = i;
	//	}
	//	for (int i = 0; i < m_iProjectionCount-1; i++) {
	//		int k = (rand() % (m_iProjectionCount - i));
	//		int t = m_piProjectionOrder[i];
	//		m_piProjectionOrder[i] = m_piProjectionOrder[i + k];
	//		m_piProjectionOrder[i + k] = t;
	//	}
	//} else if (projOrder == "custom") {
	//	vector<float32> projOrderList = _cfg.self.getOptionNumericalArray("ProjectionOrderList");
	//	m_piProjectionOrder = new int[projOrderList.size()];
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		m_piProjectionOrder[i] = static_cast<int>(projOrderList[i]);
	//	}
	//	CC.markOptionParsed("ProjectionOrderList");
	//}

	// Lambda
	m_fLambda = _cfg.self.getOptionNumerical("Lambda", m_fLambda);
	CC.markOptionParsed("Lambda");

	//// Alpha
	//m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	//CC.markOptionParsed("Alpha");

	//// Clear RaySum after each sweep. Defaults to true.
	//m_bClearRayLength = _cfg.self.getOptionBool("ClearRayLength", m_bClearRayLength);
	//CC.markOptionParsed("ClearRayLength");

	// Input volume
	XMLNode node = _cfg.self.getSingleNode("ProxInputDataId");
	ASTRA_CONFIG_CHECK(node, "ART-PROX", "No Proximal Input tag specified.");
	int id = boost::lexical_cast<int>(node.getContent());
	m_pProxInput = dynamic_cast<CFloat32VolumeData2D*>(CData2DManager::getSingleton().get(id));
	CC.markNodeParsed("ProxInputDataId");

	// create data objects
	//m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	//m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	//m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pY = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pC = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}


//----------------------------------------------------------------------------------------
bool CArtProxOperatorAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CArtAlgorithm::_check(), "ArtProxOperator", 
		"Error in ReconstructionAlgorithm2D initialization");

	return true;
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> CArtProxOperatorAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["ProjectionOrder"] = getInformation("ProjectionOrder");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any CArtProxOperatorAlgorithm::getInformation(std::string _sIdentifier) 
{
	//if (_sIdentifier == "ProjectionOrder") {
	//	vector<float32> res;
	//	for (int i = 0; i < m_iProjectionCount; i++) {
	//		res.push_back(m_piProjectionOrder[i]);
	//	}
	//	return res;
	//}
	return CAlgorithm::getInformation(_sIdentifier);
};

//----------------------------------------------------------------------------------------
// Iterate
void CArtProxOperatorAlgorithm::run(int _iNrIterations)
{
	// check initialized
	ASTRA_ASSERT(m_bIsInitialized);

	m_bShouldAbort = false;

	// data projectors
	CDataProjectorInterface* pForwardProjector;
	CDataProjectorInterface* pBackProjector;
	CDataProjectorInterface* pFirstForwardProjector;

	// Init Y = 0
	m_pY->setData(0.f);
	if (m_bClearReconstruction) {
		// Init Reconstruction = ProxInput
		m_pReconstruction->copyData(m_pProxInput->getData());
	} else {
		// Add Reconstruction to ProxInput.
		*m_pReconstruction += *m_pProxInput;
	}	
	//m_pReconstruction->updateStatistics();

	//ASTRA_INFO("Initialized from ProxInput: max=%f min=%f id=%d", 
	//	m_pReconstruction->getGlobalMax(), m_pReconstruction->getGlobalMin(),
	//	CData2DManager::getSingleton().getIndex(m_pReconstruction));
	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}	

	// constant
	float32 fSqrt2Lambda = sqrtf(m_fLambda * 2.f);
	//ASTRA_INFO("Sqrt2Lambda = %f", fSqrt2Lambda);
	//ASTRA_INFO("Alpha = %f", m_fAlpha);

	//ASTRA_INFO("UseMinConst=%d MinVal=%f", m_bUseMinConstraint, m_fMinValue);
	//ASTRA_INFO("UseMaxConst=%d MaxVal=%f", m_bUseMaxConstraint, m_fMaxValue);
	// Scale projections by sqrt(2 * lambda)
	// m_pSinogram->operator*=(fSqrt_2_lambda);

	// backprojection data projector
	pBackProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			SartProxBPPolicy(m_pReconstruction, m_pDiffSinogram, 
						 m_pTotalPixelWeight, m_pTotalRayLength, 
						 m_pY, m_pC, m_fAlpha, fSqrt2Lambda, true,
						 m_bUseMinConstraint, m_fMinValue, m_bUseMaxConstraint, m_fMaxValue),			// PSART backprojection
			m_bUseSinogramMask, m_bUseReconstructionMask, true // options on/off
		); 

	// first time forward projection data projector, unlike BICAV, don't need pixel weight
	pForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			DiffFPPolicy(m_pReconstruction, m_pDiffSinogram, m_pSinogram),								// forward projection with difference calculation
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// first time forward projection data projector,
	// computes total ray length, which is squared for ART.
	pFirstForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			TotalRayLengthPolicy(m_pTotalRayLength, true),													// calculate the total ray lengths
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// Set PixelWeight to 1, since it's unused here.
	m_pTotalPixelWeight->setData(1.0f);

	// Perform the first forward projection to compute ray lengths
	m_pTotalRayLength->setData(0.0f);
	pFirstForwardProjector->project();
	// Update with sqrt(2 lambda) because the SartProxBP multiplies by
	// teh other sqrt(2 lambda)
	// TODO: Make a new BicavProxBPPolicy instead of these workarounds.
	//*m_pTotalRayLength *= fSqrt2Lambda;

	// iteration loop, each iteration loops over all available projections
	for (int iIteration = 0; iIteration < _iNrIterations && !m_bShouldAbort; ++iIteration) {
		//ASTRA_INFO("Iteration %d", iIteration);
		// Clear RayLength before another loop over projections. This is needed so that
		// RayLength is correct, because updating RayLength with the forward projection
		// again will multiply the RayLength when processing the same ray in the next
		// iteration.
		//if (m_bClearRayLength) {
		//	m_pTotalRayLength->setData(0.f);
		//}

		// start timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount();

		// loop over rays
		for (int iR = 0; iR < m_iRayCount; ++iR) {
			// ray id
			int iProjection = m_piProjectionOrder[iR];
			int iDetector = m_piDetectorOrder[iR];
			//m_iCurrentRay = (m_iCurrentRay + 1) % m_iRayCount;
			//ASTRA_INFO(" Projection %d", iProjection);

			// forward projection and difference calculation
			pForwardProjector->projectSingleRay(iProjection, iDetector);
			// backprojection and clipping
			pBackProjector->projectSingleRay(iProjection, iDetector);

			//if (m_bUseMinConstraint)
			//	m_pReconstruction->clampMin(m_fMinValue);
			//if (m_bUseMaxConstraint)
			//	m_pReconstruction->clampMax(m_fMaxValue);
		}


		// end timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount() - m_ulTimer;

		//// Y residual
		//float32 res1 = m_pY->getNorm();
		//// current reconstruction - prox input
		//CFloat32Data2D* pResX = new CFloat32VolumeData2D(*m_pReconstruction);
		//*pResX -= *m_pProxInput;
		//float32 res2 = pResX->getNorm();
		//ASTRA_DELETE(pResX);
		//// Print residuals
		//ASTRA_INFO("norm(y) = %f & norm(z) = %f & lambda=%f", res1, res2, m_fLambda);

		// Compute metrics.
		computeIterationMetrics(iIteration, _iNrIterations);
	}

	ASTRA_DELETE(pForwardProjector);
	ASTRA_DELETE(pBackProjector);
	ASTRA_DELETE(pFirstForwardProjector);

	//for (int i=0; i < m_pReconstruction->getSize(); ++i) {
	//	ASTRA_INFO("voxel=%d val=%f", i, m_pReconstruction->getData()[i]);
	//}
}
//----------------------------------------------------------------------------------------

} // namespace astra
