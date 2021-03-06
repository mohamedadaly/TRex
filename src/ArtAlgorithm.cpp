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

#include "astra/ArtAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"


using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"

// type of the algorithm, needed to register with CAlgorithmFactory
std::string CArtAlgorithm::type = "ART";

//----------------------------------------------------------------------------------------
// Constructor
CArtAlgorithm::CArtAlgorithm() 
 : CReconstructionAlgorithm2D()
{
	m_fAlpha = 1.0f;
	m_iRayCount = 0;
	m_iCurrentRay = 0;
	m_piProjectionOrder = NULL;
	m_piDetectorOrder = NULL;
	m_bIsInitialized = false;
}

//----------------------------------------------------------------------------------------
// Destructor
CArtAlgorithm::~CArtAlgorithm() 
{
	if (m_piProjectionOrder != NULL)
		delete[] m_piProjectionOrder;
	if (m_piDetectorOrder != NULL)
		delete[] m_piDetectorOrder;
}

//---------------------------------------------------------------------------------------
// Clear - Constructors
void CArtAlgorithm::_clear()
{
	CReconstructionAlgorithm2D::_clear();
	m_piDetectorOrder = NULL;
	m_piProjectionOrder = NULL;
	m_iRayCount = 0;
	m_iCurrentRay = 0;
	m_bIsInitialized = false;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CArtAlgorithm::clear()
{
	CReconstructionAlgorithm2D::clear();
	if (m_piDetectorOrder) {
		delete[] m_piDetectorOrder;
		m_piDetectorOrder = NULL;
	}
	if (m_piProjectionOrder) {
		delete[] m_piProjectionOrder;
		m_piProjectionOrder = NULL;
	}
	m_fAlpha = 1.0f;
	m_iRayCount = 0;
	m_iCurrentRay = 0;
	m_bIsInitialized = false;
}

//---------------------------------------------------------------------------------------
// Check
bool CArtAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CReconstructionAlgorithm2D::_check(), "ART", "Error in ReconstructionAlgorithm2D initialization");

	// check ray order list
	for (int i = 0; i < m_iRayCount; i++) {
		if (m_piProjectionOrder[i] < 0 || m_piProjectionOrder[i] > m_pSinogram->getAngleCount()-1) {
			ASTRA_CONFIG_CHECK(false, "ART", "Invalid value in ray order list.");
		}
		if (m_piDetectorOrder[i] < 0 || m_piDetectorOrder[i] > m_pSinogram->getDetectorCount()-1) {
			ASTRA_CONFIG_CHECK(false, "ART", "Invalid value in ray order list.");
		}
	}

	// success
	return true;
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CArtAlgorithm::initialize(const Config& _cfg)
{
	ASTRA_ASSERT(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("ArtAlgorithm", this, _cfg);

	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}
	
	// initialization of parent class
	if (!CReconstructionAlgorithm2D::initialize(_cfg)) {
		return false;
	}

	// ray order
	string projOrder = _cfg.self.getOption("ProjectionOrder", "sequential");
	CC.markOptionParsed("ProjectionOrder");
	m_iCurrentRay = 0;
	m_iRayCount = m_pProjector->getProjectionGeometry()->getProjectionAngleCount() * 
		m_pProjector->getProjectionGeometry()->getDetectorCount();
	if (projOrder == "sequential") {
		m_piProjectionOrder = new int[m_iRayCount];
		m_piDetectorOrder = new int[m_iRayCount];
		for (int i = 0; i < m_iRayCount; i++) {
			m_piProjectionOrder[i] = (int)floor((float)i / m_pProjector->getProjectionGeometry()->getDetectorCount());
			m_piDetectorOrder[i] = i % m_pProjector->getProjectionGeometry()->getDetectorCount();
		}
	} else if (projOrder == "random") {
		srand(123);
		//
		// Shuffle the order of projections, and process rays in order in each projection
		//
		int iNumProj = m_pProjector->getProjectionGeometry()->getProjectionAngleCount();
		int* pRandProjOrder = new int[iNumProj];
		for (int i = 0; i < iNumProj; i++) {
			pRandProjOrder[i] = i;
		}
		// randomize projection order as in SART
		for (int i = 0; i < iNumProj - 1; i++) {
			int k = (rand() % (iNumProj - i));
			int t = pRandProjOrder[i];
			pRandProjOrder[i] = pRandProjOrder[i + k];
			pRandProjOrder[i + k] = t;
		}
		// fill in ray order
		int iRaysPerProj = m_pProjector->getProjectionGeometry()->getDetectorCount();
		m_piProjectionOrder = new int[m_iRayCount];
		m_piDetectorOrder = new int[m_iRayCount];
		for (int i = 0; i < iNumProj; ++i) {
			for (int d = 0; d < iRaysPerProj; ++d) {
				int iRayInd = i * iRaysPerProj + d;
				m_piProjectionOrder[iRayInd] = pRandProjOrder[i];
				m_piDetectorOrder[iRayInd] = d;
			}
		}
		delete [] pRandProjOrder;

		//
		// Shuffle the order of rays
		//
		//m_piProjectionOrder = new int[m_iRayCount];
		//m_piDetectorOrder = new int[m_iRayCount];
		//for (int i = 0; i < m_iRayCount; i++) {
		//	m_piProjectionOrder[i] = (int)floor((float)i / m_pProjector->getProjectionGeometry()->getDetectorCount());
		//	m_piDetectorOrder[i] = i % m_pProjector->getProjectionGeometry()->getDetectorCount();
		//}
		//// randomize
		//for (int i = 0; i < m_iRayCount; i++) {
		//	int k = (rand() % (m_iRayCount - i));
		//	int t1 = m_piProjectionOrder[i];
		//	m_piProjectionOrder[i] = m_piProjectionOrder[i + k];
		//	m_piProjectionOrder[i + k] = t1;

		//	int t2 = m_piDetectorOrder[i];
		//	m_piDetectorOrder[i] = m_piDetectorOrder[i + k];
		//	m_piDetectorOrder[i + k] = t2;
		//}
	} else if (projOrder == "custom") {
		vector<float32> rayOrderList = _cfg.self.getOptionNumericalArray("RayOrderList");
		m_iRayCount = rayOrderList.size() / 2;
		m_piProjectionOrder = new int[m_iRayCount];
		m_piDetectorOrder = new int[m_iRayCount];
		for (int i = 0; i < m_iRayCount; i++) {
			m_piProjectionOrder[i] = static_cast<int>(rayOrderList[2*i]);
			m_piDetectorOrder[i] = static_cast<int>(rayOrderList[2*i+1]);
		}
		CC.markOptionParsed("RayOrderList");
	} else {
		return false;
	}

	// Alpha
	m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	CC.markOptionParsed("Alpha");

	// create data objects
	m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//----------------------------------------------------------------------------------------
// Initialize - C++
bool CArtAlgorithm::initialize(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction)
{
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// required classes
	m_pProjector = _pProjector;
	m_pSinogram = _pSinogram;
	m_pReconstruction = _pReconstruction;

	// ray order
	m_iCurrentRay = 0;
	m_iRayCount = _pProjector->getProjectionGeometry()->getDetectorCount() * 
		_pProjector->getProjectionGeometry()->getProjectionAngleCount();
	m_piProjectionOrder = new int[m_iRayCount];
	m_piDetectorOrder = new int[m_iRayCount];
	for (int i = 0; i < m_iRayCount; i++) {
		m_piProjectionOrder[i] = (int)floor((float)i / _pProjector->getProjectionGeometry()->getDetectorCount());
		m_piDetectorOrder[i] = i % _pProjector->getProjectionGeometry()->getDetectorCount();
	}

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//----------------------------------------------------------------------------------------
// Set the relaxation factor.
void CArtAlgorithm::setLambda(float32 _fLambda)
{
	m_fAlpha = _fLambda;
}

//----------------------------------------------------------------------------------------
// Set the order in which the rays will be selected
void CArtAlgorithm::setRayOrder(int* _piProjectionOrder, int* _piDetectorOrder, int _iRayCount)
{
	if (m_piDetectorOrder) {
		delete[] m_piDetectorOrder;
		m_piDetectorOrder = NULL;
	}
	if (m_piProjectionOrder) {
		delete[] m_piProjectionOrder;
		m_piProjectionOrder = NULL;
	}

	m_iCurrentRay = 0;
	m_iRayCount = _iRayCount;
	m_piProjectionOrder = new int[m_iRayCount];
	m_piDetectorOrder = new int[m_iRayCount];
	for (int i = 0; i < m_iRayCount; i++) {
		m_piProjectionOrder[i] = _piProjectionOrder[i];
		m_piDetectorOrder[i] = _piDetectorOrder[i];
	}
}

//---------------------------------------------------------------------------------------
// Information - All
map<string,boost::any> CArtAlgorithm::getInformation() 
{
	map<string, boost::any> res;
	res["RayOrder"] = getInformation("RayOrder");
	res["Lambda"] = getInformation("Lambda");
	return mergeMap<string,boost::any>(CReconstructionAlgorithm2D::getInformation(), res);
};

//---------------------------------------------------------------------------------------
// Information - Specific
boost::any CArtAlgorithm::getInformation(std::string _sIdentifier) 
{
	if (_sIdentifier == "Lambda")	{ return m_fAlpha; }
	if (_sIdentifier == "RayOrder") {
		vector<float32> res;
		for (int i = 0; i < m_iRayCount; i++) {
			res.push_back(m_piProjectionOrder[i]);
		}
		for (int i = 0; i < m_iRayCount; i++) {
			res.push_back(m_piDetectorOrder[i]);
		}
		return res;
	}
	return CAlgorithm::getInformation(_sIdentifier);
};

//----------------------------------------------------------------------------------------
// Iterate
void CArtAlgorithm::run(int _iNrIterations)
{
	// check initialized
	ASTRA_ASSERT(m_bIsInitialized);

	m_bShouldAbort = false;

	// data projectors
	CDataProjectorInterface* pForwardProjector;
	CDataProjectorInterface* pBackProjector;
	CDataProjectorInterface* pFirstForwardProjector;

	m_pTotalRayLength->setData(0.0f);

	// This is initialized to 1.0 to be able to use the SIRTBPPolicy.
	m_pTotalPixelWeight->setData(1.0f);

	// Initialize m_pReconstruction to zero.
	if (m_bClearReconstruction) {
		m_pReconstruction->setData(0.f);
	}

	// forward projection data projector
	pForwardProjector = dispatchDataProjector(
		m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			DiffFPPolicy(m_pReconstruction, m_pDiffSinogram, m_pSinogram),								// forward projection with difference calculation
			m_bUseSinogramMask, m_bUseReconstructionMask, true											// options on/off
		); 

	// backprojection data projector
	pBackProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			SIRTBPPolicy(m_pReconstruction, m_pDiffSinogram, 
				m_pTotalPixelWeight, m_pTotalRayLength, m_fAlpha, NULL,
				m_bUseMinConstraint, m_fMinValue, m_bUseMaxConstraint, m_fMaxValue),  // SIRT backprojection
			m_bUseSinogramMask, m_bUseReconstructionMask, true // options on/off
		); 

	// first time forward projection data projector,
	// also computes total pixel weight and total ray length
	pFirstForwardProjector = dispatchDataProjector(
			m_pProjector, 
			SinogramMaskPolicy(m_pSinogramMask),														// sinogram mask
			ReconstructionMaskPolicy(m_pReconstructionMask),											// reconstruction mask
			TotalRayLengthPolicy(m_pTotalRayLength, true),											    // calculate the total ray lengths squared
			m_bUseSinogramMask, m_bUseReconstructionMask, true											 // options on/off
		);

	// Perform the first forward projection to compute ray lengths and pixel weights
	pFirstForwardProjector->project();

	// iteration loop, each iteration loops over all available rays
	for (int iIteration = 0; iIteration < _iNrIterations && !m_bShouldAbort; ++iIteration) {
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
			// backprojection
			pBackProjector->projectSingleRay(iProjection, iDetector);

			//if (m_bUseMinConstraint)
			//	m_pReconstruction->clampMin(m_fMinValue);
			//if (m_bUseMaxConstraint)
			//	m_pReconstruction->clampMax(m_fMaxValue);
		}

		// end timer
		m_ulTimer = CPlatformDepSystemCode::getMSCount() - m_ulTimer;

		// Compute metrics.
		computeIterationMetrics(iIteration, _iNrIterations, m_pDiffSinogram);
	}

	ASTRA_DELETE(pForwardProjector);
	ASTRA_DELETE(pBackProjector);
	ASTRA_DELETE(pFirstForwardProjector);

	//// check initialized
	//assert(m_bIsInitialized);
	//
	//// variables
	//int iIteration, iPixel;
	//int iUsedPixels, iProjection, iDetector;
	//float32 fRayForwardProj, fSumSquaredWeights;
	//float32 fProjectionDifference, fBackProjectionFactor;

	//// create a pixel buffer
	//int iPixelBufferSize = m_pProjector->getProjectionWeightsCount(0);
	//SPixelWeight* pPixels = new SPixelWeight[iPixelBufferSize];

	//// start iterations
	//for (iIteration = _iNrIterations-1; iIteration >= 0; --iIteration) {

	//	// step0: compute single weight rays
	//	iProjection = m_piProjectionOrder[m_iCurrentRay];
	//	iDetector = m_piDetectorOrder[m_iCurrentRay];
	//	m_iCurrentRay = (m_iCurrentRay + 1) % m_iRayCount;

	//	if (m_bUseSinogramMask && m_pSinogramMask->getData2D()[iProjection][iDetector] == 0) continue;	

	//	m_pProjector->computeSingleRayWeights(iProjection, iDetector, pPixels, iPixelBufferSize, iUsedPixels);

	//	// step1: forward projections
	//	fRayForwardProj = 0.0f;
	//	fSumSquaredWeights = 0.0f;
	//	for (iPixel = iUsedPixels-1; iPixel >= 0; --iPixel) {
	//		if (m_bUseReconstructionMask && m_pReconstructionMask->getDataConst()[pPixels[iPixel].m_iIndex] == 0) continue;

	//		fRayForwardProj += pPixels[iPixel].m_fWeight * m_pReconstruction->getDataConst()[pPixels[iPixel].m_iIndex];
	//		fSumSquaredWeights += pPixels[iPixel].m_fWeight * pPixels[iPixel].m_fWeight;
	//	}
	//	if (fSumSquaredWeights == 0) continue;

	//	// step2: difference
	//	fProjectionDifference = m_pSinogram->getData2D()[iProjection][iDetector] - fRayForwardProj;

	//	// step3: back projection
	//	fBackProjectionFactor = m_fAlpha * fProjectionDifference / fSumSquaredWeights;
	//	for (iPixel = iUsedPixels-1; iPixel >= 0; --iPixel) {
	//		
	//		// pixel must be loose
	//		if (m_bUseReconstructionMask && m_pReconstructionMask->getDataConst()[pPixels[iPixel].m_iIndex] == 0) continue;

	//		// update
	//		m_pReconstruction->getData()[pPixels[iPixel].m_iIndex] += fBackProjectionFactor * pPixels[iPixel].m_fWeight;
	//		
	//		// constraints
	//		if (m_bUseMinConstraint && m_pReconstruction->getData()[pPixels[iPixel].m_iIndex] < m_fMinValue) {
	//			m_pReconstruction->getData()[pPixels[iPixel].m_iIndex] = m_fMinValue;
	//		}
	//		if (m_bUseMaxConstraint && m_pReconstruction->getData()[pPixels[iPixel].m_iIndex] > m_fMaxValue) {
	//			m_pReconstruction->getData()[pPixels[iPixel].m_iIndex] = m_fMaxValue;
	//		}
	//	}

	//}
	//delete[] pPixels;

	//// update statistics
	//m_pReconstruction->updateStatistics();
}


//----------------------------------------------------------------------------------------

} // namespace astra
