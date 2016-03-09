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

#include "astra/ArtBaseAlgorithm.h"

#include <boost/lexical_cast.hpp>

#include "astra/AstraObjectManager.h"
#include "astra/DataProjectorPolicies.h"

using namespace std;

namespace astra {

#include "astra/Projector2DImpl.inl"


//---------------------------------------------------------------------------------------
// Clear - Constructors
void CArtBaseAlgorithm::_clear()
{
	CReconstructionAlgorithm2D::_clear();
	m_bIsInitialized = false;
	m_iIterationCount = 0;
	m_fAlpha = 1.0f;
	m_bClearRayLength = true;
}

//---------------------------------------------------------------------------------------
// Clear - Public
void CArtBaseAlgorithm::clear()
{
	CReconstructionAlgorithm2D::clear();
	m_bIsInitialized = false;
	m_iIterationCount = 0;
}

//----------------------------------------------------------------------------------------
// Constructor
CArtBaseAlgorithm::CArtBaseAlgorithm() 
{
	_clear();
}

//----------------------------------------------------------------------------------------
// Constructor
CArtBaseAlgorithm::CArtBaseAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction) 
{
	_clear();
	initialize(_pProjector, _pSinogram, _pReconstruction);
}

//----------------------------------------------------------------------------------------
// Constructor
CArtBaseAlgorithm::CArtBaseAlgorithm(CProjector2D* _pProjector, 
							   CFloat32ProjectionData2D* _pSinogram, 
							   CFloat32VolumeData2D* _pReconstruction,
							   int* _piProjectionOrder, 
							   int _iProjectionCount)
{
	_clear();
	initialize(_pProjector, _pSinogram, _pReconstruction, _piProjectionOrder, _iProjectionCount);
}

//----------------------------------------------------------------------------------------
// Destructor
CArtBaseAlgorithm::~CArtBaseAlgorithm() 
{
	clear();
}

//---------------------------------------------------------------------------------------
// Initialize - Config
bool CArtBaseAlgorithm::initialize(const Config& _cfg)
{
	assert(_cfg.self);
	ConfigStackCheck<CAlgorithm> CC("ArtBaseAlgorithm", this, _cfg);
	
	// if already initialized, clear first
	if (m_bIsInitialized) {
		clear();
	}

	// initialization of parent class
	if (!CReconstructionAlgorithm2D::initialize(_cfg)) {
		return false;
	}

	// Alpha
	m_fAlpha = _cfg.self.getOptionNumerical("Alpha", m_fAlpha);
	CC.markOptionParsed("Alpha");

	// Clear RaySum after each sweep. Defaults to true.
	m_bClearRayLength = _cfg.self.getOptionBool("ClearRayLength", m_bClearRayLength);
	CC.markOptionParsed("ClearRayLength");

	// create data objects
	m_pTotalRayLength = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());
	m_pTotalPixelWeight = new CFloat32VolumeData2D(m_pProjector->getVolumeGeometry());
	m_pDiffSinogram = new CFloat32ProjectionData2D(m_pProjector->getProjectionGeometry());

	// success
	m_bIsInitialized = _check();
	return m_bIsInitialized;
}

//----------------------------------------------------------------------------------------
bool CArtBaseAlgorithm::_check()
{
	// check base class
	ASTRA_CONFIG_CHECK(CReconstructionAlgorithm2D::_check(), "ART Base", "Error in ReconstructionAlgorithm2D initialization");

	return true;
}


} // namespace astra
