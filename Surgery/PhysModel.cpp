#include "PhysModel.h"
#include "collisions/CGenericCollision.h"
#include "collisions/CCollisionBrute.h"
#include "collisions/CCollisionAABB.h"
#include "shaders/CShaderProgram.h"

#include <vector_functions.hpp>

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>

#include <thrust\sort.h>
#include <thrust/execution_policy.h>

#include "GL/glu.h"


PhysModel::PhysModel(std::string filename)
{
	// create vertex array, only use color data!
	vertices = cVertexArray::create(false, false, true, false, false, false);

	// create triangle array,指针传递
	points = cPointArray::create(vertices);

	// should the frame (X-Y-Z) be displayed?
	m_showFrame = false;

	// set default collision detector
	m_collisionDetector = NULL;

	// display lists disabled by default
	m_useDisplayList = false;

	// default point size
	m_pointSize = 1.0;

	// show points
	m_showPoints = true;

	//是否把点画成球
	isDrawAsSphere = false;

	//球半径（用于渲染）
	sphereRadius = 0.01;

	//cuda
	isRegisterBuffer = false;

	//need find neighbors
	needUpdateNeighbors = true;

	//------------------------------------
	//Physical init
	//------------------------------------
	//Load point cloud
	LoadFromASC(filename);

	//find extreme particle pos
	FindExtremePoint();

	//particles count
	numParticles = vertices->m_localPos.size();

	//every particle has neighbors(26 * particles count)
	hNeighbors.resize(26 * numParticles, -1);
	hNeighborCount.resize(numParticles, 0);

	//Init hPos & hVel & hPredictPos
	double3 zero = make_double3(0, 0, 0);
	hPos.resize(numParticles, zero);
	hVel.resize(numParticles, zero);
	hPredictPos.resize(numParticles, zero);

	//Init hGridParticle Index & hash
	hGridParticleHash.resize(numParticles, 0);
	hGridParticleIndex.resize(numParticles, 0);

	//init hSortedPos & vel
	hSortedPos.resize(numParticles);
	hSortedVel.resize(numParticles);

	//Init hNearbyIndex & hNearbyCount
	hNearbyIdx.resize(2000, -1);
	hNearbyCount = 0;

	//Init hInvMass
	hInvMass.resize(numParticles, 1.0);

	//phys param setting
	physParam.Solid_1_EffectiveRadius = 0.06;

	physParam.gridSize.x = physParam.gridSize.y = physParam.gridSize.z = 64;
	physParam.numCells = physParam.gridSize.x * physParam.gridSize.y * physParam.gridSize.z;

	//init cellstart & end
	hCellStart.resize(physParam.numCells, -1);
	hCellEnd.resize(physParam.numCells, -1);

	double cellSize = physParam.Solid_1_EffectiveRadius * 2;
	physParam.cellSize = make_double3(cellSize, cellSize, cellSize);
	
	physParam.numBodies = numParticles;
	physParam.worldOrigin = make_double3(extremePoint_xMinus.x(), extremePoint_yMinus.y(), extremePoint_zMinus.z());

	//copy velocity to GPU
	AllocateArray((void**)&dVel, sizeof(double3) * hVel.size());
	CopyArrayToDevice(dVel, hVel.data(), 0, sizeof(double3) * hVel.size());

	//copy neighbors to GPU
	AllocateArray((void**)&dNeighbors, sizeof(int) * hNeighbors.size());
	CopyArrayToDevice(dNeighbors, hNeighbors.data(), 0, sizeof(int) * hNeighbors.size());

	//copy neighbor count to GPU
	AllocateArray((void**)&dNeighborCount, sizeof(uint) * hNeighborCount.size());
	CopyArrayToDevice(dNeighborCount, hNeighborCount.data(), 0, sizeof(uint) * hNeighborCount.size());

	//Allocate memsize for grid particle hash
	AllocateArray((void**)&dGridParticleHash, numParticles * sizeof(uint));

	//Allocate memsize for grid particle index
	AllocateArray((void**)&dGridParticleIndex, numParticles * sizeof(uint));

	//Allocate memsize for sorted pos & vel
	AllocateArray((void**)&dSortedPos, numParticles * sizeof(double3));
	AllocateArray((void**)&dSortedVel, numParticles * sizeof(double3));

	//Allocate memsize for cell start & end
	AllocateArray((void**)&dCellStart, physParam.numCells * sizeof(uint));
	AllocateArray((void**)&dCellEnd, physParam.numCells * sizeof(uint));

	//Allocate memsize for near tool particles indices
	AllocateArray((void**)&dNearbyIdx, sizeof(int) * hNearbyIdx.size());
	CopyArrayToDevice(dNearbyIdx, hNearbyIdx.data(), 0, sizeof(int) * hNearbyIdx.size());

	//Allocate memsize for near tool particles count
	AllocateArray((void**)&dNearbyCount, sizeof(uint));
	CopyArrayToDevice(dNearbyCount, &hNearbyCount, 0, sizeof(uint));

	//Allocate memsize for predict pos
	AllocateArray((void**)&dPredictPos, hPredictPos.size() * sizeof(double3));
	CopyArrayToDevice(dPredictPos, hPredictPos.data(), 0, sizeof(double3) * hPredictPos.size());

	//Allocate memsize for dInvMass
	AllocateArray((void**)&dInvMass, hInvMass.size() * sizeof(double));
	CopyArrayToDevice(dInvMass, hInvMass.data(), 0, sizeof(double) * hInvMass.size());

	//copy physParam to GPU
	SetParameters(&physParam);
}


PhysModel::~PhysModel()
{
	m_displayList.invalidate();

	FreeMemory();
}

//释放与物理属性、cuda相关内存
void PhysModel::FreeMemory()
{
	hPos.clear();
	hVel.clear();
	hNeighbors.clear();
	hNeighborCount.clear();
	hSortedPos.clear();
	hSortedVel.clear();
	hGridParticleHash.clear();
	hGridParticleIndex.clear();
	hCellStart.clear();
	hCellEnd.clear();
	hPredictPos.clear();
	hInvMass.clear();
	hNearbyIdx.clear();

	FreeArray(dVel);
	FreeArray(dSortedPos);
	FreeArray(dSortedVel);
	FreeArray(dCellStart);
	FreeArray(dCellEnd);
	FreeArray(dGridParticleIndex);
	FreeArray(dGridParticleHash);

	FreeArray(dNearbyIdx);
	FreeArray(dNearbyCount);

	FreeArray(dPredictPos);
	FreeArray(dInvMass);

	UnRegisterGLBufferObject(vertsPos_resource);
}

bool PhysModel::LoadFromASC(std::string filename)
{
	int pointCount = 0;
	if (!ReadASC(filename.c_str(), vertices, pointCount))
		return false;

	points->m_vertices = vertices;

	for (unsigned i = 0; i < vertices->getNumElements(); i++)
	{
		NewPoint(i);
	}

	return true;
}

void PhysModel::MarkForUpdate(const bool a_affectChildren)
{
	points->m_flagMarkForUpdate = true;

	cGenericObject::markForUpdate(a_affectChildren);
}

//--------------------------------------------------------------------------
// PUBLIC METHODS - VERTICES
//--------------------------------------------------------------------------
unsigned PhysModel::NewVertex(const double a_x, const double a_y, const double a_z)
{
	unsigned index = vertices->newVertex();

	vertices->setLocalPos(index, a_x, a_y, a_z);

	return index;
}

//! This method creates a new vertex and adds it to the vertex list.
unsigned PhysModel::NewVertex(const cVector3d& a_pos, const cColorf& a_color)
{
	unsigned index = vertices->newVertex();

	// set data
	vertices->setLocalPos(index, a_pos);
	vertices->setColor(index, a_color);

	// return vertex index
	return (index);
}

//--------------------------------------------------------------------------
// PUBLIC METHODS - POINTS
//--------------------------------------------------------------------------
unsigned PhysModel::NewPoint(const unsigned a_indexVertex0)
{
	int index = points->newPoint(a_indexVertex0);

	MarkForUpdate(true);

	return index;
}

unsigned PhysModel::NewPoint(const cVector3d& a_vertex0, const cColorf& a_colorVertex0)
{
	unsigned indexVertex0 = vertices->newVertex();

	int index = points->newPoint(indexVertex0);
	vertices->setLocalPos(indexVertex0, a_vertex0);
	vertices->setColor(indexVertex0, a_colorVertex0);

	// mark object for update
	MarkForUpdate(false);

	// return result
	return (index);
}

bool PhysModel::RemovePoint(const unsigned int a_index)
{
	points->removePoint(a_index);

	// mark object for update
	MarkForUpdate(false);

	// return success
	return (true);
}

unsigned PhysModel::GetNumPoints()
{
	return (unsigned)(points->getNumElements());
}

void PhysModel::Clear()
{
	// clear all triangles
	points->clear();

	// clear all vertices
	vertices->clear();
}

//-----------------------------------------------------------------------
// PUBLIC METHODS - GRAPHIC PROPERTIES:
//-----------------------------------------------------------------------

void PhysModel::setTransparencyLevel(const float a_level, const bool a_applyToVertices, const bool a_applyToTextures, const bool a_affectChildren)
{
	// if the transparency level is equal to 1.0, then do not apply transparency
	// otherwise enable it.
	if (a_level < 1.0)
	{
		setUseTransparency(true);
	}
	else
	{
		setUseTransparency(false);
	}

	// apply new value to material
	if (m_material != nullptr)
	{
		m_material->setTransparencyLevel(a_level);
	}

	// apply new value to texture
	if (m_texture != nullptr)
	{
		if (m_texture->m_image != nullptr)
		{
			unsigned char level = (unsigned char)(255.0 * a_level);
			m_texture->m_image->setTransparency(level);
		}
	}

	// apply change to children
	if (a_affectChildren)
	{
		std::vector<cGenericObject*>::iterator it;
		for (it = m_children.begin(); it < m_children.end(); it++)
		{
			(*it)->setTransparencyLevel(a_level, a_applyToVertices, a_applyToTextures, true);
		}
	}
}

//--------------------------------------------------------------------------
// PUBLIC METHODS - POINTS GRAPHIC PROPERTIES
//--------------------------------------------------------------------------

void PhysModel::SetPointColor(const cColorf& a_color)
{
	// apply color to all vertex colors
	int numVertices = vertices->getNumElements();

	for (int i = 0; i<numVertices; i++)
	{
		vertices->setColor(i, a_color);
	}
}

//--------------------------------------------------------------------------
// PUBLIC METHODS - COLLISION DETECTION:
//--------------------------------------------------------------------------

void PhysModel::createBruteForceCollisionDetector()
{
	// delete previous collision detector
	if (m_collisionDetector != NULL)
	{
		delete m_collisionDetector;
		m_collisionDetector = NULL;
	}

	// create brute collision detector
	m_collisionDetector = new cCollisionBrute(points);
}


void PhysModel::createAABBCollisionDetector(const double a_radius)
{
	// delete previous collision detector
	if (m_collisionDetector != NULL)
	{
		delete m_collisionDetector;
		m_collisionDetector = NULL;
	}

	// create AABB collision detector
	cCollisionAABB* collisionDetector = new cCollisionAABB();
	collisionDetector->initialize(points, a_radius);

	// assign new collision detector
	m_collisionDetector = collisionDetector;
}

//--------------------------------------------------------------------------
// PUBLIC METHODS - GEOMETRY:
//--------------------------------------------------------------------------

//! This method scales this point cloud by using different scale factors along X, Y and Z axes.
void PhysModel::ScaleXYZ(const double a_scaleX, const double a_scaleY, const double a_scaleZ)
{
	int numVertices = vertices->getNumElements();

	for (int i = 0; i<numVertices; i++)
	{
		vertices->m_localPos[i].mul(a_scaleX, a_scaleY, a_scaleZ);
	}

	m_boundaryBoxMax.mul(a_scaleX, a_scaleY, a_scaleZ);
	m_boundaryBoxMin.mul(a_scaleX, a_scaleY, a_scaleZ);
}

//! This method shifts all vertex positions by the specified amount.
void PhysModel::offsetVertices(const cVector3d& a_offset, const bool a_updateCollisionDetector)
{
	// offset all vertices
	int numVertices = vertices->getNumElements();

	for (int i = 0; i<numVertices; i++)
	{
		vertices->m_localPos[i].add(a_offset);
	}

	// update boundary box
	m_boundaryBoxMin += a_offset;
	m_boundaryBoxMax += a_offset;
}

//! This method computes the center of mass of this point cloud, based on vertex positions.
cVector3d PhysModel::getCenterOfMass()
{
	cVector3d centerOfMass(0, 0, 0);
	int numVertices = vertices->getNumElements();

	if (numVertices > 0)
	{
		for (int i = 0; i<numVertices; i++)
		{
			cVector3d pos = vertices->getLocalPos(i);
			centerOfMass += pos;
		}
		centerOfMass.mul(1.0 / numVertices);
	}

	return (centerOfMass);
}


//--------------------------------------------------------------------------
// PROTECTED METHODS:
//--------------------------------------------------------------------------

//! This method renders this object graphically using OpenGL.
void PhysModel::render(cRenderOptions& a_options)
{
	/////////////////////////////////////////////////////////////////////////
	// Render parts that are always opaque
	/////////////////////////////////////////////////////////////////////////
	if (SECTION_RENDER_OPAQUE_PARTS_ONLY(a_options))
	{
		// render points
		if (m_showPoints && !isDrawAsSphere)
		{
			renderPoints(a_options);
		}
		else if (m_showPoints && isDrawAsSphere)
		{
			RenderSpheres(a_options);
		}
	}

	/////////////////////////////////////////////////////////////////////////
	// Render parts that use material properties
	/////////////////////////////////////////////////////////////////////////
	if (SECTION_RENDER_PARTS_WITH_MATERIALS(a_options, m_useTransparency))
	{
	}
}

void PhysModel::RenderSpheres(cRenderOptions& a_options)
{
#ifdef C_USE_OPENGL
	//--------------------------------------------------------------------------
	// INITIALIZATION
	//--------------------------------------------------------------------------

	unsigned int numVertices = vertices->getNumElements();
	unsigned int numPoints = points->getNumElements();

	// check if object contains any points or vertices
	if ((numVertices == 0) || (numPoints == 0))
	{
		return;
	}

	// initialize some OpenGL flags
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_INDEX_ARRAY);
	glDisableClientState(GL_EDGE_FLAG_ARRAY);
	glDisable(GL_COLOR_MATERIAL);

	//--------------------------------------------------------------------------
	// RENDER MATERIAL
	//--------------------------------------------------------------------------

	// render material properties if enabled
	if (m_useMaterialProperty && a_options.m_render_materials)
	{
		m_material->render(a_options);
	}

	//--------------------------------------------------------------------------
	// RENDER TEXTURE
	//--------------------------------------------------------------------------

	// check if texture is available
	if (m_texture == nullptr)
	{
		m_useTextureMapping = false;
	}

	// render texture if enabled
	if ((m_texture != nullptr) && (m_useTextureMapping) && (a_options.m_render_materials))
	{
		m_texture->renderInitialize(a_options);
	}

	//--------------------------------------------------------------------------
	// SETTINGS
	//--------------------------------------------------------------------------

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);

	//--------------------------------------------------------------------------
	// RENDER OBJECT (OBJECT BUFFER + SHADERS)
	//--------------------------------------------------------------------------
	if ((m_shaderProgram != nullptr) && (!a_options.m_creating_shadow_map))
	{
		// enable shader
		m_shaderProgram->use(this, a_options);

		// render normal texture if enabled
		if (m_normalMap != nullptr)
		{
			if (m_shaderProgram->isUsed())
			{
				m_normalMap->renderInitialize(a_options);
			}
		}

		// render VBO
		points->render();

		// disable normal map is available
		if (m_normalMap != nullptr)
		{
			m_normalMap->renderFinalize(a_options);
		}

		// disable shader
		m_shaderProgram->disable();
	}
	//--------------------------------------------------------------------------
	// RENDER OBJECT (OLD METHOD)
	//--------------------------------------------------------------------------
	else if (!m_displayList.render(m_useDisplayList))
	{
		// get texture unit
		GLenum textureUnit = GL_TEXTURE1_ARB;
		if (m_texture != nullptr)
		{
			textureUnit = m_texture->getTextureUnit();
		}

		// if requested, begin creating display list
		m_displayList.begin(m_useDisplayList);

		//-------------------------------------------------------------------
		// RENDER ALL Spheres
		//-------------------------------------------------------------------

		/////////////////////////////////////////////////////////////////////
		// RENDER SEGMENTS USING CLASSIC OPENGL COMMANDS
		/////////////////////////////////////////////////////////////////////
		{
			// begin rendering
			for (unsigned i = 0; i < vertices->m_localPos.size(); i++)
			{
				glPushMatrix();
				cVector3d tmpPos = vertices->m_localPos[i];
				glTranslated(tmpPos.x(), tmpPos.y(), tmpPos.z());
				cDrawSphere(sphereRadius, 5, 5);
				glPopMatrix();
			}
			
		}

		//-------------------------------------------------------------------
		// FINALIZE DISPLAY LIST
		//-------------------------------------------------------------------

		// if being created, finalize display list
		m_displayList.end(true);
	}

	//--------------------------------------------------------------------------
	// RESTORE OPENGL
	//--------------------------------------------------------------------------

	// turn off texture rendering
	if ((m_texture != nullptr) && (m_useTextureMapping))
	{
		m_texture->renderFinalize(a_options);
	}

	// restore OpenGL variables
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
#endif
}

void PhysModel::renderPoints(cRenderOptions& a_options)
{
#ifdef C_USE_OPENGL

	//--------------------------------------------------------------------------
	// INITIALIZATION
	//--------------------------------------------------------------------------

	unsigned int numVertices = vertices->getNumElements();
	unsigned int numPoints = points->getNumElements();

	// check if object contains any points or vertices
	if ((numVertices == 0) || (numPoints == 0))
	{
		return;
	}

	// initialize some OpenGL flags
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_INDEX_ARRAY);
	glDisableClientState(GL_EDGE_FLAG_ARRAY);
	glDisable(GL_COLOR_MATERIAL);

	//--------------------------------------------------------------------------
	// RENDER MATERIAL
	//--------------------------------------------------------------------------

	// render material properties if enabled
	if (m_useMaterialProperty && a_options.m_render_materials)
	{
		m_material->render(a_options);
	}

	//--------------------------------------------------------------------------
	// RENDER TEXTURE
	//--------------------------------------------------------------------------

	// check if texture is available
	if (m_texture == nullptr)
	{
		m_useTextureMapping = false;
	}

	// render texture if enabled
	if ((m_texture != nullptr) && (m_useTextureMapping) && (a_options.m_render_materials))
	{
		m_texture->renderInitialize(a_options);
	}

	//--------------------------------------------------------------------------
	// SETTINGS
	//--------------------------------------------------------------------------

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);

	//--------------------------------------------------------------------------
	// RENDER OBJECT (OBJECT BUFFER + SHADERS)
	//--------------------------------------------------------------------------
	if ((m_shaderProgram != nullptr) && (!a_options.m_creating_shadow_map))
	{
		// enable shader
		m_shaderProgram->use(this, a_options);

		// render normal texture if enabled
		if (m_normalMap != nullptr)
		{
			if (m_shaderProgram->isUsed())
			{
				m_normalMap->renderInitialize(a_options);
			}
		}

		// render VBO
		points->render();//存在无效操作！已解决

		//CUDA: register buffer for first time
		if (!isRegisterBuffer)
		{
			RegisterCudaBuffer(&vertsPos_resource, vertices->m_positionBuffer);

			isRegisterBuffer = true;
		}

		// disable normal map is available
		if (m_normalMap != nullptr)
		{
			m_normalMap->renderFinalize(a_options);
		}

		// disable shader
		m_shaderProgram->disable();
	}
	//--------------------------------------------------------------------------
	// RENDER OBJECT (OLD METHOD)
	//--------------------------------------------------------------------------
	else if (!m_displayList.render(m_useDisplayList))
	{
		// get texture unit
		GLenum textureUnit = GL_TEXTURE1_ARB;
		if (m_texture != nullptr)
		{
			textureUnit = m_texture->getTextureUnit();
		}
	
		// if requested, begin creating display list
		m_displayList.begin(m_useDisplayList);

		//-------------------------------------------------------------------
		// RENDER ALL POINTS
		//-------------------------------------------------------------------

		// set point size
		glPointSize((GLfloat)m_pointSize);

		/////////////////////////////////////////////////////////////////////
		// RENDER SEGMENTS USING CLASSIC OPENGL COMMANDS
		/////////////////////////////////////////////////////////////////////
		{
			// begin rendering
			glBegin(GL_POINTS);

			for (unsigned int i = 0; i<numPoints; i++)
			{
				if (points->m_allocated[i])
				{
					unsigned int index0 = points->getVertexIndex0(i);

					// render vertex 0
					glColor4fv(vertices->m_color[index0].getData());
					glVertex3dv(&vertices->m_localPos[index0](0));
				}
			}

			// finalize rendering 
			glEnd();
		}

		//-------------------------------------------------------------------
		// FINALIZE DISPLAY LIST
		//-------------------------------------------------------------------

		// if being created, finalize display list
		m_displayList.end(true);
	}

	//--------------------------------------------------------------------------
	// RESTORE OPENGL
	//--------------------------------------------------------------------------

	// turn off texture rendering
	if ((m_texture != nullptr) && (m_useTextureMapping))
	{
		m_texture->renderFinalize(a_options);
	}

	// restore OpenGL variables
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);

#endif
}

//! This method updates the global position of each vertex.
void PhysModel::updateGlobalPositions(const bool a_frameOnly)
{
	if (a_frameOnly) return;

	int numVertices = vertices->getNumElements();

	for (int i = 0; i<numVertices; i++)
	{
		vertices->computeGlobalPosition(i, m_globalPos, m_globalRot);
	}
}

//! This method updates the boundary box dimensions based on the vertices.
void PhysModel::updateBoundaryBox()
{
	unsigned int numPoints = points->getNumElements();

	if (numPoints == 0)
	{
		m_boundaryBoxMin.zero();
		m_boundaryBoxMax.zero();
		m_boundaryBoxEmpty = true;
		return;
	}

	double xMin = C_LARGE;
	double yMin = C_LARGE;
	double zMin = C_LARGE;
	double xMax = -C_LARGE;
	double yMax = -C_LARGE;
	double zMax = -C_LARGE;
	bool flag = false;

	// loop over all points
	for (unsigned int i = 0; i<numPoints; i++)
	{
		if (points->m_allocated[i])
		{
			cVector3d tVertex0 = vertices->getLocalPos(points->getVertexIndex0(i));
			xMin = cMin(tVertex0(0), xMin);
			yMin = cMin(tVertex0(1), yMin);
			zMin = cMin(tVertex0(2), zMin);
			xMax = cMax(tVertex0(0), xMax);
			yMax = cMax(tVertex0(1), yMax);
			zMax = cMax(tVertex0(2), zMax);

			flag = true;
		}
	}

	if (flag)
	{
		m_boundaryBoxMin.set(xMin, yMin, zMin);
		m_boundaryBoxMax.set(xMax, yMax, zMax);
		m_boundaryBoxEmpty = false;
	}
	else
	{
		m_boundaryBoxMin.zero();
		m_boundaryBoxMax.zero();
		m_boundaryBoxEmpty = true;
	}
}


//---------------------------------------------------------------
bool AllisNum(std::string str)
{
	for (int i = 0; i < str.size(); i++)
	{
		int tmp = (int)str[i];
		if (tmp >= 48 && tmp <= 57)
		{
			continue;
		}
	    else
		{
		    return false;
		}
	}

    return true;
}

bool ReadASC(const char *filename, cVertexArrayPtr _vertices, int &pointCount)
{
	std::ifstream dataFile;

	dataFile.open(filename);

	if (!dataFile)
	{
		std::cout << "文件路径有误，读取失败，按任意键退出。" << std::endl;
		getchar();
		return false;
	}

	std::string str;
	_vertices->clear();
	unsigned index = 0;

	while (dataFile.peek() != EOF)
	{
		dataFile >> str;
		if (str == "particles")
		{
			dataFile >> pointCount;
		}
		else if (AllisNum(str) && index < pointCount)
		{
			double posX, posY, posZ;
			dataFile >> posX >> posY >> posZ;
			//pointPos[index++] = cVector3d(posX, posY, posZ);

			unsigned index = _vertices->newVertex();

			_vertices->setLocalPos(index, cVector3d(posX, posY, posZ));
			_vertices->setNormal(index, cVector3d(1, 0, 0));
			_vertices->setTexCoord(index, cVector3d(0, 0, 0));
			_vertices->setColor(index, cColorf(1.0, 1.0, 1.0));
		}
		std::getline(dataFile, str);
	}

	dataFile.close();

	return true;
}

//------------------------------------------------------------------------------
// CUDA releated:
//------------------------------------------------------------------------------
void PhysModel::FindExtremePoint()
{
	extremePoint_xMinus = cVector3d(FLT_MAX, 0, 0);
	extremePoint_xPlus = cVector3d(-FLT_MAX, 0, 0);
	extremePoint_yMinus = cVector3d(0, FLT_MAX, 0);
	extremePoint_yPlus = cVector3d(0, -FLT_MAX, 0);
	extremePoint_zMinus = cVector3d(0, 0, FLT_MAX);
	extremePoint_zPlus = cVector3d(0, 0, -FLT_MAX);

	int pointNum = vertices->m_localPos.size();
	cVector3d *data = vertices->m_localPos.data();

	for (int i = 0; i < pointNum; i++)
	{
		if (data[i].x() < extremePoint_xMinus.x())
			extremePoint_xMinus = data[i];
		else if (data[i].x() > extremePoint_xPlus.x())
			extremePoint_xPlus = data[i];

		if (data[i].y() < extremePoint_yMinus.y())
			extremePoint_yMinus = data[i];
		else if (data[i].y() > extremePoint_yPlus.y())
			extremePoint_yPlus = data[i];

		if (data[i].z() < extremePoint_zMinus.z())
			extremePoint_zMinus = data[i];
		else if (data[i].z() > extremePoint_zPlus.z())
			extremePoint_zPlus = data[i];
	}
}

double PhysModel::FindNearstPoint(int srcIndex, int& dstIndex)
{
	cVector3d srcPoint = vertices->m_localPos[srcIndex];
	int pointNum = vertices->m_localPos.size();
	cVector3d *data = vertices->m_localPos.data();

	double miniDist = FLT_MAX;
	int targIndex = 0;

	for (int i = 0; i < pointNum; i++)
	{
		if (i == srcIndex)
			continue;

		if (cDistanceSq(srcPoint, data[i]) < miniDist)
		{
			targIndex = i;
			miniDist = cDistanceSq(srcPoint, data[i]);
		}
	}

	dstIndex = targIndex;
	return sqrt(miniDist);
}

int PhysModel::FindPointCountWithinDst(int srcIndex, double dst)
{
	cVector3d srcPoint = vertices->m_localPos[srcIndex];
	int pointNum = vertices->m_localPos.size();
	cVector3d *data = vertices->m_localPos.data();

	int count = 0;
	for (int i = 0; i < pointNum; i++)
	{
		if (i == srcIndex)
			continue;

		if (cDistanceSq(srcPoint, data[i]) - dst * dst <= FLT_MIN)
			count++;
	}

	return count;
}

//debug
void PhysModel::OutputInfo(ArrayType type, uint start, uint count)
{
	switch (type)
	{
	case POSITION:
	{
		CopyArrayFromDevice(hPos.data(), 0, &vertsPos_resource, sizeof(double3) * numParticles);
		for (unsigned i = start; i < start + count; i++)
		{
			std::cout << "pos" << i << ": " << hPos[i].x << " " << hPos[i].y << " " << hPos[i].z << std::endl;
		}
		break;
	}
	case VELOCITY:
	{
		CopyArrayFromDevice(hVel.data(), dVel, 0, sizeof(double3) * numParticles);
		for (unsigned i = start; i < start + count; i++)
		{
			std::cout << "vel" << i << ": " << hVel[i].x << " " << hVel[i].y << " " << hVel[i].z << std::endl;
		}
		break;
	}
	case NEIGHBOR:
	{
		CopyArrayFromDevice(hNeighbors.data(), dNeighbors, 0, sizeof(int) * 26 * numParticles);
		CopyArrayFromDevice(hNeighborCount.data(), dNeighborCount, 0, sizeof(uint) * numParticles);
		for (unsigned i = start; i < start + count; i++)
		{
			uint neighborCount = hNeighborCount[i];
			std::cout << "particle" << i << " has " << neighborCount << " neighbors: ";

			for (unsigned j = 0; j < neighborCount; j++)
			{
				std::cout << hNeighbors[i * 26 + j] << " ";
			}
			std::cout << std::endl;
		}
		break;
	}
	case CELLSTART:
	{
		CopyArrayFromDevice(hCellStart.data(), dCellStart, 0, sizeof(uint) * physParam.numCells);
		for (unsigned i = start; i < start + count; i++)
		{
			std::cout << "cell" << i << " start index: " << hCellStart[i] << std::endl;
		}
		break;
	}
	case CELLEND:
	{
		CopyArrayFromDevice(hCellEnd.data(), dCellEnd, 0, sizeof(uint) * physParam.numCells);
		for (unsigned i = start; i < start + count; i++)
		{
			std::cout << "cell" << i << " end index: " << hCellEnd[i] << std::endl;
		}
		break;
	}
	default:
		break;
	}
}

void PhysModel::UpdateNeighbors()
{
	double3* dPos = (double3*)MapGLBufferObject(&vertsPos_resource);

	// calculate grid hash
	CalcHash(dGridParticleHash, dGridParticleIndex, dPos, numParticles);

	// sort particles index based on hash
	SortParticles(dGridParticleHash, dGridParticleIndex, numParticles);

	// find start and end of each cell
	FindCellStartEnd(dCellStart, dCellEnd, dGridParticleHash, numParticles, physParam.numCells);

	// reorder particle arrays into sorted order
	ReorderData(dSortedPos, dSortedVel, dGridParticleIndex, dPos, dVel, numParticles);

	CopyArrayFromDevice(hSortedPos.data(), dSortedPos, 0, sizeof(double3) * numParticles);

	//find neighbors
	FindNeighborsWithinDst(dNeighbors, dNeighborCount, dGridParticleIndex, dSortedPos, 0.02, dCellStart, dCellEnd, numParticles);

	UnMapGLBufferObject(vertsPos_resource);

	CopyArrayFromDevice(hNeighbors.data(), dNeighbors, 0, sizeof(int) * 26 * numParticles);
	CopyArrayFromDevice(hNeighborCount.data(), dNeighborCount, 0, sizeof(uint) * numParticles);
}


//------------------------------------------------------------------------
// host find all particles neighbors
//------------------------------------------------------------------------

int3 PhysModel::HostCalcGridPos(double3 pos)
{
	int3 gridPos;
	gridPos.x = floor((pos.x - physParam.worldOrigin.x) / physParam.cellSize.x);
	gridPos.y = floor((pos.y - physParam.worldOrigin.y) / physParam.cellSize.y);
	gridPos.z = floor((pos.z - physParam.worldOrigin.z) / physParam.cellSize.z);

	return gridPos;
}

uint PhysModel::HostCalcGridHash(int3 gridPos)
{
	gridPos.x = gridPos.x & (physParam.gridSize.x - 1);  // wrap grid, assumes size is power of 2
	gridPos.y = gridPos.y & (physParam.gridSize.y - 1);
	gridPos.z = gridPos.z & (physParam.gridSize.z - 1);

	return gridPos.z * physParam.gridSize.y * physParam.gridSize.x + gridPos.y * physParam.gridSize.x + gridPos.x;
}

void PhysModel::HostCalcHash()
{
	glBindBuffer(GL_ARRAY_BUFFER, vertices->m_positionBuffer);
	double3 *hPos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
	for (unsigned i = 0; i < numParticles; i++)
	{
		int3 gridPos = HostCalcGridPos(hPos[i]);
		uint hash = HostCalcGridHash(gridPos);

		hGridParticleHash[i] = hash;
		hGridParticleIndex[i] = i;
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
}

void PhysModel::HostSortParticles()
{
	/*std::vector<std::pair<uint, uint>> pairs;
	pairs.resize(numParticles);
	for (unsigned i = 0; i < numParticles; i++)
	{
		pairs[i] = std::make_pair(hGridParticleHash[i], hGridParticleIndex[i]);
	}

	std::sort(pairs.data(), pairs.data() + numParticles, [](const std::pair<uint, uint>&firstElem, const std::pair<uint, uint>& secElem) {
		return firstElem.first < secElem.first;
	});

	for (unsigned i = 0; i < numParticles; i++)
	{
		hGridParticleHash[i] = pairs[i].first;
		hGridParticleIndex[i] = pairs[i].second;
	}*/
	thrust::sort_by_key(thrust::host, hGridParticleHash.data(), hGridParticleHash.data() + numParticles, hGridParticleIndex.data());
}

void PhysModel::HostFindCellStartEnd()
{
	for (unsigned i = 0; i < numParticles; i++)
	{
		uint hash = hGridParticleHash[i];

		if (i == 0 || (i > 0 && hash != hGridParticleHash[i - 1]))
		{
			hCellStart[hash] = i;
			if (i > 0)
				hCellEnd[hGridParticleHash[i - 1]] = i;
		}
		if (i == numParticles - 1)
		{
			hCellEnd[hash] = i + 1;
		}
	}
}

void PhysModel::HostReorderData()
{
	glBindBuffer(GL_ARRAY_BUFFER, vertices->m_positionBuffer);
	double3 *hPos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

	for (unsigned i = 0; i < numParticles; i++)
	{
		uint sortedIndex = hGridParticleIndex[i];

		hSortedPos[i] = hPos[sortedIndex];
		hSortedVel[i] = hVel[sortedIndex];
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
}

void PhysModel::HostFindNeighborsWithinDst(double dst)
{
	for (unsigned i = 0; i < numParticles; i++)
	{
		uint originalIndex = hGridParticleIndex[i];
		double3 pos = hSortedPos[i];
		int3 gridPos = HostCalcGridPos(pos);

		uint neighborId = 0;
		for (int z = -1; z <= 1; z++)
		{
			for (int y = -1; y <= 1; y++)
			{
				for (int x = -1; x <= 1; x++)
				{
					int3 neighbourGridPos = gridPos + make_int3(x, y, z);
					uint gridHash = HostCalcGridHash(neighbourGridPos);

					uint startIndex = hCellStart[gridHash];
					if (startIndex != 0xffffffff) //cell is not empty
					{
						uint endIndex = hCellEnd[gridHash];
						for (int j = startIndex; j < endIndex; j++)
						{
							if (j != i)
							{
								double3 targPos = hSortedPos[j];
								if (YH::Distance(pos, targPos) <= dst * dst)
								{
									if (neighborId >= 26)
										break;

									uint originalTargIndex = hGridParticleIndex[j];

									hNeighbors[originalIndex * 26 + neighborId] = originalTargIndex; //默认每个particle周围最多有26个邻居
									neighborId++;
								}
							}
						}
						hNeighborCount[originalIndex] = neighborId;
					}
				}
			}
		}
	}
}

void PhysModel::HostUpdateNeighbors()
{
	HostCalcHash();
	HostSortParticles();
	HostFindCellStartEnd();
	HostReorderData();
	HostFindNeighborsWithinDst(physParam.Solid_1_EffectiveRadius * 5); //*5目的是找到26个邻居
}