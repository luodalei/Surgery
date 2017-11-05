#include "PhysModel.h"
#include "collisions/CGenericCollision.h"
#include "collisions/CCollisionBrute.h"
#include "collisions/CCollisionAABB.h"
#include "shaders/CShaderProgram.h"

#include <iostream>
#include <fstream>
#include <string>

#include "GL/glu.h"

//四种状态粒子初始质量
double PhysConstant::SOLID_1_Mass;
double PhysConstant::SOLID_2_Mass;
double PhysConstant::LIQUID_Mass;
double PhysConstant::GAS_Mass;

double PhysConstant::Solid_Density;
double PhysConstant::Liquid_Density;
double PhysConstant::Gas_Density;

double PhysConstant::Air_Temperature;

double PhysConstant::Split_Temperature; //SOLID_1 -> SOLID_2
double PhysConstant::Fusion_Temperature;//SOLID_2 -> LIQUID
double PhysConstant::Boil_Temperature; //LIQUID -> GAS
double PhysConstant::LatentHeat_Solid2Liquid;//SOLID_2 -> LIQUID
double PhysConstant::LatentHeat_Liquid2Gas;

double PhysConstant::Thermal_Conductivity; //h: Qi=h(Tair - Ti)*A
double PhysConstant::Thermal_Diffusion; //Cd： in heat transfer between particle
double PhysConstant::Heat_Capacity_Solid; //C: delta_T=Qi/C*m
double PhysConstant::Heat_Capacity_Liquid;
double PhysConstant::Heat_Capacity_Gas;

//光子参数
double PhysConstant::Emissivity;
double PhysConstant::Boltzmann_Constant;
double PhysConstant::Source_Temperature;
double PhysConstant::Source_Area;
int PhysConstant::PhotonsNum_TimeStep;//单位时间射出的光子数

//各种半径参数
double PhysConstant::Solid_1_EffectiveRadius; //re
double PhysConstant::Solid_2_EffectiveRadius;


PhysModel::PhysModel()
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
}


PhysModel::~PhysModel()
{
	m_displayList.invalidate();
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
			OutputVertsPos();

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