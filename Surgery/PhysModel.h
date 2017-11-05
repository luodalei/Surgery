#pragma once
#include "math/CVector3d.h"
#include "graphics/CPointArray.h"
#include "world/CGenericObject.h"
#include <vector>
#include <iostream>

//CUDA
#include "Deformation.cuh"

using namespace chai3d;

enum ParticleType { SOLID_1, SOLID_2, LIQUID, GAS };

typedef struct _PhysConstant
{
	//四种状态粒子初始质量
	static double SOLID_1_Mass;
	static double SOLID_2_Mass;
	static double LIQUID_Mass;
	static double GAS_Mass;

	static double Solid_Density;
	static double Liquid_Density;
	static double Gas_Density;

	static double Air_Temperature;

	static double Split_Temperature; //SOLID_1 -> SOLID_2
	static double Fusion_Temperature;//SOLID_2 -> LIQUID
	static double Boil_Temperature; //LIQUID -> GAS
	static double LatentHeat_Solid2Liquid;//SOLID_2 -> LIQUID
	static double LatentHeat_Liquid2Gas;

	static double Thermal_Conductivity; //h: Qi=h(Tair - Ti)*A
	static double Thermal_Diffusion; //Cd： in heat transfer between particle
	static double Heat_Capacity_Solid; //C: delta_T=Qi/C*m
	static double Heat_Capacity_Liquid;
	static double Heat_Capacity_Gas;

	//光子参数
	static double Emissivity;
	static double Boltzmann_Constant;
	static double Source_Temperature;
	static double Source_Area;
	static int PhotonsNum_TimeStep;//单位时间射出的光子数

	//各种半径参数
	static double Solid_1_EffectiveRadius; //re
	static double Solid_2_EffectiveRadius;


} PhysConstant;


class PhysModel : public cGenericObject
{
public:
	PhysModel();
	virtual ~PhysModel();
	bool LoadFromASC(std::string filename);
	void MarkForUpdate(const bool a_affectChildren);

	//--------------------------------------------------------------------------
	// PUBLIC METHODS - VERTICES
	//--------------------------------------------------------------------------
	//! This method creates a new vertex and adds it to the vertex list.
	unsigned NewVertex(const double a_x = 0.0, const double a_y = 0.0, const double a_z = 0.0);

	//! This method creates a new vertex and adds it to the vertex list.
	unsigned NewVertex(const cVector3d& a_pos, const cColorf& a_color = cColorf(0, 0, 0, 1));

	//! This method returns the number of stored vertices.
	inline unsigned GetNumVertices() const { return (unsigned int)(vertices->getNumElements()); }

	//--------------------------------------------------------------------------
	// PUBLIC METHODS - POINTS
	//--------------------------------------------------------------------------
	//! This method creates a new point by passing a vertex index.
	unsigned NewPoint(const unsigned a_indexVertex0);

	//! This method creates a new point by passing a vertex position and color.
	unsigned NewPoint(const cVector3d& a_vertex0 = cVector3d(0, 0, 0), const cColorf& a_colorVertex0 = cColorf(0, 0, 0, 1));

	//! This method removed a selected point.
	bool RemovePoint(const unsigned int a_index);

	//! This method returns the number of stored points.
	unsigned GetNumPoints();

	//! This method clears all points and vertices.
	void Clear();

	//-----------------------------------------------------------------------
	// PUBLIC METHODS - GRAPHIC PROPERTIES:
	//-----------------------------------------------------------------------

	//! This method sets the alpha value at each vertex and in all of my material colors.
	virtual void setTransparencyLevel(const float a_level, const bool a_applyToVertices = false, const bool a_applyToTextures = false, const bool a_affectChildren = false);


	//--------------------------------------------------------------------------
	// PUBLIC METHODS - POINTS GRAPHIC PROPERTIES
	//--------------------------------------------------------------------------

	//! This method sets a color for all points.
	void SetPointColor(const cColorf& a_color);

	//! This method sets the point size that is used to graphically render the points.
	inline void SetPointSize(const double a_pointSize) { m_pointSize = fabs(a_pointSize); }

	//! This method returns the point size.
	inline double GetPointSize() const { return (m_pointSize); }

	//--------------------------------------------------------------------------
	// PUBLIC METHODS - COLLISION DETECTION:
	//--------------------------------------------------------------------------

	//! This method builds a brute force collision detector for this mesh.
	virtual void createBruteForceCollisionDetector();

	//! This method builds an AABB collision detector for this mesh.
	virtual void createAABBCollisionDetector(const double a_radius);

	//--------------------------------------------------------------------------
	// PUBLIC METHODS - GEOMETRY:
	//--------------------------------------------------------------------------

	//! This method scales this point cloud by using different scale factors along X, Y and Z axes.
	void ScaleXYZ(const double a_scaleX, const double a_scaleY, const double a_scaleZ);

	//! This method shifts all vertex positions by the specified amount.
	virtual void offsetVertices(const cVector3d& a_offset, const bool a_updateCollisionDetector = true);

	//! This method computes the center of mass of this point cloud, based on vertex positions.
	virtual cVector3d getCenterOfMass();

	void DrawAsSphere(bool _flag) { isDrawAsSphere = _flag; }


	//--------------------------------------------------------------------------
	// PROTECTED METHODS:
	//--------------------------------------------------------------------------

protected:

	//! This method renders this object graphically using OpenGL.
	virtual void render(cRenderOptions& a_options);

	//! This method renders all points.
	virtual void renderPoints(cRenderOptions& a_options);

	//这个方法把点渲染成球
	void GenerateSphereVertices(std::vector<cVector3d>& _vertices, std::vector<int>& _indices, int _lats, int _longs);
	void RenderSpheres(cRenderOptions& a_options);

	//! This method updates the global position of each vertex.
	virtual void updateGlobalPositions(const bool a_frameOnly);

	//! This method updates the boundary box dimensions based on the vertices.
	virtual void updateBoundaryBox();

	//! This method scales this object by a scale factor.
	virtual void scaleObject(const double& a_scaleFactor) { ScaleXYZ(a_scaleFactor, a_scaleFactor, a_scaleFactor); }

	//------------------------------------------------------------------------------
	// CUDA releated:
	//------------------------------------------------------------------------------
public:
	bool isRegisterBuffer;
	
public:
	cPointArrayPtr points;
	cVertexArrayPtr vertices;	
	bool isDrawAsSphere;
	double sphereRadius;

protected:
	bool m_showPoints; //! If __true__, then segments are displayed.
	double m_pointSize; //! Display size of point.

public:
	std::vector<double> masses;
	std::vector<cVector3d> velocities;
	std::vector<double> temperatures;
	std::vector<double> lastTemperatures;
	std::vector<ParticleType> particleTypes;
	std::vector<double> heatQ;

};


bool ReadASC(const char *filename, cVertexArrayPtr _vertices, int &pointCount);