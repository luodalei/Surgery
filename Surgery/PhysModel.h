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




class PhysModel : public cGenericObject
{
public:
	PhysModel(std::string filename);
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
	enum ArrayType
	{
		POSITION,
		VELOCITY,
		NEIGHBOR,
		CELLSTART,
		CELLEND
	};

	//find xyz-minus xyz-plus extreme point
	void FindExtremePoint();
	double FindNearstPoint(int srcIndex, int& dstIndex);
	int FindPointCountWithinDst(int srcIndex, double dst);

	void OutputInfo(ArrayType type, uint start, uint count);//debug
	void FreeMemory(); //释放与物理属性、cuda相关内存
	void UpdateNeighbors();

	int3 HostCalcGridPos(double3 p);
	uint HostCalcGridHash(int3 gridPos);
	void HostCalcHash();
	void HostSortParticles();
	void HostFindCellStartEnd();
	void HostReorderData();
	void HostFindNeighborsWithinDst(double dst);
	void HostUpdateNeighbors();

public:
	bool isRegisterBuffer;

	cVector3d extremePoint_xMinus;
	cVector3d extremePoint_xPlus;
	cVector3d extremePoint_yMinus;
	cVector3d extremePoint_yPlus;
	cVector3d extremePoint_zMinus;
	cVector3d extremePoint_zPlus;
	
public:
	cPointArrayPtr points;
	cVertexArrayPtr vertices;

	bool isDrawAsSphere;
	double sphereRadius;
	PhysParam physParam;
	uint numParticles;
	bool needUpdateNeighbors;

	//CPU data
	struct
	{
		std::vector<double3> hPos; //store vertices->m_localPos, 用作输出，注意m_localpos中的值会被缩放影响，hPos可能需要相应更新
		std::vector<double3> hVel; //store velocity
		struct cudaGraphicsResource *vertsPos_resource;

		//data for neighbor find
		std::vector<int> hNeighbors; //保存邻居索引，length = 26 * particles count
		std::vector<uint> hNeighborCount;

		std::vector<double3> hSortedPos;
		std::vector<double3> hSortedVel;

		std::vector<uint> hGridParticleHash;//暂时没用
		std::vector<uint> hGridParticleIndex;
		std::vector<uint> hCellStart;
		std::vector<uint> hCellEnd;

		//data for haptic force
		std::vector<int> hNearbyIdx;
		uint hNearbyCount;

		//data for pbd
		std::vector<double3> hPredictPos;
		std::vector<double> hInvMass;
	};

	//GPU data
	struct
	{
		double3 *dVel; //velocity on GPU

		//data for neighbor find
		int* dNeighbors; //length = 26 * particles count
		uint* dNeighborCount;

		double3 *dSortedPos;
		double3 *dSortedVel;

		uint *dGridParticleHash; // grid hash value for each particle
		uint *dGridParticleIndex;// particle index for each particle
		uint *dCellStart;        // index of start of each cell in sorted list
		uint *dCellEnd;          // index of end of cell

		//data for haptic force
		int *dNearbyIdx;
		uint *dNearbyCount;

		//data for pbd
		double3 *dPredictPos;
		double *dInvMass;
	};
protected:
	bool m_showPoints; //! If __true__, then segments are displayed.
	double m_pointSize; //! Display size of point.

public:
	
	std::vector<double> temperatures;
	std::vector<double> lastTemperatures;
	std::vector<ParticleType> particleTypes;
	std::vector<double> heatQ;

};


bool ReadASC(const char *filename, cVertexArrayPtr _vertices, int &pointCount);