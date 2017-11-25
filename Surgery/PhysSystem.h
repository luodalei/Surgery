#pragma once

#include "math/CVector3d.h"
#include "world/CMultiPoint.h"
#include "world/CMultiSegment.h"
#include "world/CGenericObject.h"
#include <vector>
#include <iostream>

//CUDA
#include "Deformation.cuh"

using namespace chai3d;

class PhysSystem : public cGenericObject
{
public:
	PhysSystem(std::string filename);
	virtual ~PhysSystem();

	bool LoadFromASC(std::string filename);

	


public:
	cMultiPoint *pointCloud;
	cMultiSegment *stretches;
	cVertexArrayPtr vertices;

	bool isShowPointCloud;
};

