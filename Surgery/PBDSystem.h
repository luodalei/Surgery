#pragma once

#include "PhysModel.h"

using namespace chai3d;

class PBDSystem
{
public:
	PBDSystem(PhysModel* _physModel);
	~PBDSystem();
	void FreeMemory();

	void ComputeDistanceConstraints();
	void ComputeVolumeConstraints();

	void DoPBD(cVector3d toolPosInModelSpace, cVector3d toolForceInModelSpace);
	void HostDoPBD(cVector3d toolPosInModelSpace, cVector3d toolForceInModelSpace, bool applyForceOnPhysModel);

	void HostFindParticlesWithinDstFromTarget(double3 targPos, double dst);

	void HostUpdateVelocities(double dt);
	void HostIntegrateExplicitWithDamping(double dt);
	bool HostUpdateDistanceConstraints(unsigned idx);
	bool HostUpdateVolumeConstraints(unsigned idx);

	bool HostInitShapeMatchingConstraints();
	bool HostUpdateShapeMatchingConstraints(bool allowStretch, Mat3d *rot = nullptr);

	void HostInitCubeVolumeConstraints();
	bool HostUpdateCubeVolumeConstraints(unsigned idx);

	void HostUpdateAllConstraints();
	void HostIntegrate(double dt);

public:

	struct int7
	{
		int a, b, c, d, e, f, g;
	};
	
	//CPU data
	struct
	{
		std::vector<int2> hDistanceConstraints;
		std::vector<double> hRestLength;
		std::vector<double3> hDistanceCorrSum;
		uint hNumDistanceConstraints;

		std::vector<int4> hVolumeConstraints;
		std::vector<double> hRestVolume;
		uint hNumVolumeConstraints;

		std::vector<int7> hCubeVolumeConstraints;
		std::vector<double> hRestCubeVolume;
		std::vector<double3> hCubeVolumeCorrSum;
		uint hNumCubeVolumeConstraints;

		std::vector<unsigned> hNumClusters;
	};
	
	//GPU data
	struct
	{
		int2 *dDistanceConstraints;
		int4 *dVolumeConstraints;
		double *dRestLength;
		uint *dNumDistanceConstraint;
	};

	PhysModel* physModel;

	double3 toolPos;
	double3 toolForce;
	double3 environmentForce;

	//distance constraints
	bool needRebuildDistanceConstraintArray;
	double compressionStiffness;
	double stretchStiffness;
	double fkDamping;

	//volume constraints
	bool needRebuildVolumeConstraintArray;
	double posVolumeStiffness;
	double negVolumeStiffness;

	//shape matching constraints
	bool needRecomputeShapeMatching;
	double3 restCm;
	Mat3d invRestMat;
	double shapeMatchingStiffness;

	//cube volume constraints
	bool needRebuildCubeVolumeArray;
	double posCubeVolumeStiffness;
	double negCubeVolumeStiffness;

	double dt;
};

