#include "PBDSystem.h"

const double eps = 1e-6;

PBDSystem::PBDSystem(PhysModel* _physModel)
{
	physModel = _physModel;

	//distance constraint
	AllocateArray((void**)&dDistanceConstraints, sizeof(int2) * physModel->numParticles * 26 / 2); //不重复包含索引

	//Allocate memsize for dRestLength
	AllocateArray((void**)&dRestLength, sizeof(double) * physModel->numParticles * 26 / 2);

	//Allocate memsize for numconstraint
	AllocateArray((void**)&dNumDistanceConstraint, sizeof(uint));

	toolPos.x = toolPos.y = toolPos.z = 0.0;
	toolForce.x = toolForce.y = toolForce.z = 0.0;
	environmentForce.x = environmentForce.y = environmentForce.z = 0.0;

	//distance constraint const
	needRebuildDistanceConstraintArray = true;
	compressionStiffness = /*0.8*/1.0;
	stretchStiffness = /*0.8*/1.0;
	fkDamping = 0.0;
	hDistanceCorrSum.resize(physModel->numParticles);

	//volume constraint const
	needRebuildVolumeConstraintArray = true;
	posVolumeStiffness = 0.8;
	negVolumeStiffness = 0.8;

	//shape matching constraint const
	needRecomputeShapeMatching = true;
	invRestMat.SetZero();
	shapeMatchingStiffness = 0.8;

	//cube volume constraint const
	needRebuildCubeVolumeArray = true;
	posCubeVolumeStiffness = 1.0;
	negCubeVolumeStiffness = 1.0;
	hCubeVolumeCorrSum.resize(physModel->numParticles);

	dt = 0.01;
}


PBDSystem::~PBDSystem()
{
	FreeMemory();
}

void PBDSystem::FreeMemory()
{
	hDistanceConstraints.clear();
	hRestLength.clear();

	FreeArray(dDistanceConstraints);
	FreeArray(dRestLength);
	FreeArray(dNumDistanceConstraint);
}

//--------------------------------------------------------------------------
//Compute constraints
//--------------------------------------------------------------------------
void PBDSystem::ComputeDistanceConstraints()
{
	if (needRebuildDistanceConstraintArray)
	{
		glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
		double3 *hPos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

		int *neighbors = physModel->hNeighbors.data();
		uint *neighborCount = physModel->hNeighborCount.data();
		uint numParticles = physModel->numParticles;

		for (unsigned h = 0; h < numParticles; h++)
		{
			for (unsigned j = 0; j < neighborCount[h]; j++)
			{
				int neighborIndex = neighbors[h * 26 + j];
				//std::cout << h << "-" << j << ": " << neighborIndex<< std::endl;
				if (neighborIndex > h)
				{
					//vertex pair index
					int2 pair = make_int2(h, neighborIndex);
					hDistanceConstraints.push_back(pair);

					//distance of constraint
					double3 dis;
					dis = hPos[h] - hPos[neighborIndex];
					hRestLength.push_back(YH::Length(dis));
				}
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		hNumDistanceConstraints = hDistanceConstraints.size();

		CopyArrayToDevice(dDistanceConstraints, hDistanceConstraints.data(), 0, sizeof(int2) * hDistanceConstraints.size());
		CopyArrayToDevice(dRestLength, hRestLength.data(), 0, sizeof(double) * hRestLength.size());
		
		/*double3 *dPos = (double3*)MapGLBufferObject(&physModel->vertsPos_resource);
		UpdateDistanceConstraints(dDistanceConstraints, dRestLength, dNumDistanceConstraint, dPos, physModel->dNeighbors, physModel->dNeighborCount, physModel->numParticles);

		CopyArrayFromDevice(&hNumDistanceConstraints, dNumDistanceConstraint, 0, sizeof(uint));

		UnMapGLBufferObject(physModel->vertsPos_resource);*/

		needRebuildDistanceConstraintArray = false;
	}
}

void PBDSystem::ComputeVolumeConstraints()
{
	if (needRebuildVolumeConstraintArray)
	{
		glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
		double3 *hPos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

		int *neighbors = physModel->hNeighbors.data();
		uint *neighborCount = physModel->hNeighborCount.data();
		uint numParticles = physModel->numParticles;
		std::vector<uint> firstClass;
		std::vector<uint> secondClass;

		for (unsigned i = 0; i < numParticles; i++)
		{
			if (neighborCount[i] == 6)
			{
				for (unsigned j = 0; j < 6; j++)
				{
					uint neighborIdx = neighbors[i * 26 + j];
					if (hPos[neighborIdx].x - hPos[i].x > 0.01 || hPos[neighborIdx].y - hPos[i].y > 0.01 || hPos[neighborIdx].z - hPos[i].z < -0.01)
					{
						firstClass.push_back(neighborIdx);
					}
					else if (hPos[neighborIdx].x - hPos[i].x < -0.01 || hPos[neighborIdx].y - hPos[i].y < -0.01 || hPos[neighborIdx].z - hPos[i].z > 0.01)
					{
						secondClass.push_back(neighborIdx);
					}
				}
				if (firstClass.size() == 3)
				{
					int4 tet = make_int4(i, firstClass[0], firstClass[1], firstClass[2]);
					hVolumeConstraints.push_back(tet);

					//compute rest volume
					double3 p0 = hPos[i];
					double3 p1 = hPos[firstClass[0]];
					double3 p2 = hPos[firstClass[1]];
					double3 p3 = hPos[firstClass[2]];

					double restVolume = fabs((1.0 / 6.0) * YH::Dot(p3 - p0, YH::Cross(p2 - p0, p1 - p0)));
					hRestVolume.push_back(restVolume);
				}
				if (secondClass.size() == 3)
				{
					int4 tet = make_int4(i, secondClass[0], secondClass[1], secondClass[2]);
					hVolumeConstraints.push_back(tet);

					//compute rest volume
					double3 p0 = hPos[i];
					double3 p1 = hPos[secondClass[0]];
					double3 p2 = hPos[secondClass[1]];
					double3 p3 = hPos[secondClass[2]];

					double restVolume = fabs((1.0 / 6.0) * YH::Dot(p3 - p0, YH::Cross(p2 - p0, p1 - p0)));
					hRestVolume.push_back(restVolume);
				}
				firstClass.clear();
				secondClass.clear();
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		hNumVolumeConstraints = hVolumeConstraints.size();
		needRebuildVolumeConstraintArray = false;
	}
}

bool PBDSystem::HostInitShapeMatchingConstraints()
{
	if (needRecomputeShapeMatching)
	{
		//physModel->vertices->m_localPos中的数据可以作为初始不变数据
		cVector3d *pos0 = physModel->vertices->m_localPos.data();

		invRestMat.Identity();

		//center of mass
		double wsum = 0.0;
		cVector3d restCm3d; restCm3d.zero();
		for (unsigned i = 0; i < physModel->numParticles; i++)
		{
			double wi = 1.0 / (physModel->hInvMass[i] + eps);
			//std::cout << "wi: " << wi << std::endl;
			//std::cout << "pos0: " << pos0[i].x() << " " << pos0[i].y() << " " << pos0[i].z() << std::endl;
			restCm3d += pos0[i] * wi;
			wsum += wi;
		}
		if (wsum == 0.0)
			return false;
		restCm.x = restCm3d.x(); restCm.y = restCm3d.y(); restCm.z = restCm3d.z();

		//A
		Mat3d A; A.SetZero();
		for (unsigned i = 0; i < physModel->numParticles; i++)
		{
			const cVector3d qi = pos0[i] - restCm3d;

			double wi = 1.0 / (physModel->hInvMass[i] + eps);
			double x2 = wi * qi.x() * qi.x();
			double y2 = wi * qi.y() * qi.y();
			double z2 = wi * qi.z() * qi.z();
			double xy = wi * qi.x() * qi.y();
			double xz = wi * qi.x() * qi.z();
			double yz = wi * qi.y() * qi.z();

			A(0, 0) += x2; A(0, 1) += xy; A(0, 2) += xz;
			A(1, 0) += xy; A(1, 1) += y2; A(1, 2) += yz;
			A(2, 0) += xz; A(2, 1) += yz; A(2, 2) += z2;
		}

		double det = A.Det();
		if (fabs(det) > eps)
		{
			A.Invertr(invRestMat);

			needRecomputeShapeMatching = false;
			return true;
		}
		return false;
	}
}

void PBDSystem::HostInitCubeVolumeConstraints()
{
	if (needRebuildCubeVolumeArray)
	{
		glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
		double3 *hPos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

		int *neighbors = physModel->hNeighbors.data();
		uint *neighborCount = physModel->hNeighborCount.data();
		uint numParticles = physModel->numParticles;

		for (unsigned i = 0; i < numParticles; i++)
		{
			if (neighborCount[i] == 6)
			{
				int7 cubeIndex;
				for (unsigned j = 0; j < 6; j++)
				{
					uint neighborIdx = neighbors[i * 26 + j];

					cubeIndex.a = i; //center
					if (hPos[neighborIdx].x - hPos[i].x > 0.01) //+x
						cubeIndex.b = neighborIdx;
					else if (hPos[neighborIdx].x - hPos[i].x < -0.01) //-x
						cubeIndex.c = neighborIdx;
					else if (hPos[neighborIdx].y - hPos[i].y > 0.01) //+y
						cubeIndex.d = neighborIdx;
					else if (hPos[neighborIdx].y - hPos[i].y < -0.01)//-y
						cubeIndex.e = neighborIdx;
					else if (hPos[neighborIdx].z - hPos[i].z > 0.01)//+z
						cubeIndex.f = neighborIdx;
					else if (hPos[neighborIdx].z - hPos[i].z < -0.01)//-z
						cubeIndex.g = neighborIdx;
				}
				hCubeVolumeConstraints.push_back(cubeIndex);

				//rest volume
				double cubeVolume;
				double3 p_xplus = hPos[cubeIndex.b];
				double3 p_xminus = hPos[cubeIndex.c];
				double3 p_yplus = hPos[cubeIndex.d];
				double3 p_yminus = hPos[cubeIndex.e];
				double3 p_zplus = hPos[cubeIndex.f];
				double3 p_zminus = hPos[cubeIndex.g];

				cubeVolume = YH::Dot(p_xplus, YH::Cross(p_yplus, p_zplus)) + YH::Dot(p_xplus, YH::Cross(p_yminus, p_zminus)) +
					YH::Dot(p_xplus, YH::Cross(p_yplus, p_zminus)) + YH::Dot(p_xminus, YH::Cross(p_yminus, p_zplus)) +
					YH::Dot(p_zplus, YH::Cross(p_yplus, p_xminus)) + YH::Dot(p_zplus, YH::Cross(p_yminus, p_xplus)) +
					YH::Dot(p_zminus, YH::Cross(p_yplus, p_xplus)) + YH::Dot(p_zminus, YH::Cross(p_yminus, p_xminus));

				hRestCubeVolume.push_back(cubeVolume);
			}
		}
		glUnmapBuffer(GL_ARRAY_BUFFER);

		hNumCubeVolumeConstraints = hCubeVolumeConstraints.size();
		needRebuildCubeVolumeArray = false;
	}
}

//------------------------------------------------------------------------------
//CUDA PBD
//------------------------------------------------------------------------------
void PBDSystem::DoPBD(cVector3d toolPosInModelSpace, cVector3d toolForceInModelSpace)
{
	toolPos.x = toolPosInModelSpace.x(); toolPos.y = toolPosInModelSpace.y(); toolPos.z = toolPosInModelSpace.z();
	toolForce.x = toolForceInModelSpace.x(); toolForce.y = toolForceInModelSpace.y(); toolForce.z = toolForceInModelSpace.z();
	toolForce.x = 0.0; toolForce.y = -5; toolForce.z = 0;

	FindParticlesWithinDstFromTarget(physModel->dNearbyIdx, physModel->dNearbyCount, physModel->dGridParticleIndex, physModel->dSortedPos, toolPos, 0.1, physModel->dCellStart, physModel->dCellEnd);
	
	//测试是否找对了tool周围粒子！成功
	/*CopyArrayFromDevice(physModel->hNearbyIdx.data(), physModel->dNearbyIdx, 0, sizeof(int) * physModel->hNearbyIdx.size());
	CopyArrayFromDevice(&physModel->hNearbyCount, physModel->dNearbyCount, 0, sizeof(uint));

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		physModel->vertices->setColor(i, cColorb(255, 255, 255));
	}

	for (unsigned j = 0; j < physModel->hNearbyCount; j++)
	{
		physModel->vertices->setColor(physModel->hNearbyIdx[j], cColorb(255, 0, 0));
	}*/

	double3 *dPos = (double3*)MapGLBufferObject(&physModel->vertsPos_resource);

	//UpdateVelocities(physModel->dVel, dt, physModel->dNearbyIdx, physModel->dNearbyCount, toolForce, environmentForce, physModel->dInvMass, physModel->numParticles);

	SemiImplicitEuler(physModel->dPredictPos, dPos, physModel->dVel, dt, physModel->dNearbyIdx, physModel->dNearbyCount, toolForce, environmentForce, physModel->dInvMass, physModel->numParticles);

/*	IntegrateExplicitWithDamping(physModel->dPredictPos, //output
								dPos, physModel->dVel, 
								0.9, //kDamping
								dt,
								physModel->numParticles);
*/
	UpdateAllConstraints(physModel->dPredictPos, physModel->dInvMass, dDistanceConstraints, compressionStiffness, stretchStiffness, dRestLength, hNumDistanceConstraints, 5);

	DetermineFinalVel(physModel->dPredictPos, dPos, physModel->dVel, dt, physModel->numParticles);

	CopyArrayFromDeviceToDevice(dPos, physModel->dPredictPos, 0, sizeof(double3) * physModel->numParticles);

	UnMapGLBufferObject(physModel->vertsPos_resource);
}

//---------------------------------------------------------------------------
//host find particles near tool pos
//---------------------------------------------------------------------------

void PBDSystem::HostFindParticlesWithinDstFromTarget(double3 targPos, double dst)
{

	int3 gridPos = physModel->HostCalcGridPos(targPos);
	uint nearbyId = 0;

	//当dst取得比较大时，这里可以放宽搜索范围
	for (int z = -2; z <= 2; z++)
	{
		for (int y = -2; y <= 2; y++)
		{
			for (int x = -2; x <= 2; x++)
			{
				int3 neighbourGridPos = gridPos + make_int3(x, y, z);
				uint gridHash = physModel->HostCalcGridHash(neighbourGridPos);

				uint startIndex = physModel->hCellStart[gridHash];
				if (startIndex != 0xffffffff) //cell is not empty
				{
					uint endIndex = physModel->hCellEnd[gridHash];
					for (int j = startIndex; j < endIndex; j++)
					{
						double3 testPos = physModel->hSortedPos[j];
						if (YH::Distance(targPos, testPos) <= dst * dst)
						{
							uint originalTestPosIndex = physModel->hGridParticleIndex[j];
							physModel->hNearbyIdx[nearbyId] = originalTestPosIndex;
							nearbyId++;
						}
					}
					physModel->hNearbyCount = nearbyId;
				}
			}
		}
	}

}

//---------------------------------------------------------------------------
// host pbd for test
//---------------------------------------------------------------------------

void PBDSystem::HostUpdateVelocities(double dt)
{
	glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
	double3 *pos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		//physModel->hVel[i].x = physModel->hVel[i].y = physModel->hVel[i].z = 0.0;//感觉不应该加这一句

		physModel->hVel[i] += environmentForce * physModel->hInvMass[i] * dt;
		for (unsigned j = 0; j < physModel->hNearbyCount; j++) //在tool附近会受到力
		{
			if (i == physModel->hNearbyIdx[j])
			{
				physModel->hVel[i] += toolForce * physModel->hInvMass[i] * dt;
			}
		}
		physModel->hPredictPos[i] = pos[i] + physModel->hVel[i] * dt;
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
}

void PBDSystem::HostIntegrateExplicitWithDamping(double dt)
{
	double3 averPos, averVel;

	glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
	double3 *pos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		averPos += pos[i];
		averVel += physModel->hVel[i];
	}
	averPos /= physModel->numParticles;
	averVel /= physModel->numParticles;

	if (averVel.x < 0.0000001 && averVel.x > -0.0000001) averVel.x = 0;
	if (averVel.y < 0.0000001 && averVel.y > -0.0000001) averVel.y = 0;
	if (averVel.z < 0.0000001 && averVel.z > -0.0000001) averVel.z = 0;

	cMatrix3d I; I.identity();
	double3 L;
	cVector3d W3d;

	std::vector<cVector3d> R;

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		double3 ri = pos[i] - averPos;
		L += YH::Cross(ri, physModel->hVel[i]);

		cMatrix3d tmp(0, ri.z, -ri.y,
					 -ri.z, 0, ri.x,
					 ri.y, -ri.x, 0);

		cMatrix3d tmpTran;
		tmp.transr(tmpTran);

		I += tmp * tmpTran;

		R.push_back(cVector3d(ri.x, ri.y, ri.z));
	}
	cVector3d L3d(L.x, L.y, L.z);
	W3d = I.invert() * L3d;

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		cVector3d crossVec3d = cCross(W3d, R[i]);
		double3 crossVec; crossVec.x = crossVec3d.x(); crossVec.y = crossVec3d.y(); crossVec.z = crossVec3d.z();
		double3 deltVi = averVel + crossVec - physModel->hVel[i];

		physModel->hVel[i] += fkDamping * deltVi;
		physModel->hPredictPos[i] = pos[i] + physModel->hVel[i] * dt;

		/*if (YH::Length(physModel->hVel[i]) < 0.0000001)
			physModel->hVel[i].x = physModel->hVel[i].y = physModel->hVel[i].z = 0;*/
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
}

bool PBDSystem::HostUpdateDistanceConstraints(unsigned idx)
{
	int2 pair = hDistanceConstraints[idx];
	int p1 = pair.x;
	int p2 = pair.y;

	double w1 = physModel->hInvMass[p1];
	double w2 = physModel->hInvMass[p2];
	double invMass = w1 + w2;

	if (invMass < FLT_MIN)
		return false;

	double3 p12 = physModel->hPredictPos[p1] - physModel->hPredictPos[p2];
	double len = YH::Length(p12);
	double3 p12_norm = YH::Normalize(p12);

	/*if (len < FLT_MIN)
		return false;*/

	double3 corr, corr0, corr1;
	if (len < hRestLength[idx])
		corr = compressionStiffness * p12_norm * (len - hRestLength[idx]) / invMass;
	else
		corr = stretchStiffness * p12_norm * (len - hRestLength[idx]) / invMass;

	//double3 deltP;
	//deltP = (len - hRestLength[idx]) * (p12 / len) * compressionStiffness / invMass;

	corr0 = -w1 * corr;
	corr1 = w2 * corr;

	if (w1 > 0.0)
	{
		physModel->hPredictPos[p1] += corr0;
		//hDistanceCorrSum[p1] += corr0;
	}
	if (w2 > 0.0)
	{
		physModel->hPredictPos[p2] += corr1;
		//hDistanceCorrSum[p2] += corr1;
	}

	return true;
}

bool PBDSystem::HostUpdateVolumeConstraints(unsigned idx)
{
	int4 tet = hVolumeConstraints[idx];
	int id1 = tet.x;
	int id2 = tet.y;
	int id3 = tet.z;
	int id4 = tet.w;

	double3 &vert1 = physModel->hPredictPos[id1];
	double3 &vert2 = physModel->hPredictPos[id2];
	double3 &vert3 = physModel->hPredictPos[id3];
	double3 &vert4 = physModel->hPredictPos[id4];

	double invMass1 = physModel->hInvMass[id1];
	double invMass2 = physModel->hInvMass[id2];
	double invMass3 = physModel->hInvMass[id3];
	double invMass4 = physModel->hInvMass[id4];

	double3 corr1, corr2, corr3, corr4;

	double3 d1 = vert2 - vert1;
	double3 d2 = vert3 - vert2;
	double3 d3 = vert4 - vert3;
	double volume = (1.0 / 6.0) * YH::Dot(vert4 - vert1, YH::Cross(vert3 - vert1, vert2 - vert1));

	if (posVolumeStiffness == 0.0 && volume > 0.0)
		return false;

	if (negVolumeStiffness == 0.0 && volume < 0.0)
		return false;

	double3 grad1 = YH::Cross(vert2 - vert3, vert4 - vert3);
	double3 grad2 = YH::Cross(vert3 - vert1, vert4 - vert1);
	double3 grad3 = YH::Cross(vert1 - vert2, vert4 - vert2);
	double3 grad4 = YH::Cross(vert2 - vert1, vert3 - vert1);

	double lambda = invMass1 * YH::SquareLength(grad1) +
					invMass2 * YH::SquareLength(grad2) +
					invMass3 * YH::SquareLength(grad3) +
					invMass4 * YH::SquareLength(grad4);

	if (fabs(lambda) < eps)
		return false;

	if (volume < 0.0)
		lambda = negVolumeStiffness * (volume - hRestVolume[idx]) / lambda;
	else
		lambda = posVolumeStiffness * (volume - hRestVolume[idx]) / lambda;

	corr1 = -lambda * invMass1 * grad1;
	corr2 = -lambda * invMass2 * grad2;
	corr3 = -lambda * invMass3 * grad3;
	corr4 = -lambda * invMass4 * grad4;

	if (invMass1 > 0.0)
		vert1 += corr1;
	if (invMass2 > 0.0)
		vert2 += corr2;
	if (invMass3 > 0.0)
		vert3 += corr3;
	if (invMass4 > 0.0)
		vert4 += corr4;

	return true;
}



bool PBDSystem::HostUpdateShapeMatchingConstraints(bool allowStretch, Mat3d *rot)
{
	double3 *pos = physModel->hPredictPos.data();
	
	double3 *corr = new double3[physModel->numParticles];

	//center of mass
	double3 cm;
	double wsum = 0.0;
	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		double wi = 1.0 / (physModel->hInvMass[i] + eps);
		cm += pos[i] * wi;
		wsum += wi;
	}
	if (wsum == 0.0)
		return false;
	cm /= wsum;

	//A
	Mat3d mat; mat.SetZero();
	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		double3 q = pos[i] - restCm;
		double3 p = pos[i] - cm;
		double w = 1.0 / (physModel->hInvMass[i] + eps);
		p *= w;

		mat(0, 0) += p.x * q.x; mat(0, 1) += p.x * q.y; mat(0, 2) += p.x * q.z;
		mat(1, 0) += p.y * q.x; mat(1, 1) += p.y * q.y; mat(1, 2) += p.y * q.z;
		mat(2, 0) += p.z * q.x; mat(2, 1) += p.z * q.y; mat(2, 2) += p.z * q.z;
	}
	mat = mat * invRestMat;

	Mat3d R, U, D;
	R = mat;
	if (allowStretch)
		R = mat;
	else
		YH::PolarDecompositionStable(mat, eps, R);

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		double3 goal = cm + R * (pos[i] - restCm);
		std::cout << "goal: " << goal.x << " " << goal.y << " " << goal.z << std::endl;
		corr[i] = (goal - pos[i]) * shapeMatchingStiffness;
	}

	if (rot)
		*rot = R;


	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		if (physModel->hInvMass[i] != 0.0)
			pos[i] += corr[i]; //现在没有使用cluster
	}

	return true;
}

bool PBDSystem::HostUpdateCubeVolumeConstraints(unsigned idx)
{
	int i_center = hCubeVolumeConstraints[idx].a;
	int i_xplus = hCubeVolumeConstraints[idx].b;
	int i_xminus = hCubeVolumeConstraints[idx].c;
	int i_yplus = hCubeVolumeConstraints[idx].d;
	int i_yminus = hCubeVolumeConstraints[idx].e;
	int i_zplus = hCubeVolumeConstraints[idx].f;
	int i_zminus = hCubeVolumeConstraints[idx].g;

	double3 &vert_center = physModel->hPredictPos[i_center];
	double3 &vert_xplus = physModel->hPredictPos[i_xplus];
	double3 &vert_xminus = physModel->hPredictPos[i_xminus];
	double3 &vert_yplus = physModel->hPredictPos[i_yplus];
	double3 &vert_yminus = physModel->hPredictPos[i_yminus];
	double3 &vert_zplus = physModel->hPredictPos[i_zplus];
	double3 &vert_zminus = physModel->hPredictPos[i_zminus];

	double invMass_center = physModel->hInvMass[i_center];
	double invMass_xplus = physModel->hInvMass[i_xplus];
	double invMass_xminus = physModel->hInvMass[i_xminus];
	double invMass_yplus = physModel->hInvMass[i_yplus];
	double invMass_yminus = physModel->hInvMass[i_yminus];
	double invMass_zplus = physModel->hInvMass[i_zplus];
	double invMass_zminus = physModel->hInvMass[i_zminus];

	double cubeVolume = YH::Dot(vert_xplus, YH::Cross(vert_yplus, vert_zplus)) + YH::Dot(vert_xplus, YH::Cross(vert_yminus, vert_zminus)) +
		YH::Dot(vert_xplus, YH::Cross(vert_yplus, vert_zminus)) + YH::Dot(vert_xminus, YH::Cross(vert_yminus, vert_zplus)) +
		YH::Dot(vert_zplus, YH::Cross(vert_yplus, vert_xminus)) + YH::Dot(vert_zplus, YH::Cross(vert_yminus, vert_xplus)) +
		YH::Dot(vert_zminus, YH::Cross(vert_yplus, vert_xplus)) + YH::Dot(vert_zminus, YH::Cross(vert_yminus, vert_xminus));
	
	double3 grad_zminus = YH::Cross(-vert_xplus, vert_yminus) - YH::Cross(vert_xminus, vert_yplus) + YH::Cross(vert_yplus, vert_xplus) + YH::Cross(vert_yminus, vert_xminus);
	double3 grad_zplus = YH::Cross(-vert_xplus, vert_yplus) - YH::Cross(vert_xminus, vert_yminus) + YH::Cross(vert_yplus, vert_xminus) + YH::Cross(vert_yminus, vert_xplus);
	double3 grad_yminus = YH::Cross(vert_xplus, vert_zminus) + YH::Cross(vert_xminus, vert_zplus) + YH::Cross(vert_zplus, vert_xplus) + YH::Cross(vert_zminus, vert_xminus);
	double3 grad_yplus = YH::Cross(vert_xplus, vert_zplus) + YH::Cross(vert_xminus, vert_zminus) + YH::Cross(vert_zplus, vert_xminus) + YH::Cross(vert_zminus, vert_xplus);
	double3 grad_xminus = YH::Cross(vert_yplus, vert_zminus) + YH::Cross(vert_yminus, vert_zplus) - YH::Cross(vert_zplus, vert_yplus) - YH::Cross(vert_zminus, vert_yminus);
	double3 grad_xplus = -(grad_zminus + grad_zplus + grad_yminus + grad_yplus + grad_xminus);

	double lambda = invMass_xplus * YH::SquareLength(grad_xplus) + invMass_xminus * YH::SquareLength(grad_xminus) +
					invMass_yplus * YH::SquareLength(grad_yplus) + invMass_yminus * YH::SquareLength(grad_yminus) +
					invMass_zplus * YH::SquareLength(grad_zplus) + invMass_zminus * YH::SquareLength(grad_zminus);

	if (fabs(lambda) < eps)
		return false;

	if(cubeVolume < 0.0)
		lambda = negCubeVolumeStiffness * (cubeVolume - hRestCubeVolume[idx]) / lambda;
	else
		lambda = posCubeVolumeStiffness * (cubeVolume - hRestCubeVolume[idx]) / lambda;

	double3 corr_xplus, corr_xminus, corr_yplus, corr_yminus, corr_zplus, corr_zminus;

	corr_xplus = -invMass_xplus * lambda * grad_xplus;
	corr_xminus = -invMass_xminus * lambda * grad_xminus;
	corr_yplus = -invMass_yplus * lambda * grad_yplus;
	corr_yminus = -invMass_yminus * lambda * grad_yminus;
	corr_zplus = -invMass_zplus * lambda * grad_zplus;
	corr_zminus = -invMass_zminus * lambda * grad_zminus;

	if (invMass_xplus > 0.0)
	{
		//vert_xplus += corr_xplus;
		hCubeVolumeCorrSum[i_xplus] += corr_xplus;
	}
	if (invMass_xminus > 0.0)
	{
		vert_xminus += corr_xminus;
		//hCubeVolumeCorrSum[i_xminus] += corr_xminus;
	}
	if (invMass_yplus > 0.0)
	{
		vert_yplus += corr_yplus;
		//hCubeVolumeCorrSum[i_yplus] += corr_yplus;
	}
	if (invMass_yminus > 0.0)
	{
		vert_yminus += corr_yminus;
		//hCubeVolumeCorrSum[i_yminus] += corr_yminus;
	}
	if (invMass_zplus > 0.0)
	{
		vert_zplus += corr_zplus;
		//hCubeVolumeCorrSum[i_zplus] += corr_zplus;
	}
	if (invMass_zminus > 0.0)
	{
		vert_zminus += corr_zminus;
		//hCubeVolumeCorrSum[i_zminus] += corr_zminus;
	}

	return true;
}

void PBDSystem::HostUpdateAllConstraints()
{
	for (unsigned i = 0; i < hNumDistanceConstraints; i++)
	{
		HostUpdateDistanceConstraints(i);
	}
	/*for (unsigned i = 0; i > physModel->numParticles; i++)
	{
		uint neighborCount = physModel->hNeighborCount[i];
		if (neighborCount > 0)
			physModel->hPredictPos[i] += (1.0 / neighborCount) * hDistanceCorrSum[i];
	}*/
	/*for (unsigned i = 0; i < hNumVolumeConstraints; i++)
	{
		HostUpdateVolumeConstraints(i);
	}*/
	/*for (unsigned i = 0; i < hNumCubeVolumeConstraints; i++)
	{
		HostUpdateCubeVolumeConstraints(i);
	}*/
	/*for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		physModel->hPredictPos[i] += (1.0 / 6.0) * hCubeVolumeCorrSum[i];
	}*/

	//HostUpdateShapeMatchingConstraints(false);
}

void PBDSystem::HostIntegrate(double dt)
{
	glBindBuffer(GL_ARRAY_BUFFER, physModel->vertices->m_positionBuffer);
	double3 *pos = (double3*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE);

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		physModel->hVel[i] = (physModel->hPredictPos[i] - pos[i]) / dt;
		pos[i] = physModel->hPredictPos[i];
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
}

void PBDSystem::HostDoPBD(cVector3d toolPosInModelSpace, cVector3d toolForceInModelSpace, bool applyForceOnPhysModel)
{
	toolPos.x = toolPosInModelSpace.x(); toolPos.y = toolPosInModelSpace.y(); toolPos.z = toolPosInModelSpace.z();
	//toolForce.x = toolForceInModelSpace.x(); toolForce.y = toolForceInModelSpace.y(); toolForce.z = toolForceInModelSpace.z();
	if (applyForceOnPhysModel)
	{
		toolForce.x = 0; toolForce.y = -9.8; toolForce.z = 0;
	}
	else
	{
		toolForce.x = 0; toolForce.y = 0; toolForce.z = 0;
	}

	//host find nearby particles
	HostFindParticlesWithinDstFromTarget(toolPos, physModel->physParam.Solid_1_EffectiveRadius * 5);

	for (unsigned i = 0; i < physModel->numParticles; i++)
	{
		physModel->vertices->setColor(i, cColorb(255, 255, 255));
	}
	for (unsigned j = 0; j < physModel->hNearbyCount; j++)
	{
		physModel->vertices->setColor(physModel->hNearbyIdx[j], cColorb(255, 0, 0));
	}
	

	//host update velocity
	HostUpdateVelocities(dt);

	//host update all constraints
	for (unsigned i = 0; i < 1; i++)
	{
		HostUpdateAllConstraints();
	}

	//host Integrate
	HostIntegrate(dt);

	physModel->hNearbyIdx.clear();
	physModel->hNearbyCount = 0;
}