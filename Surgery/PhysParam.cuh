#ifndef PHYSPARAM_H
#define PHYSPARAM_H

#include "helper_math.h"
#include <cuda_runtime.h>

#include <float.h>

typedef struct _PhysParam
{
	//四种状态粒子初始质量
	double SOLID_1_Mass;
	double SOLID_2_Mass;
	double LIQUID_Mass;
	double GAS_Mass;

	double Solid_Density;
	double Liquid_Density;
	double Gas_Density;

	double Air_Temperature;

	double Split_Temperature; //SOLID_1 -> SOLID_2
	double Fusion_Temperature;//SOLID_2 -> LIQUID
	double Boil_Temperature; //LIQUID -> GAS
	double LatentHeat_Solid2Liquid;//SOLID_2 -> LIQUID
	double LatentHeat_Liquid2Gas;

	double Thermal_Conductivity; //h: Qi=h(Tair - Ti)*A
	double Thermal_Diffusion; //Cd： in heat transfer between particle
	double Heat_Capacity_Solid; //C: delta_T=Qi/C*m
	double Heat_Capacity_Liquid;
	double Heat_Capacity_Gas;

	//光子参数
	double Emissivity;
	double Boltzmann_Constant;
	double Source_Temperature;
	double Source_Area;
	int PhotonsNum_TimeStep;//单位时间射出的光子数

	//各种半径参数
	double Solid_1_EffectiveRadius; //re
	double Solid_2_EffectiveRadius;

	//cuda neighbor find
	uint3 gridSize; //64
	uint numCells; 
	double3 worldOrigin;
	double3 cellSize; // cell size equal to particle effective radius * 2
	uint numBodies;
	uint maxParticlesPerCell;//没用上

} PhysParam;


template <typename T>
struct _Mat3 {

public:

	__host__ __device__ _Mat3() {}

	_Mat3(T t1, T t2, T t3)
	{
		m_data[0][0] = t1; m_data[0][1] = 0; m_data[0][2] = 0;//row0
		m_data[1][0] = 0; m_data[1][1] = t2; m_data[1][2] = 0;
		m_data[2][0] = 0; m_data[2][1] = 0; m_data[2][2] = t3;
	}

	__host__ __device__ _Mat3(T t1, T t2, T t3,
							T t4, T t5, T t6,
							T t7, T t8, T t9)
	{
		m_data[0][0] = t1; m_data[0][1] = t2; m_data[0][2] = t3;
		m_data[1][0] = t4; m_data[1][1] = t5; m_data[1][2] = t6;
		m_data[2][0] = t7; m_data[2][1] = t8; m_data[2][2] = t9;
	}

	__host__ __device__ inline void SetZero()
	{
		m_data[0][0] = 0; m_data[0][1] = 0; m_data[0][2] = 0;
		m_data[1][0] = 0; m_data[1][1] = 0; m_data[1][2] = 0;
		m_data[2][0] = 0; m_data[2][1] = 0; m_data[2][2] = 0;
	}

	//行列式
	__host__ __device__ inline double Det() const
	{
		return (+(*this)(0, 0) * (*this)(1, 1) * (*this)(2, 2)
			+ (*this)(0, 1) * (*this)(1, 2) * (*this)(2, 0)
			+ (*this)(0, 2) * (*this)(1, 0) * (*this)(2, 1)
			- (*this)(2, 0) * (*this)(1, 1) * (*this)(0, 2)
			- (*this)(2, 1) * (*this)(1, 2) * (*this)(0, 0)
			- (*this)(2, 2) * (*this)(1, 0) * (*this)(0, 1));
	}

	__host__ __device__ inline void Mulr(const double3& vec, double3& res)
	{
		res.x = (*this)(0, 0) * vec.x + (*this)(0, 1) * vec.y + (*this)(0, 1) * vec.z;
		res.y = (*this)(1, 0) * vec.x + (*this)(1, 1) * vec.y + (*this)(1, 2) * vec.z;
		res.z = (*this)(2, 0) * vec.x + (*this)(2, 1) * vec.y + (*this)(2, 2) * vec.z;
	}

	__host__ __device__ inline void Mulr(const _Mat3& a_matrix, _Mat3& a_result) const
	{
		// compute multiplication between both matrices
		a_result(0, 0) = (*this)(0, 0) * a_matrix(0, 0) + (*this)(0, 1) * a_matrix(1, 0) + (*this)(0, 2) * a_matrix(2, 0);
		a_result(0, 1) = (*this)(0, 0) * a_matrix(0, 1) + (*this)(0, 1) * a_matrix(1, 1) + (*this)(0, 2) * a_matrix(2, 1);
		a_result(0, 2) = (*this)(0, 0) * a_matrix(0, 2) + (*this)(0, 1) * a_matrix(1, 2) + (*this)(0, 2) * a_matrix(2, 2);
		a_result(1, 0) = (*this)(1, 0) * a_matrix(0, 0) + (*this)(1, 1) * a_matrix(1, 0) + (*this)(1, 2) * a_matrix(2, 0);
		a_result(1, 1) = (*this)(1, 0) * a_matrix(0, 1) + (*this)(1, 1) * a_matrix(1, 1) + (*this)(1, 2) * a_matrix(2, 1);
		a_result(1, 2) = (*this)(1, 0) * a_matrix(0, 2) + (*this)(1, 1) * a_matrix(1, 2) + (*this)(1, 2) * a_matrix(2, 2);
		a_result(2, 0) = (*this)(2, 0) * a_matrix(0, 0) + (*this)(2, 1) * a_matrix(1, 0) + (*this)(2, 2) * a_matrix(2, 0);
		a_result(2, 1) = (*this)(2, 0) * a_matrix(0, 1) + (*this)(2, 1) * a_matrix(1, 1) + (*this)(2, 2) * a_matrix(2, 1);
		a_result(2, 2) = (*this)(2, 0) * a_matrix(0, 2) + (*this)(2, 1) * a_matrix(1, 2) + (*this)(2, 2) * a_matrix(2, 2);
	}

	__host__ __device__ inline T& operator() (const unsigned idx1, const unsigned idx2)
	{
		return m_data[idx1][idx2];
	}

	__host__ __device__ inline const T& operator() (const unsigned idx1, const unsigned idx2) const
	{
		return m_data[idx1][idx2];
	}

	__host__ __device__ inline double3 operator* (const double3& vec)
	{
		double3 res;
		this->Mulr(vec, res);
		return res;
	}

	__host__ __device__ inline _Mat3 operator* (const _Mat3& mat)
	{
		_Mat3 res;
		this->Mulr(mat, res);
		return res;
	}

	__host__ __device__ inline void operator+= (const _Mat3<T>& src)
	{
		(*this)(0, 0) += src(0, 0);
		(*this)(0, 1) += src(0, 1);
		(*this)(0, 2) += src(0, 2);

		(*this)(1, 0) += src(1, 0);
		(*this)(1, 1) += src(1, 1);
		(*this)(1, 2) += src(1, 2);

		(*this)(2, 0) += src(2, 0);
		(*this)(2, 1) += src(2, 1);
		(*this)(2, 2) += src(2, 2);
	}

	__host__ __device__ inline bool Invertr(_Mat3<T>& res) const
	{
		// compute determinant
		T d = (+(*this)(0, 0) * (*this)(1, 1) * (*this)(2, 2)
				+ (*this)(0, 1) * (*this)(1, 2) * (*this)(2, 0)
				+ (*this)(0, 2) * (*this)(1, 0) * (*this)(2, 1)
				- (*this)(2, 0) * (*this)(1, 1) * (*this)(0, 2)
				- (*this)(2, 1) * (*this)(1, 2) * (*this)(0, 0)
				- (*this)(2, 2) * (*this)(1, 0) * (*this)(0, 1));

		// check if determinant null.
		if ((d < FLT_MIN) && (d > -FLT_MIN))
		{
			// determinant null, matrix inversion can not be performed
			return (false);
		}

		// compute inverted matrix
		res(0, 0) = ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1)*(*this)(1, 2)) / d;
		res(0, 1) = -((*this)(0, 1) * (*this)(2, 2) - (*this)(2, 1)*(*this)(0, 2)) / d;
		res(0, 2) = ((*this)(0, 1) * (*this)(1, 2) - (*this)(1, 1)*(*this)(0, 2)) / d;

		res(1, 0) = -((*this)(1, 0) * (*this)(2, 2) - (*this)(2, 0)*(*this)(1, 2)) / d;
		res(1, 1) = ((*this)(0, 0) * (*this)(2, 2) - (*this)(2, 0)*(*this)(0, 2)) / d;
		res(1, 2) = -((*this)(0, 0) * (*this)(1, 2) - (*this)(1, 0)*(*this)(0, 2)) / d;

		res(2, 0) = ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0)*(*this)(1, 1)) / d;
		res(2, 1) = -((*this)(0, 0) * (*this)(2, 1) - (*this)(2, 0)*(*this)(0, 1)) / d;
		res(2, 2) = ((*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0)*(*this)(0, 1)) / d;

		// return success
		return (true);
	}

	__host__ __device__ inline void Transr(_Mat3<T>& res) const
	{
		res(0, 0) = (*this)(0, 0);
		res(0, 1) = (*this)(1, 0);
		res(0, 2) = (*this)(2, 0);

		res(1, 0) = (*this)(0, 1);
		res(1, 1) = (*this)(1, 1);
		res(1, 2) = (*this)(2, 1);

		res(2, 0) = (*this)(0, 2);
		res(2, 1) = (*this)(1, 2);
		res(2, 2) = (*this)(2, 2);
	}

	__host__ __device__ inline void Identity()
	{
		(*this)(0, 0) = 1.0;  (*this)(0, 1) = 0.0;  (*this)(0, 2) = 0.0;
		(*this)(1, 0) = 0.0;  (*this)(1, 1) = 1.0;  (*this)(1, 2) = 0.0;
		(*this)(2, 0) = 0.0;  (*this)(2, 1) = 0.0;  (*this)(2, 2) = 1.0;
	}

	__host__ __device__ inline double3 GetRow(const unsigned i)
	{
		double3 res;
		res.x = (*this)(i, 0);
		res.y = (*this)(i, 1);
		res.z = (*this)(i, 2);

		return res;
	}

	__host__ __device__ inline void SetRow(const unsigned i, const double3& newRow)
	{
		(*this)(i, 0) = newRow.x;
		(*this)(i, 1) = newRow.y;
		(*this)(i, 2) = newRow.z;
	}

public:
	T m_data[3][3];
};

typedef _Mat3<double> Mat3d;

__host__ __device__ inline double3 operator- (const double3& vec1, const double3& vec2)
{
	double3 res;

	res.x = vec1.x - vec2.x;
	res.y = vec1.y - vec2.y;
	res.z = vec1.z - vec2.z;

	return res;
}

__host__ __device__ inline double3 operator+ (const double3& vec1, const double3& vec2)
{
	double3 res;

	res.x = vec1.x + vec2.x;
	res.y = vec1.y + vec2.y;
	res.z = vec1.z + vec2.z;

	return res;
}

__host__ __device__ inline void operator+= (double3& vec1, const double3& vec2)
{
	vec1.x += vec2.x;
	vec1.y += vec2.y;
	vec1.z += vec2.z;
}

__host__ __device__ inline void operator-= (double3& vec1, const double3& vec2)
{
	vec1.x -= vec2.x;
	vec1.y -= vec2.y;
	vec1.z -= vec2.z;
}


__host__ __device__ inline double3 operator* (const double scale, const double3& vec)
{
	double3 res;
	res.x = vec.x * scale;
	res.y = vec.y * scale;
	res.z = vec.z * scale;

	return res;
}

__host__ __device__ inline double3 operator* (const double3& vec, const double scale)
{
	double3 res;
	res.x = vec.x * scale;
	res.y = vec.y * scale;
	res.z = vec.z * scale;

	return res;
}

__host__ __device__ inline void operator*= (double3& vec, const double scale)
{
	vec.x = vec.x * scale;
	vec.y = vec.y * scale;
	vec.z = vec.z * scale;
}

__host__ __device__ inline void operator/= (double3& vec, const double scale)
{
	vec.x = vec.x / scale;
	vec.y = vec.y / scale;
	vec.z = vec.z / scale;
}

__host__ __device__ inline double3 operator/ (const double3& vec, const double scale)
{
	double3 res;
	res.x = vec.x / scale;
	res.y = vec.y / scale;
	res.z = vec.z / scale;

	return res;
}

__host__ __device__ inline double3 operator-(const double3& vec)
{
	double3 res;
	res.x = -vec.x;
	res.y = -vec.y;
	res.z = -vec.z;
	return res;
}

namespace YH
{
	__host__ __device__ inline Mat3d Invert(const Mat3d& src)
	{
		Mat3d res;
		src.Invertr(res);
		return res;
	}

	
	__host__ __device__ inline Mat3d Transpose(const Mat3d& src)
	{
		Mat3d res;
		src.Transr(res);
		return res;
	}

	__host__ __device__ inline Mat3d Identity()
	{
		Mat3d res(1.0, 1.0, 1.0);
		return res;
	}

	__host__ __device__ inline double Dot(const double3& vec1, const double3& vec2)
	{
		return (vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z);
	}

	__host__ __device__ inline double3 Cross(const double3& vec1, const double3& vec2)
	{
		double3 res;

		// compute cross product
		double a = (vec1.y * vec2.z) - (vec1.z * vec2.y);
		double b = -(vec1.x * vec2.z) + (vec1.z * vec2.x);
		double c = (vec1.x * vec2.y) - (vec1.y * vec2.x);

		// store result in current vector
		res.x = a;
		res.y = b;
		res.z = c;

		return res;
	}

	__host__ __device__ inline double Length(const double3& vec)
	{
		return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}

	__host__ __device__ inline double SquareLength(const double3& vec)
	{
		return (vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
	}

	__host__ __device__ inline double Distance(const double3& vec1, const double3& vec2)
	{
		return (vec1.x - vec2.x) * (vec1.x - vec2.x) + (vec1.y - vec2.y) * (vec1.y - vec2.y) + (vec1.z - vec2.z) * (vec1.z - vec2.z);
	}

	__host__ __device__ inline double Distanceq(const double3& vec1, const double3& vec2)
	{
		return sqrt(Distance(vec1, vec2));
	}

	__host__ __device__ inline double3 Normalize(const double3& vec)
	{
		double3 res;

		double len = Length(vec);
		if (len == 0.0)
		{
			res.x = res.y = res.z = FLT_MAX;
		}
		else
		{
			double factor = 1 / len;
			res = vec * factor;
		}
		return res;
	}

	// Return the one norm of the matrix.
	__host__ __device__ inline double OneNorm(const Mat3d &A)
	{
		const double sum1 = fabs(A(0, 0)) + fabs(A(1, 0)) + fabs(A(2, 0));
		const double sum2 = fabs(A(0, 1)) + fabs(A(1, 1)) + fabs(A(2, 1));
		const double sum3 = fabs(A(0, 2)) + fabs(A(1, 2)) + fabs(A(2, 2));

		double maxSum = sum1;
		if (sum2 > maxSum)
			maxSum = sum2;
		if (sum3 > maxSum)
			maxSum = sum3;
		return maxSum;
	}

	// Return the inf norm of the matrix.
	__host__ __device__ inline double InfNorm(const Mat3d &A)
	{
		const double sum1 = fabs(A(0, 0)) + fabs(A(0, 1)) + fabs(A(0, 2));
		const double sum2 = fabs(A(1, 0)) + fabs(A(1, 1)) + fabs(A(1, 2));
		const double sum3 = fabs(A(2, 0)) + fabs(A(2, 1)) + fabs(A(2, 2));

		double maxSum = sum1;
		if (sum2 > maxSum)
			maxSum = sum2;
		if (sum3 > maxSum)
			maxSum = sum3;
		return maxSum;
	}

	// Perform a polar decomposition of matrix M and return the rotation matrix R. This method handles the degenerated cases.
	__host__ __device__ inline void PolarDecompositionStable(const Mat3d &M, const double  tolerance, Mat3d &R)
	{
		Mat3d Mt; M.Transr(Mt);
		double Mone = OneNorm(M);
		double Minf = InfNorm(M);
		double Eone;
		Mat3d MadjTt, Et;

		do
		{
			MadjTt.SetRow(0, YH::Cross(Mt.GetRow(1), Mt.GetRow(2)));
			MadjTt.SetRow(1, YH::Cross(Mt.GetRow(2), Mt.GetRow(0)));
			MadjTt.SetRow(2, YH::Cross(Mt.GetRow(0), Mt.GetRow(1)));

			double det = Mt(0, 0) * MadjTt(0, 0) + Mt(0, 1) * MadjTt(0, 1) + Mt(0, 2) * MadjTt(0, 2);

			if (fabs(det) < 1.0e-12)
			{
				double3 len;
				unsigned index = 0xffffffff;
				
				len.x = YH::SquareLength(MadjTt.GetRow(0));
				len.y = YH::SquareLength(MadjTt.GetRow(1));
				len.z = YH::SquareLength(MadjTt.GetRow(2));

				if (len.x > 1.0e-12)
					index = 0;
				else if (len.y > 1.0e-12)
					index = 1;
				else if (len.z > 1.0e-12)
					index = 2;
				
				if (index == 0xffffffff)
				{
					R.Identity();
					return;
				}
				else
				{
					Mt.SetRow(index, YH::Cross(Mt.GetRow((index + 1) % 3), Mt.GetRow((index + 2) % 3)));
					MadjTt.SetRow((index + 1) % 3, YH::Cross(Mt.GetRow((index + 2) % 3), Mt.GetRow(index)));
					MadjTt.SetRow((index + 2) % 3, YH::Cross(Mt.GetRow(index), Mt.GetRow((index + 1) % 3)));
					
					Mat3d M2; Mt.Transr(M2);
					Mone = OneNorm(M2);
					Minf = InfNorm(M2);
					det = Mt(0, 0) * MadjTt(0, 0) + Mt(0, 1) * MadjTt(0, 1) + Mt(0, 2) * MadjTt(0, 2);
				}
			}

			const double MadjTone = OneNorm(MadjTt);
			const double MadjTinf = InfNorm(MadjTt);

			const double gamma = sqrt(sqrt((MadjTone * MadjTinf) / (Mone * Minf)) / fabs(det));

			const double g1 = gamma * 0.5;
			const double g2 = 0.5 / (gamma * det);

			for (unsigned i = 0; i < 3; i++)
			{
				for (unsigned j = 0; j < 3; j++)
				{
					Et(i, j) = Mt(i, j);
					Mt(i, j) = g1 * Mt(i, j) + g2 * MadjTt(i, j);
					Et(i, j) -= Mt(i, j);
				}
			}

			Eone = OneNorm(Et);

			Mone = OneNorm(Mt);
			Minf = InfNorm(Mt);

		} while (Eone > Mone * tolerance);

		//Q = Mt^T
		Mt.Transr(R);
	}
}


#endif