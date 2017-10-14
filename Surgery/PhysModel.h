#pragma once
#include "Utilities.h"

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


class PhysModel
{
public:
	PhysModel();
	~PhysModel();

public:
	std::vector<double> masses;
	std::vector<glm::vec3> positions;
	std::vector<glm::vec3> velocities;
	std::vector<double> temperatures;
	std::vector<double> lastTemperatures;
	std::vector<ParticleType> particleTypes;
	std::vector<double> heatQ;
};

