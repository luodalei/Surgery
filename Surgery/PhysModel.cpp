#include "PhysModel.h"

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
}


PhysModel::~PhysModel()
{
}
