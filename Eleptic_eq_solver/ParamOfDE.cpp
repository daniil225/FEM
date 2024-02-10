#include "ParamOfDE.h"

Func_3D GetLambda(const ParamOfDE& param, int num) { return param.Lambda[num]; }


Func_3D GetGamma(const ParamOfDE& param, int num) { return param.Gamma[num]; }

Func_3D GetF(const ParamOfDE& param, int num) { return param.F[num]; }

Func_3D GetBoundCond12(const ParamOfDE& param, int num)
{
	std::pair<int, int> idx = param.Map[num];

	if (idx.second == 1)
	{
		return param.FirstBoundCond[idx.first];
	}
	else
		return param.SecondBoundCond[idx.first];
}

PairFunc_3D GetBoundCond3(const ParamOfDE& param, int num)
{
	std::pair<int, int> idx = param.Map[num];

	//if (idx.second == 3)
	return param.ThirdBoundCond[idx.first];
}

ParamOfDE InitTest1()
{
	ParamOfDE Test;
	Test.Lambda.resize(1);
	Test.Lambda[0] = [](double x, double y, double z) -> double { return 1; };

	Test.Gamma.resize(1);
	Test.Gamma[0] = [](double x, double y, double z) -> double { return 0; };

	Test.F.resize(1);
	Test.F[0] = [](double x, double y, double z) -> double {return -2-2*sin(z) + y*y*sin(z); };

	Test.FirstBoundCond.resize(3);
	Test.SecondBoundCond.resize(3);

	Test.SecondBoundCond[0] = [](double x, double y, double z) -> double { return -2*x; };
	Test.FirstBoundCond[0] = [](double x, double y, double z) -> double { return 1 + y*y*sin(z); };
	Test.FirstBoundCond[1] = [](double x, double y, double z) -> double { return x*x; };
	Test.SecondBoundCond[1] = [](double x, double y, double z) -> double { return 2*y*sin(z); };
	Test.FirstBoundCond[2] = [](double x, double y, double z) -> double { return x*x; };
	Test.SecondBoundCond[2] = [](double x, double y, double z) -> double { return y*y*cos(z); };


	Test.Map.resize(6);

	Test.Map[0].first = 0; // Индекс формулы в соответсвующем массиве 
	Test.Map[0].second = 2; // Массив первых краевых 

	Test.Map[1].first = 0; // Индекс формулы в соответсвующем массиве 
	Test.Map[1].second = 1; // Массив первых краевых

	Test.Map[2].first = 1; // Индекс формулы в соответсвующем массиве 
	Test.Map[2].second = 1; // Массив первых краевых
	
	Test.Map[3].first = 1; // Индекс формулы в соответсвующем массиве 
	Test.Map[3].second = 2; // Массив первых краевых
	
	Test.Map[4].first = 2; // Индекс формулы в соответсвующем массиве 
	Test.Map[4].second = 1; // Массив первых краевых
	
	Test.Map[5].first = 2; // Индекс формулы в соответсвующем массиве 
	Test.Map[5].second = 2; // Массив первых краевых

	return Test;
}