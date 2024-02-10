#pragma once

#include <functional>
#include <tuple>
#include <vector>

/* Парамеьры ДУ упакованные в класс 
	Для каждой области настраиваем свой экземпляр 
*/

typedef std::function < double(double, double, double)> Func_3D;
typedef std::vector<Func_3D> ArrayFunc_3D;
typedef std::pair<Func_3D, Func_3D> PairFunc_3D;
typedef std::vector<PairFunc_3D> ArrayPairFunc_3D;
typedef std::vector<std::pair<int,int>> VecPairInt;
/* Объединение параметров в одно*/
struct ParamOfDE
{
	ArrayFunc_3D Lambda; // Значение параметра Lambda
	ArrayFunc_3D Gamma; // Значение параметре Гамма
	ArrayFunc_3D F; // Правая часть 
	
	/* Кривоватое решение но введем массив соответсвий для КУ */

	// Это вспомогательная структура в ней сопоставляются номера формул и номера формул в массивах 1-х КУ 2-х КУ 3-х КУ 
	// индекс это набор формул значение в ячейке индекс в соответсвующем массиве а так же массиву в котором вся информация есть 
	VecPairInt Map;  

	ArrayFunc_3D FirstBoundCond;
	ArrayFunc_3D SecondBoundCond;
	ArrayPairFunc_3D ThirdBoundCond;
};

/* получить функцию лямбда 
* @param
* const ParamOfDE &param - параметры ДУ
* int num - номер формул
*/
Func_3D GetLambda(const ParamOfDE &param,int num);

/* получить функцию Гамма
* @param
* const ParamOfDE &param - параметры ДУ
* int num - номер формул
*/
Func_3D GetGamma(const ParamOfDE& param, int num);

/* получить функцию правой части 
* @param
* const ParamOfDE &param - параметры ДУ
* int num - номер формул
*/
Func_3D GetF(const ParamOfDE& param, int num);

/* получить функцию КУ 1 и 2 ого типов 
* @param
* const ParamOfDE &param - параметры ДУ
* int num - номер формул
*/
Func_3D GetBoundCond12(const ParamOfDE& param, int num);


/* получить пару функций 3-х КУ
* @param
* const ParamOfDE &param - параметры ДУ
* int num - номер формул
*/
PairFunc_3D GetBoundCond3(const ParamOfDE& param, int num);

/* Инициализация первого теста
* @param
* void
* ret ParamOfDE - параметры ДУ
*/
ParamOfDE InitTest1();
