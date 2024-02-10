#pragma once

/* Вроде как тоже просматривается общая структура для построения матриц МКЭ возможно потом перепишу в ООП стилистике */

/* Функционал для сборки и обработки матрицы конечно элментной СЛАУ */
#include "SlauSolver.h"
#include "Grid.h"
#include "ParamOfDE.h"


namespace FormMatrix
{
	using namespace SLAUSolvers::IterationSolvers;
	using namespace Grid;


	/*	Генерация конечноэлементной СЛАУ
	* @param
	* SLAU& slau - СЛАУ
	* const GridS &grid - Конечноэлементная сетка
	* const ParamOfDE & Param - Параметры определяющие решение ДУ 
	* ret: void
	* 
	* Сгенерированная матрица хранится в разряженном строчностолбцовом формате как !!! НЕ СИММЕТРИЧНАЯ - соответсвенно нужно выбирать соответсвующий решатель  
	*/
	void GenerateSLAU(SLAU& slau, const GridS& grid, const ParamOfDE & Param);

	namespace PrivateFormMatrix
	{

		/* Составные части матриц */
		struct LocalPartMatrix
		{
			/* Коэффициенты для вычисления якобиана и обратного преобразования функций */
			double sign_a0;
			double a0;
			double a1;
			double a2;

			/* Коэффициенты появляются в матрице Жесткости нужный для вычисления ее компонент*/
			double b1;
			double b2;
			double b3;
			double b4;
			double b5;
			double b6;

			/* Произфодные базисных функций */
			std::function<double(double)> dPhi_dKsi[4] =
			{
				[](double Eth) -> double { return Eth - 1; },
				[](double Eth) -> double { return 1 - Eth; },
				[](double Eth) -> double { return -Eth;    },
				[](double Eth) -> double { return Eth;     }
			};

			std::function<double(double)> dPhi_dEth[4] = 
			{
				[](double Ksi) -> double { return Ksi - 1; },
				[](double Ksi) -> double { return -Ksi;    },
				[](double Ksi) -> double { return 1 - Ksi; },
				[](double Ksi) -> double { return Ksi;     }
			};


			std::function<double(double, double)> J; // Якобиан

			/* Составные части матрицы массы */
			double M_xz0[4][4] =
			{ 
				{4.0,2.0,2.0,1.0},
				{2.0,4.0,1.0,2.0},
				{2.0,1.0,4.0,2.0},
				{1.0,2.0,2.0,4.0} 
			};

			double M_xz1[4][4] = {
				{2.0,2.0,1.0,1.0},
				{2.0,6.0,1.0,3.0},
				{1.0,1.0,2.0,2.0},
				{1.0,3.0,2.0,6.0}
			};

			double M_xz2[4][4] = {
				{2.0,1.0,2.0,1.0},
				{1.0,2.0,1.0,2.0},
				{2.0,1.0,6.0,3.0},
				{1.0,2.0,3.0,6.0}
			};

			// Собранная матрица массы из частей M = sign(a0)*(a0/36*M_xz0 + a1/72*M_xz1 + a2/72*M_xz2)
			// a0, a1, a2 зависят от матрицы 
			double M_xz[4][4];

			/* Генерация матрицы производится численным вычислением интеграла  */
			double C_xz[4][4];

			double M_y[2][2] = {
				{2,1},
				{1,2}
			};

			double C_y[2][2] = {
				{1,-1},
				{-1,1}
			};
		};


		/* Локальная матрица и вектор F */
		struct Local
		{
			double LocalMassMatrix[Element::FinitElementSize][Element::FinitElementSize]; // Матрица массы для заданного конечного элемента по ней можно расчатать все необходимые вклады 
			double LocalStiffnessMatrix[Element::FinitElementSize][Element::FinitElementSize]; // Матрица жесткости 

			double LocalMatrix[Element::FinitElementSize][Element::FinitElementSize]; // Результирующая матрица Которую нужно добавлять в Глобальную 
			double LocalF[Element::FinitElementSize]; // Локальный вектор F в нем сразу учет 2-х и 3-х КУ + учет правой части в ДУ
			int LocalToGlobal[Element::FinitElementSize]; // Отображение локальной индексации в Глобальную 
		};
		
		/*	Генерация портрета матрицы СЛАУ
		* @param
		* SLAU& slau - СЛАУ
		* const GridS &grid - Конечноэлементная сетка
		* ret: void
		*/
		void GeneratePortrait(SLAU& slau, const GridS& grid);

		/* Добавление локальной матрицы и вектора в Глобальную
		* @param
		* SLAU& slau - Глобальная матрица и ветор правой части 
		* const Local& local - локальный элемент 
		*  ret: void 
		* Результат добавили локальную матрицу и вектор f в Глобальную 
		*/
		void AddLocal(SLAU& slau, const Local& local);
		
		/* Занесение элемента элемента а в глобальную матрицу с индеком Aij
		* @param
		* SLAU& slau - Матрица куда заносим 
		* int i - индекс 
		* int j - индекс
		* double a - заносимое значение 
		* ret: void 
		* Результат Элемент занесен в матрицу СЛАУ 
		*/
		void AddElement(SLAU& slau, int i, int j, double a);

		/* Устанавливает элемент в в глобальную матрицу с индексом Aij
		* @param
		* SLAU& slau - Матрица куда заносим 
		* int i - индекс 
		* int j - индекс
		* double a - заносимое значение 
		* ret: void 
		* Результат Элемент установлен в матрицу СЛАУ 
		*/
		void SetElement(SLAU& slau, int i, int j, double a);

		/* Генерация локального элемента
		* Local& local - локальный элемент 
		* const Element &element - Конечный элемент содержит в себе всю инфомацию о границах + всю информацию об области
		* const ParamOfDE& Param - параметры ДУ 
		* ret: void 
		* Результат Сгенерированная Локальная матрица и ветор правой части по Условиям уравнения + КУ 2 и 3 его рода 
		*/
		void GenerateLocal(Local& local, const Element &element,const ParamOfDE& Param);

		/* Учет 1-х КУ */
		void AddMainCond(SLAU& slau, const Element& element, const ParamOfDE& Param);

		/* Генерация составных частей локальных матриц для заданного конечного элемента
		* @param
		* const Element& element - конечный элемент 
		* ret: LocalPartMatrix - Локальные части матриц для сборки общих локальных матриц 
		*/
		LocalPartMatrix GenerateLocalPartMatrix(const Element& element);
	}
};