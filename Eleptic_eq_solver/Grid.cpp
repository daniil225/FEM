#include "Grid.h"
#include <fstream>
#include <iostream>
#include <exception>
#include <cmath>
#include <iomanip>
#include <algorithm> // !!! Из библиотеки тяну лишь одну функцию не совсем хорошее решение для того что бы тянуть целую библиотеку 

namespace Grid
{

	/* Public */

	void Load(std::string filename, BaseGrid& baseGrid)
	{

		try
		{
			std::ifstream in;
			in.open(filename);
			if (!in.is_open())
			{
				std::cout << "File: " << filename << " can not open\n";
				exit(-1);
			}

			/* Базовая сетка по XZ */
			in >> baseGrid.Nx >> baseGrid.Nz;
			baseGrid.BaseGridXZ = new BaseGrid::PointXZS * [baseGrid.Nz];
			for (int i = 0; i < baseGrid.Nz; i++)
				baseGrid.BaseGridXZ[i] = new BaseGrid::PointXZS[baseGrid.Nx];

			for (int i = 0; i < baseGrid.Nz; i++)
			{
				for (int j = 0; j < baseGrid.Nx; j++)
				{
					in >> baseGrid.BaseGridXZ[i][j].x >> baseGrid.BaseGridXZ[i][j].z;
				}
			}
			/***********************/

			/* Базовая сетка по Y */
			in >> baseGrid.Ny;
			baseGrid.BaseGridY = new double[baseGrid.Ny];
			for (int i = 0; i < baseGrid.Ny; i++)
				in >> baseGrid.BaseGridY[i];
			/***********************/

			/* Расчетные подобласти */
			in >> baseGrid.L;
			baseGrid.CalculationArea = new int* [baseGrid.L];
			for (int i = 0; i < baseGrid.L; i++)
				baseGrid.CalculationArea[i] = new int[baseGrid.SizeOfCalculationAreaElemet];

			for (int i = 0; i < baseGrid.L; i++)
			{
				for (int j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
				{
					in >> baseGrid.CalculationArea[i][j];
					baseGrid.CalculationArea[i][j]--; // приведение нумераии с нуля 
				}
			}
			/***********************/

			/* Описание Границ */
			in >> baseGrid.P;
			baseGrid.BoundsArea = new int* [baseGrid.P];
			for (int i = 0; i < baseGrid.P; i++)
				baseGrid.BoundsArea[i] = new int[baseGrid.SizeOfBoundsAreaElement];

			for (int i = 0; i < baseGrid.P; i++)
			{
				for (int j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
				{
					in >> baseGrid.BoundsArea[i][j];
					baseGrid.BoundsArea[i][j]--;  // приведение нумераии с нуля
				}
			}
			/***********************/

			/* Правила дробления базовой сетки */
			baseGrid.DivideParam = new BaseGrid::DivideParamS * [baseGrid.SizeOfDivideParam];
			baseGrid.DivideParam[0] = new BaseGrid::DivideParamS[baseGrid.Nx - 1];
			baseGrid.DivideParam[1] = new BaseGrid::DivideParamS[baseGrid.Nz - 1];
			baseGrid.DivideParam[2] = new BaseGrid::DivideParamS[baseGrid.Ny - 1];

			for (int i = 0; i < baseGrid.Nx - 1; i++)
				in >> baseGrid.DivideParam[0][i].num >> baseGrid.DivideParam[0][i].coef;


			for (int i = 0; i < baseGrid.Nz - 1; i++)
				in >> baseGrid.DivideParam[1][i].num >> baseGrid.DivideParam[1][i].coef;


			for (int i = 0; i < baseGrid.Ny - 1; i++)
				in >> baseGrid.DivideParam[2][i].num >> baseGrid.DivideParam[2][i].coef;
			/***********************/

			in.close();
		}
		catch (std::bad_alloc& e)
		{
			std::cout << e.what();
			exit(-1);
		}
		catch (...)
		{
			std::cout << "Unknown error\n";
			exit(-1);
		}
	}

	void DivideGrid(BaseGrid& baseGrid, int coef)
	{
		baseGrid.CountOfDivision *= coef;

		for (int i = 0; i < baseGrid.Nx - 1; i++)
		{
			baseGrid.DivideParam[0][i].num *= coef;
			std::cout << baseGrid.DivideParam[0][i].num << "\n";
		}
		for (int i = 0; i < baseGrid.Nz - 1; i++)
		{
			baseGrid.DivideParam[1][i].num *= coef;
		}
		for (int i = 0; i < baseGrid.Ny - 1; i++)
		{
			baseGrid.DivideParam[2][i].num *= coef;
		}
	}

	void GenerateGrid(BaseGrid& baseGrid, GridS& grid)
	{
		// Расчет общего количества узлов в сетке 
		// Для этого пройдемся по массиву разбиения каждого отрезка и вычислим общее число узлов 
		Grid_private::GetTotalNumberOfNodes(baseGrid, grid);

		// Генерация базовой сетки 
		Grid_private::GenerateBaseGrid(baseGrid, grid);

		// Учет фиктивных узлов 
		Grid_private::DivisionIntoSubAreas(baseGrid, grid);
		// Учет КУ и расстановка границ 
		Grid_private::DivisionIntoSubBounds(baseGrid, grid);

	}

	void ReGenerateGrid(BaseGrid& baseGrid, GridS& grid)
	{
		DeallocateGrid(grid);
		grid.Dim = 0;
		grid.Nx = 0;
		grid.Ny = 0;
		grid.Nz = 0;
		GenerateGrid(baseGrid, grid);
	}

	Element GetFinitElement(const GridS& grid, int idx)
	{

		auto ZeroInitArray = [](int* arr, const int size)
		{
			for (int i = 0; i < size; i++)
				arr[i] = 0;
		};

		auto ZeroInitArrayPair = [](std::pair<int, int>* arr, const int size)
		{
			for (int i = 0; i < size; i++)
			{
				arr[i].first = 0;
				arr[i].second = 0;
			}
		};

		/* Возвращает индекс массива где лежит max значение  */
		auto GetArgMax = [](int* arr, const int size) -> int
		{
			int res = -1;
			for (int i = 0; i < size; i++)
			{
				if (arr[i] == 8)
				{
					res = i;
					break;
				}
			}

			return res;
		};

		auto GetArgMaxPair = [](std::pair<int, int>* arr, const int size) -> std::pair<int, int>
		{
			std::pair<int, int> res = std::make_pair(-1, -1);
			for (int i = 0; i < size; i++)
			{
				if (arr[i].first == 4)
				{
					res.first = i;
					res.second = arr[i].second;
					break;
				}
			}

			return res;
		};

		/* Вернет по номеру конечного элемента значение его праой нижней границе (локальный номер 1)  */
		int Nx = grid.Nx;
		int Ny = grid.Ny;
		int Nz = grid.Nz;
		auto K = [Nx, Ny, Nz](int idx) -> int
		{
			int shiftXZ = 0;
			int projidx = idx % (((Nx - 1) * (Nz - 1)));
			if (projidx < Nx - 1)
			{
				shiftXZ = projidx;
			}
			else
			{
				int level = std::floor(projidx / (Nx - 1));
				int start = level * Nx;
				int shift = projidx - (Nx - 1) * level;
				shiftXZ = start + shift;
			}

			int levelY = std::floor((idx) / ((Nx - 1) * (Nz - 1)));
			int shiftY = levelY * Nx * Nz;
			return shiftY + shiftXZ;
		};

		Element res;
		/* Расчет глобальных индексов */
		res.GlobalIdx[0] = K(idx);
		res.GlobalIdx[1] = res.GlobalIdx[0] + 1;
		res.GlobalIdx[2] = res.GlobalIdx[0] + Nx;
		res.GlobalIdx[3] = res.GlobalIdx[1] + Nx;

		res.GlobalIdx[4] = res.GlobalIdx[0] + Nx * Nz;
		res.GlobalIdx[5] = res.GlobalIdx[4] + 1;
		res.GlobalIdx[6] = res.GlobalIdx[4] + Nx;
		res.GlobalIdx[7] = res.GlobalIdx[5] + Nx;

		for (int i = 0; i < res.FinitElementSize; i++)
		{
			res.e[i] = grid.Grid[res.GlobalIdx[i]];

			/* Проверим сразу на фиктивность  */
			if (!PointInfo::IsFiFictitious(res.e[i].info))
			{
				/* Элемент фиктивный устанавливаем этот факт и прекращаем цикл */
				res.isFictive = false; // Да фиктивный 
				return res; // Возврат из функции
			}
		}
		/* Опредеелим теперь всю информацию о конечном элементе */
		/* Если мы здесь значит элемент не фиктивный */
		res.isFictive = true; // НЕ ФИКТИВНЫЙ 

		/* Идем по подобласти и определяем принадлежность к области */

		int Map[PointInfo::Comand::RightBoundValAreaInfo]; // Массив всевозможных элементов подобласти 
		ZeroInitArray(Map, PointInfo::Comand::RightBoundValAreaInfo);

		for (int i = 0; i < res.FinitElementSize; i++)
		{
			PointInfo::AreaInfo Area = PointInfo::GetAreaInfo(res.e[i].info);
			/* Идем по структурке и расставляем все типы что есть  */
			for (int j = 0; j < Area.size; j++)
				Map[Area.Cond[j]]++;
		}
		/* В карте лежит информация о количестве всевозможных типах */
		/* Нам нужна тот индекс где лежит 8 это будет max*/
		res.AreaInfo = GetArgMax(Map, PointInfo::Comand::RightBoundValAreaInfo);

		/* Разбираемся с краевыми */
		/* Установим локальные номера */
		// {1, 2, 3, 4}, { 5,6,7,8 }, { 1,3,5,7 }, { 2,4,6,8 }, { 1,2,5,6 }, { 3,4,7,8 }
		res.Bound[0].LocalIdx[0] = 0; res.Bound[0].LocalIdx[1] = 1; res.Bound[0].LocalIdx[2] = 2; res.Bound[0].LocalIdx[3] = 3;
		res.Bound[1].LocalIdx[0] = 4; res.Bound[1].LocalIdx[1] = 5; res.Bound[1].LocalIdx[2] = 6; res.Bound[1].LocalIdx[3] = 7;
		res.Bound[2].LocalIdx[0] = 0; res.Bound[2].LocalIdx[1] = 2; res.Bound[2].LocalIdx[2] = 4; res.Bound[2].LocalIdx[3] = 6;
		res.Bound[3].LocalIdx[0] = 1; res.Bound[3].LocalIdx[1] = 3; res.Bound[3].LocalIdx[2] = 5; res.Bound[3].LocalIdx[3] = 7;
		res.Bound[4].LocalIdx[0] = 0; res.Bound[4].LocalIdx[1] = 1; res.Bound[4].LocalIdx[2] = 4; res.Bound[4].LocalIdx[3] = 5;
		res.Bound[5].LocalIdx[0] = 2; res.Bound[5].LocalIdx[1] = 3; res.Bound[5].LocalIdx[2] = 6; res.Bound[5].LocalIdx[3] = 7;


		/* Даем глобальные */
		for (int i = 0; i < res.BoundCount; i++)
		{
			for (int j = 0; j < PointInfo::Comand::RightBoundBoundCount; j++)
			{
				res.Bound[i].GlobalIdx[j] = res.GlobalIdx[res.Bound[i].LocalIdx[j]];
			}
		}

		/* По все границам идем и устанавливаем информацию о границе */
		for (int i = 0; i < res.BoundCount; i++)
		{
			std::pair<int, int> MapB[PointInfo::Comand::RightBoundValBoundInfo];
			ZeroInitArrayPair(MapB, PointInfo::Comand::RightBoundValBoundInfo);

			for (int j = 0; j < PointInfo::Comand::RightBoundBoundCount; j++)
			{
				PointInfo::BoundInfo bnd = PointInfo::GetBoundInfo(res.e[res.Bound[i].LocalIdx[j]].info); // Получить информацию о границах на границе 

				/* Идем по структуре */
				for (int k = 0; k < bnd.size; k++)
				{
					MapB[bnd.Cond[k]].first++;
					MapB[bnd.Cond[k]].second = bnd.TypeCond[k];
				}
			}
			/* Получаем максимальный аргумент и устанавливаем значения для границы */
			std::pair<int, int> Cond = GetArgMaxPair(MapB, PointInfo::Comand::RightBoundValBoundInfo);
			res.Bound[i].BoundInfo = Cond.first;
			res.Bound[i].BoundType = Cond.second;
		}

		/* Расстановка фиктивных границ */
		for (int i = 0; i < res.BoundCount; i++)
		{
			if (res.Bound[i].BoundInfo != -1)
				res.Bound[i].IsBound = true;
		}

		return res;
	}

	void DeallocateGrid(GridS& grid)
	{
		if (grid.Grid != nullptr)
			delete[] grid.Grid;
	}

	void DeallocateBaseGrid(BaseGrid& baseGrid)
	{
		if (baseGrid.DivideParam != nullptr)
		{
			for (int i = 0; i < baseGrid.SizeOfDivideParam; i++)
			{
				if (baseGrid.DivideParam[i] != nullptr)
					delete[] baseGrid.DivideParam[i];
			}
			delete[] baseGrid.DivideParam;
		}

		if (baseGrid.BaseGridXZ != nullptr)
		{
			for (int i = 0; i < baseGrid.Nz; i++)
			{
				if (baseGrid.BaseGridXZ[i] != nullptr)
					delete[] baseGrid.BaseGridXZ[i];
			}

			delete[] baseGrid.BaseGridXZ;
		}

		if (baseGrid.BaseGridY != nullptr)
			delete[] baseGrid.BaseGridY;

		if (baseGrid.BoundsArea != nullptr)
		{
			for (int i = 0; i < baseGrid.P; i++)
			{
				if (baseGrid.BoundsArea[i] != nullptr)
					delete[] baseGrid.BoundsArea[i];
			}

			delete[] baseGrid.BoundsArea;
		}

		if (baseGrid.CalculationArea != nullptr)
		{
			for (int i = 0; i < baseGrid.L; i++)
			{
				if (baseGrid.CalculationArea[i] != nullptr)
					delete[] baseGrid.CalculationArea[i];
			}

			delete[] baseGrid.CalculationArea;
		}
	}

	/* Debug functions */
	void PrintElement(const Element& element)
	{
		std::cout << "Finit Element: \n";
		for (int i = 0; i < 8; i++)
			std::cout << "Global idx = " << element.GlobalIdx[i] << " Cord = (" << element.e[i].x << ";" << element.e[i].z << ";" << element.e[i].y << ")\n";
		std::cout << "\n";

		std::cout << "Info: \n";
		std::cout << "Is Fictitios: " << element.isFictive << "\n";
		std::cout << "Type Area: " << element.AreaInfo << "\n";

		std::cout << "\n";
		std::cout << "Info About Bound: \n";
		for (int i = 0; i < element.BoundCount; i++)
		{
			std::cout << "Is Bound: " << element.Bound[i].IsBound << "\n";
			std::cout << "Bound Formula: " << element.Bound[i].BoundInfo << "  ";
			std::cout << "Bound Type: " << element.Bound[i].BoundType << "\n";
			std::cout << "Local Idx: {";
			for (int j = 0; j < 4; j++)
			{
				std::cout << element.Bound[i].LocalIdx[j] << ";";
			}
			std::cout << "} \n";

			std::cout << "Global Idx: {";
			for (int j = 0; j < 4; j++)
			{
				std::cout << element.Bound[i].GlobalIdx[j] << ";";
			}
			std::cout << "} \n\n";
		}
	}

	void PrintListGrid(const GridS& grid, int list)
	{
		int Nx = grid.Nx;
		int Nz = grid.Nz;

		int idx = list * Nx * Nz;
		std::cout << "Start idx: " << idx << "  End idx: " << idx + Nx * Nz - 1 << " Step Row: " << Nx << "\n";

		for (int i = 0; i < Nz; i++)
		{
			for (int j = 0; j < Nx; j++)
			{
				std::cout << std::fixed << std::setprecision(2) << "(" << grid.Grid[idx].x << ";" << grid.Grid[idx].z << ";" << grid.Grid[idx].y << ") ";
				idx++;
			}
			std::cout << "\n";
		}
	}

	/* Private  */
	namespace Grid_private
	{

		void GetTotalNumberOfNodes(const BaseGrid& baseGrid, GridS& grid)
		{
			for (int i = 0; i < baseGrid.Nx - 1; i++)
				grid.Nx += baseGrid.DivideParam[0][i].num;

			for (int i = 0; i < baseGrid.Nz - 1; i++)
				grid.Nz += baseGrid.DivideParam[1][i].num;

			for (int i = 0; i < baseGrid.Ny - 1; i++)
				grid.Ny += baseGrid.DivideParam[2][i].num;

			grid.Nx++;
			grid.Ny++;
			grid.Nz++;
		}

		void GenerateBaseGrid(const BaseGrid& baseGrid, GridS& grid)
		{
			struct SettingForDivide
			{
				double step; // Шаг на отрезке
				double coef; // Коэффициент увеличения шага 
				int num; // Количество интервалов идем то num-1 и потом явно вставляем элемент 
			};

			/* Расчитываем шаг для сетки  */
			/*
				@param
				int i - Номер массива от 0 до 2
				int j - Номер элемента в массиве
				double left - левая грани отрезка
				double right - правая граница отрезка
				ret: SettingForDivide -  структура с вычесленными параметрами деления сетки
			*/
			auto CalcSettingForDivide = [&](int i, int j, double left, double right) -> SettingForDivide
			{
				SettingForDivide res;
				int num = baseGrid.DivideParam[i][j].num;
				double coef = baseGrid.DivideParam[i][j].coef;

				if (coef > 1.0)
				{
					double coefStep = 1.0 + (coef * (std::pow(coef, num - 1) - 1.0)) / (coef - 1.0);

					res.step = (right - left) / coefStep;
				}
				else
				{
					res.step = (right - left) / num;
				}

				//  Убираем погрешность 
				if (std::abs(res.step) < eps)
					res.step = 0.0;

				res.num = num;
				res.coef = coef;
				return res;
			};

			/* Генерация разбиения по X или Z( гороизонтальная линия или вертикальная ) с учетом разбиения */
			/*
				@param
				SettingForDivide &param - параметр разбиения
				double left - левая граница отрезка
				double right - правая граница отрезка
				double *Line - генерируемый массив
				int &idx - индекс в массиве на какую позицию ставить элемент
			*/
			auto GenerateDivide = [](SettingForDivide& param, double left, double right, double* Line, int& idx) -> void
			{
				int num = param.num;
				double coef = param.coef;
				double step = param.step;

				Line[idx] = left;
				idx++;
				double ak = left;
				for (int k = 0; k < num - 1; k++)
				{
					ak = ak + step * std::pow(coef, k);
					Line[idx] = ak;
					idx++;
				}
				Line[idx] = right;
			};


			/* Вернет число соответсвующее стартовой позиции по сути это скачок  */
			/*
				@param
				int i - номер
				int axis - соответсвующая ость 0 - x, 1 - z 2 - y
				ret int: Величина скачка в сетке
			*/
			auto Getlevel = [&baseGrid](int i, int axis) -> int
			{
				int res = 0;
				for (int k = 0; k < i; k++)
					res += baseGrid.DivideParam[axis][k].num;
				return res;
			};

			try
			{
				int total = grid.Nx * grid.Ny * grid.Nz;
				grid.Dim = total; // Инициализация параметра размерности 
				grid.Grid = new Point[total];
				// Псевдоним для быстрого обращения 
				Point* Grid = grid.Grid;
				int Nx = baseGrid.Nx;
				int Nz = baseGrid.Nz;
				int Ny = baseGrid.Ny;
				BaseGrid::PointXZS** BaseGridXZ = baseGrid.BaseGridXZ;

				int GlobalNx = grid.Nx;
				int GlobalNz = grid.Nz;
				int GlobalNy = grid.Ny;
				double* LineX = new double[GlobalNx]; // Массив элементов в строке по Х
				double* LineZ = new double[GlobalNz]; // Массив элементов в строке по Z
				double* LineY = new double[GlobalNy]; // Массив элементов в строке по Y
				// Сгенерируем одну плоскость xz, а потом ее растиражируем по y

				/* Разбиение по x и z */
				/* Расстановка элементов основных линий с учетом их расположения в сетке ( Опорная сетка ) */

				/* Расстановка элементов по x ( Опорных ) */

				/* Нужно значения х расставить в соответсвующие строки они соответсвуют разбиению по z */

				for (int i = 0; i < Nz; i++)
				{
					int idx = 0;

					for (int j = 0; j < Nx - 1; j++)
					{
						double left = baseGrid.BaseGridXZ[i][j].x;
						double right = baseGrid.BaseGridXZ[i][j + 1].x;
						SettingForDivide param = CalcSettingForDivide(0, j, left, right);
						GenerateDivide(param, left, right, LineX, idx);
					}
					/* Заносим соответствующие значения x на свои позиции */

					int startIdx = Getlevel(i, 1) * GlobalNx;
					int endIdx = startIdx + GlobalNx;
					for (int k = startIdx, kk = 0; k < endIdx; k++, kk++)
						Grid[k].x = LineX[kk];
				}


				/* Расстановка элементов по z ( Опорных ) */
				for (int i = 0; i < Nx; i++)
				{
					int idx = 0;
					for (int j = 0; j < Nz - 1; j++)
					{
						double left = baseGrid.BaseGridXZ[j][i].z;
						double right = baseGrid.BaseGridXZ[j + 1][i].z;
						SettingForDivide param = CalcSettingForDivide(1, j, left, right);
						GenerateDivide(param, left, right, LineZ, idx);
					}

					/* Процедура расстановки узлов в глобальный массив */
					int startIdx = Getlevel(i, 0); // Стартовый индекс для прохода по массиву 
					for (int k = 0; k < GlobalNz; k++)
					{
						// Скачки будут ровно на величину GlobalNx
						Grid[startIdx].z = LineZ[k];
						startIdx += GlobalNx;
					}
				}


				/****************************************************/

				/* Генерация вспомогательных вертикальных линий  */
				/*
					Кратко Алгоритм:
					Работаем с осью Z соответственно индексация будет происходить по этой оси
					в цикле идем по всем столбцам массива сетки

					Нужно получить левую и правую границу на каждом интервале
					Сформировать массив отрезков по  данной координате
					Занести полученный массив в Глобальную сетку
				*/
				/* Цикл по всем горизонтальным линиям */
				for (int i = 0; i < GlobalNx; i++)
				{
					int idx = 0;
					/* Цикл по интеравалам оси Z */
					for (int j = 0; j < Nz - 1; j++)
					{
						int startIdx = i + GlobalNx * Getlevel(j, 1);
						int endIdx = i + GlobalNx * Getlevel(j + 1, 1);
						double left = Grid[startIdx].x; // Левая граница по х
						double right = Grid[endIdx].x; // Правая граница по х 

						// Разбиение интервала подчиняется разбиению по оси z
						SettingForDivide param = CalcSettingForDivide(1, j, left, right);
						GenerateDivide(param, left, right, LineZ, idx);
					}
					/* Занесение результата в Итоговый массив */

					int startIdx = i; // Стартовая позиция
					for (int k = 0; k < GlobalNz; k++)
					{
						Grid[startIdx].x = LineZ[k];
						startIdx += GlobalNx;
					}
				}
				/****************************************************/

				/* Генерация вспомогательных горизонтальных линий  */
				/* Цикл по всем горизонтальным линиям */
				for (int i = 0; i < GlobalNz; i++)
				{
					int idx = 0;
					for (int j = 0; j < Nx - 1; j++)
					{
						int startIdx = Getlevel(j, 0) + i * GlobalNx;
						int endIdx = Getlevel(j + 1, 0) + i * GlobalNx;
						double left = Grid[startIdx].z;
						double right = Grid[endIdx].z;
						// Разбиение интервала подчиняется разбиению по оси x
						SettingForDivide param = CalcSettingForDivide(0, j, left, right);
						GenerateDivide(param, left, right, LineX, idx);

					}
					/* Занесение результатов в Глобальную сетку */
					int  startIdx = i * GlobalNx;
					for (int k = 0; k < GlobalNx; k++)
					{
						Grid[startIdx].z = LineX[k];
						startIdx++;
					}
				}
				/****************************************************/



				/* Разбиение по y теражируем сечение с учетом сечения плоскостью */
				// Создадим разбиение элементов массива Y
				int idxY = 0;
				for (int i = 0; i < Ny - 1; i++)
				{
					double left = baseGrid.BaseGridY[i];
					double right = baseGrid.BaseGridY[i + 1];
					SettingForDivide param = CalcSettingForDivide(2, i, left, right);
					GenerateDivide(param, left, right, LineY, idxY);
				}

				// Вставляем элемент Y0 в базовом сечении 
				int idx = 0;
				for (int j = 0; j < GlobalNz; j++)
				{
					for (int k = 0; k < GlobalNx; k++)
					{
						PointInfo::ClearInfo(Grid[idx].info); // Инициализация битовых полей без этой команды будет Ошибка !!!
						Grid[idx].y = LineY[0];
						idx++;
					}
				}

				// Тиражирование сетки по всем узлам 
				// i отвечает за индексацию в массиве по Y
				// Idx за индексацию по элементам 


				for (int i = 1; i < GlobalNy; i++)
				{
					int ListIdx = 0;
					for (int j = 0; j < GlobalNz; j++)
					{
						for (int k = 0; k < GlobalNx; k++)
						{
							PointInfo::ClearInfo(Grid[idx].info); // Инициализация битовых полей без этой команды будет Ошибка !!!
							Grid[idx].y = LineY[i];
							Grid[idx].x = Grid[ListIdx].x;
							Grid[idx].z = Grid[ListIdx].z;
							ListIdx++;
							idx++;
						}
					}
				}

				/* Очистка памяти */
				delete[] LineX;
				delete[] LineZ;
				delete[] LineY;

			}
			catch (std::bad_alloc& e)
			{
				std::cout << e.what();
				exit(-1);
			}
			catch (...)
			{
				std::cout << "Unknown error\n";
				exit(-1);
			}
		}

		void DivisionIntoSubAreas(const BaseGrid& baseGrid, GridS& grid)
		{

			struct BoundArea
			{
				int AreaNum = -1; // номер подобласти (определяет набор формул отвечающих параметрам ДУ)  
				int PlaneXZSize = 0; // Размер массива PlaneXZ
				int* PlaneXZ; // Массив точек многоугольника (номера)

				double StartY = 0;
				double EndY = 0;
			};

			auto Getlevel = [&baseGrid](int i, int axis) -> int
			{
				int res = 0;
				for (int k = 0; k < i; k++)
					res += baseGrid.DivideParam[axis][k].num;
				return res;
			};

			try
			{
				/* Размер массива берем с запасом его размер равен Nx*Nz
				   контур номеров точек многоугольника соттветствующих каким то координатам
				*/
				int* BoundGrid = new int[grid.Nx * grid.Nz];
				/* Получает подобласть в соответсвии с ее номером */
				/*
					@param
					int i - Номер подобласти (Порядковый) в массиве определяется порядок
					ret BoundArea - сформированный массив подобласти
				*/
				auto GetBound = [&](int i) -> BoundArea
				{
					int GlobalNx = grid.Nx;
					int GlobalNz = grid.Nz;
					BoundArea Bound;

					Bound.AreaNum = baseGrid.CalculationArea[i][0]; // Выставили номер подобласти 

					/* Вычислим все номера принадлежащие данной области */
					/* Инедексы границ многоугольника в глобальной нумерации */
					/* XZ */

					int leftStartX = Getlevel(baseGrid.CalculationArea[i][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][3], 1);
					int leftEndX = Getlevel(baseGrid.CalculationArea[i][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][4], 1);
					int rightStartX = Getlevel(baseGrid.CalculationArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][3], 1);
					int rightEndX = Getlevel(baseGrid.CalculationArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[i][4], 1);

					int nX = (rightEndX - rightStartX) / GlobalNx;
					int nZ = (rightStartX - leftStartX);

					int idxBoundGrid = 0;
					for (int i = 0; i <= nX; i++)
					{
						int Idx = leftStartX + i * GlobalNx;
						for (int j = 0; j <= nZ; j++)
						{
							BoundGrid[idxBoundGrid] = Idx;
							idxBoundGrid++;
							Idx++;
						}
					}

					/* Y */
					/* В силу линейности по оси можно взять только первое и последнее значение */
					Bound.StartY = baseGrid.BaseGridY[baseGrid.CalculationArea[i][5]];
					Bound.EndY = baseGrid.BaseGridY[baseGrid.CalculationArea[i][6]];
					Bound.PlaneXZSize = idxBoundGrid;
					Bound.PlaneXZ = BoundGrid;

					return Bound;
				};

				/* Бинарный поиск по массиву */
				/*
					@param
					const BoundArea& Bound - Массив границ
					int numPointGlobal - Значение для поиска
					ret bool: true - элемент в массиве есть false в противном случае

				*/
				auto BinarySerch = [](const BoundArea& Bound, int numPointGlobal) -> bool
				{
					int* arr = Bound.PlaneXZ;
					int left = 0;
					int right = Bound.PlaneXZSize - 1;
					int midd = 0;
					while (1)
					{
						midd = (left + right) / 2;

						if (numPointGlobal < arr[midd])       // если искомое меньше значения в ячейке
							right = midd - 1;      // смещаем правую границу поиска
						else if (numPointGlobal > arr[midd])  // если искомое больше значения в ячейке
							left = midd + 1;    // смещаем левую границу поиска
						else                       // иначе (значения равны)
							break;           // функция возвращает индекс ячейки

						// если границы сомкнулись
						if (left > right)
						{
							midd = -1;
							break;
						}
					}
					if (midd == -1)
						return false;
					else
						return true;
				};

				/* Проверяет принадлежит ли точка Заданной области */
				/*
					@param
					const BoundArea &Bound - Заданная граница
					Point &point - Точка для проверки
					int numPointGlobal - индекс точки в глобальной нумерации
					ret bool: true - Точка принадлежит заданной области, false в противном случае
				*/
				auto IsInArea = [&](const BoundArea& Bound, Point& point, int numPointGlobal) -> bool
				{
					auto le = [](double x1, double x2) -> bool
					{
						if (x1 < x2) return true;
						if (std::abs(x1 - x2) < eps) return true;
						return false;
					};

					auto ge = [&](double x1, double x2) -> bool
					{
						return le(x2, x1);
					};

					numPointGlobal %= grid.Nx * grid.Nz; // Проекция точки на плоскость 
					bool arg1 = BinarySerch(Bound, numPointGlobal);
					bool arg2 = false;
					if (le(point.z, Bound.EndY) && ge(point.z, Bound.StartY)) arg2 = true;

					return arg1 && arg2;

				};


				/* В цикле по всем областям */
				int N = grid.Dim;
				for (int i = 0; i < baseGrid.L; i++)
				{
					BoundArea Bound = GetBound(i);

					/* По той части области где распологаются элементы (включая фиктивный)*/
					/* Для этого получим минимальный и максимальный значений индексов */
					// int leftAreaBound = ;
					// int RightAreaBound = ;

					/* Этот цикл тупой он по всем элементам идет */
					for (int j = 0; j < N; j++)
					{
						if (IsInArea(Bound, grid.Grid[j], j))
						{
							/* Мы в заданной области устанавливаем нужные параметры */
							PointInfo::SetFictitiousBit(grid.Grid[j].info, 1);
							PointInfo::SetAreaInfo(grid.Grid[j].info, Bound.AreaNum);
						}
					}
				}

				/* Очистка локальной переменной  */
				delete[] BoundGrid;

			}
			catch (std::bad_alloc& e)
			{
				std::cout << e.what();
				exit(-1);
			}
			catch (...)
			{
				std::cout << "Unknown error\n";
				exit(-1);
			}


			//В цикле по всем элементам будем определять принадлежность какой-либо расчетной области 

		}

		void DivisionIntoSubBounds(const BaseGrid& baseGrid, GridS& grid)
		{

			struct Bound
			{
				int Size = 0;
				int* Plane = nullptr; // Массив точек границы 
				int BoundType = 0;
				int BoundFormula = -1;
			};

			auto Getlevel = [&baseGrid](int i, int axis) -> int
			{
				int res = 0;
				for (int k = 0; k < i; k++)
					res += baseGrid.DivideParam[axis][k].num;
				return res;
			};

			auto GetMaxPair = [](int a, int b, int c) -> std::pair<int, int>
			{

				int arr[3]{ a,b,c };
				std::sort(arr, arr + 3);
				//std::cout << arr[0] << " " << arr[1] << " " << arr[2] << "\n";

				return std::make_pair(arr[1], arr[2]);
			};

			/* Бинарный поиск по массиву */
				/*
					@param
					const BoundArea& Bound - Массив границ
					int numPointGlobal - Значение для поиска
					ret bool: true - элемент в массиве есть false в противном случае

				*/
			auto BinarySerch = [](const Bound& bound, int numPointGlobal) -> bool
			{
				int* arr = bound.Plane;
				int left = 0;
				int right = bound.Size - 1;
				int midd = 0;
				while (1)
				{
					midd = (left + right) / 2;

					if (numPointGlobal < arr[midd])       // если искомое меньше значения в ячейке
						right = midd - 1;      // смещаем правую границу поиска
					else if (numPointGlobal > arr[midd])  // если искомое больше значения в ячейке
						left = midd + 1;    // смещаем левую границу поиска
					else                       // иначе (значения равны)
						break;           // функция возвращает индекс ячейки

					// если границы сомкнулись
					if (left > right)
					{
						midd = -1;
						break;
					}
				}
				if (midd == -1)
					return false;
				else
					return true;
			};

			auto IsInBound = [&](const Bound& bound, int numPointGlobal) { return BinarySerch(bound, numPointGlobal); };

			try
			{
				std::pair<int, int> max2 = GetMaxPair(grid.Nx, grid.Nz, grid.Ny);
				int* MemoryPool = new int[max2.first * max2.second];
				/* Краевые условия задаются практически так же как и области с тем лишь исключением, что одна из координат фиксируется */
				auto GetBound = [&](int i) -> Bound
				{
					int GlobalNx = grid.Nx;
					int GlobalNz = grid.Nz;
					int GlobalNy = grid.Ny;

					Bound bound;
					bound.BoundType = baseGrid.BoundsArea[i][1];
					bound.BoundFormula = baseGrid.BoundsArea[i][0];

					/* Это определяется равенством  соответствующих точек*/
					/* y - фиксирован */
					if (baseGrid.BoundsArea[i][6] == baseGrid.BoundsArea[i][7])
					{
						/*  Определим смещение для оси y */
						int startBaseIdx = GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // Узнаем край для 

						/* Определяем границы областей */
						int leftStartX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1);
						int leftEndX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1);
						int rightStartX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1);
						int rightEndX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1);

						/* Сколько шагать  */
						int nX = (rightEndX - rightStartX) / GlobalNx;
						int nZ = (rightStartX - leftStartX);

						int idxBoundGrid = 0;
						for (int i = 0; i <= nX; i++)
						{
							int Idx = leftStartX + i * GlobalNx + startBaseIdx; // Смещенный индекс 
							for (int j = 0; j <= nZ; j++)
							{
								MemoryPool[idxBoundGrid] = Idx;
								idxBoundGrid++;
								Idx++;
							}
						}

						bound.Plane = MemoryPool;
						bound.Size = idxBoundGrid;
					}
					/* х - фиксирован */
					else if (baseGrid.BoundsArea[i][2] == baseGrid.BoundsArea[i][3])
					{
						/* Определим базовые узлы на прямой OZ а потом растиражируем узлы по границе ( Дробление фиксированное )  шаг по массиву GlobalNx*/
						int StartPozitionZ = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // Стартовая позиция для  Z
						int EndPozitionZ = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1); // Конечная точка для оси Z

						int StartPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // Стартовая позиция для Y
						int EndPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][7], 2); // Конечня позиция для Y

						int nZ = (EndPozitionZ - StartPozitionZ) / GlobalNx; //  Количество узлов по оси Z
						int nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // Количество шагов для Y

						int idxBoundGrid = 0;
						for (int i = 0; i <= nY; i++)
						{
							int Idx = StartPositionY + i * GlobalNx * GlobalNz;
							for (int j = 0; j <= nZ; j++)
							{
								MemoryPool[idxBoundGrid] = Idx;
								idxBoundGrid++;
								Idx += GlobalNx;
							}
						}
						bound.Plane = MemoryPool;
						bound.Size = idxBoundGrid;
					}
					/* z - фиксирован */
					else if (baseGrid.BoundsArea[i][4] == baseGrid.BoundsArea[i][5])
					{
						/* Определяем базовые узлы по оси OX, а потом аналогично растиражируем узлы по границе (Дробление фиксированное) шаг по массиву + 1 */
						int StartPositionX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // Стартовая позиция по оси Х
						int EndPositionX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // Конечная позиция по оси Х

						int StartPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // Стартовая позиция для Y
						int EndPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][7], 2); // Конечня позиция для Y правй конец (если смотреть со стороны положительного направления) 

						int nX = EndPositionX - StartPositionX; //  Количество узлов по оси X
						int nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // Количество шагов для Y 

						int idxBoundGrid = 0;
						for (int i = 0; i <= nY; i++)
						{
							int Idx = StartPositionY + StartPositionX + i * GlobalNx * GlobalNz;
							for (int j = 0; j <= nX; j++)
							{
								MemoryPool[idxBoundGrid] = Idx;
								idxBoundGrid++;
								Idx++;
							}
						}
						bound.Plane = MemoryPool;
						bound.Size = idxBoundGrid;
					}

					return bound;
				};


				/* Цикл по всем границам */
				int N = grid.Dim;
				for (int i = 0; i < baseGrid.P; i++)
				{
					/* Получить список точек границы  */
					Bound bound = GetBound(i);
					for (int j = 0; j < N; j++)
					{
						/* Расставляем нужные значения в соответствующие точки сетки */
						if (IsInBound(bound, j))
						{
							/* Мы в заданной области устанавливаем нужные параметры */
							PointInfo::SetBoundInfo(grid.Grid[j].info, bound.BoundFormula, bound.BoundType + 1);
						}

					}
				}

				delete[] MemoryPool;
			}
			catch (std::bad_alloc& e)
			{
				std::cout << e.what();
				exit(-1);
			}
			catch (...)
			{
				std::cout << "Unknown error\n";
				exit(-1);
			}
		}
	};

	/* Секция с функциями и структурами определяющими дополнительную информацию о Точке в области  */
	namespace PointInfo
	{
		/* Функции управления битовыми структурами */
		/* Private Section */

		void ClearInfo(Info& info)
		{
			info.BaseInfo &= Comand::BaseInfoClear;
			info.AreaInfo &= Comand::AreaInfoClear;
			info.BoundInfo &= Comand::BoundInfoClear;
			info.TypeBoundCond &= Comand::TypeBoundCondClear;
		}

		void SetFictitiousBit(Info& info, uint8_t val)
		{
			if (val >= 0 && val <= 1)
			{
				if (IsFiFictitious(info)) info.BaseInfo &= Comand::SetZeroFiFictitious;
				info.BaseInfo |= val;
			}
			else
				throw "Error Bit Operation\n";
		}

		void SetAreaCountBits(Info& info, uint32_t val)
		{
			if (val >= 0 && val <= Comand::RightBoundAreaCount)
			{
				if (GetAreaCount(info) != 0) info.BaseInfo &= Comand::SetZeroAreaCount;
				info.BaseInfo |= val << Comand::ShiftAreaCountBits;
			}
			else
				throw "Error Bit Operation\n";
		}

		void SetBoundCountBits(Info& info, uint32_t val)
		{
			if (val >= 0 && val <= Comand::RightBoundBoundCount)
			{
				if (GetBoundCount(info) != 0) info.BaseInfo &= Comand::SetZeroBountCount;
				info.BaseInfo |= val << Comand::ShiftBoundCountBits;
			}
			else
				throw "Error Bit Operation\n";
		}

		uint32_t GetAreaCount(const Info& info)
		{
			return (info.BaseInfo & Comand::GetAreaCount) >> Comand::ShiftAreaCountBits;
		}


		uint32_t GetBoundCount(const Info& info)
		{
			return (info.BaseInfo & Comand::GetBoundCount) >> Comand::ShiftBoundCountBits;
		}
		/*********************************************/

		/* Public section */

		void SetBoundInfo(Info& info, uint32_t val, uint32_t boundType)
		{
			if (val >= 0 && val <= Comand::RightBoundValBoundInfo && boundType >= 0 && boundType <= Comand::RightBoundValBoundType)
			{
				int Bnum = GetBoundCount(info);
				/* если не было значений */
				if (Bnum == 0)
				{
					info.BoundInfo |= val;
					info.TypeBoundCond |= boundType;
					SetBoundCountBits(info, 1);
				}
				else
				{
					BoundInfo Bound = GetBoundInfo(info);
					for (int i = 0; i < Bnum; i++)
					{
						// Проверяем дубликаты 
						if ((Bound.Cond[i] ^ val) == 0) // Есть дубликат 
							return;
					}

					/* Новый элемент */
					info.BoundInfo |= val << Comand::BaseBoundShift * Bnum;
					info.TypeBoundCond |= boundType << Comand::BaseTypeBoundShift * Bnum;
					SetBoundCountBits(info, Bnum + 1);
				}
			}
			else
				throw "Error Bit Operation\n";

		}


		void SetAreaInfo(Info& info, uint32_t val)
		{
			if (val >= 0 && val <= Comand::RightBoundValAreaInfo)
			{

				int Anum = GetAreaCount(info);
				if (Anum == 0)
				{
					info.AreaInfo |= val;
					SetAreaCountBits(info, 1);
				}
				else
				{
					AreaInfo Area = GetAreaInfo(info);
					for (int i = 0; i < Anum; i++)
					{
						// Проверяем дубликаты 
						if ((Area.Cond[i] ^ val) == 0) // Есть дубликат 
							return;
					}

					/* Новый элемент */
					info.AreaInfo |= val << Comand::BaseAreaShift * Anum;
					SetAreaCountBits(info, Anum + 1);
				}
			}
			else
				throw "Error Bit Operation\n";
		}

		bool IsFiFictitious(const Info& info)
		{
			return info.BaseInfo & Comand::GetFiFictitious;
		}


		bool IsBound(const Info& info)
		{
			return GetBoundCount(info) == 0 ? false : true;
		}


		BoundInfo GetBoundInfo(const Info& info)
		{
			BoundInfo Bound;
			int Bnum = GetBoundCount(info);

			Bound.size = Bnum;

			for (int i = 0; i <= Bnum - 1; i++)
			{
				Bound.Cond[i] = (info.BoundInfo & (Comand::GetBoundInfoBaseComand << i * Comand::BaseBoundShift)) >> i * Comand::BaseBoundShift;
				Bound.TypeCond[i] = (info.TypeBoundCond & (Comand::GetTypeBoundBaseComand << i * Comand::BaseTypeBoundShift)) >> i * Comand::BaseTypeBoundShift;
			}

			return Bound;
		}

		AreaInfo GetAreaInfo(const Info& info)
		{
			AreaInfo Area;
			int Anum = GetAreaCount(info);
			Area.size = Anum;

			for (int i = 0; i <= Anum - 1; i++)
				Area.Cond[i] = (info.AreaInfo & (Comand::GetAreaInfoBaseComand << Comand::BaseAreaShift * i)) >> Comand::BaseAreaShift * i;

			return Area;
		}

		void PrintInfo(const Info& info)
		{
			std::cout << "Info Struct\n";
			std::cout << "IsFiFictitious: " << IsFiFictitious(info) << "\n";
			std::cout << "AreaCount: " << GetAreaCount(info) << "\n";
			std::cout << "BoundCount: " << GetBoundCount(info) << "\n";
			PrintAreaInfo(GetAreaInfo(info));
			PrintBoundInfo(GetBoundInfo(info));

		}

		void PrintBoundInfo(const BoundInfo& Bound)
		{
			std::cout << "\nBound Info\n";
			for (int i = 0; i < Bound.size; i++)
				std::cout << "K = " << i + 1 << " Num formula: " << (uint32_t)Bound.Cond[i] << " Type Bound: " << (uint32_t)Bound.TypeCond[i] << "\n";
		}

		void PrintAreaInfo(const AreaInfo& Area)
		{
			std::cout << "\nArea Info\n";
			for (int i = 0; i < Area.size; i++)
				std::cout << "K = " << i + 1 << " Num formula: " << (uint32_t)Area.Cond[i] << "\n";
		}
	};
};