#include "Grid.h"
#include <fstream>
#include <iostream>
#include <exception>
#include <cmath>
#include <iomanip>
#include <algorithm> // !!! �� ���������� ���� ���� ���� ������� �� ������ ������� ������� ��� ���� ��� �� ������ ����� ���������� 

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

			/* ������� ����� �� XZ */
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

			/* ������� ����� �� Y */
			in >> baseGrid.Ny;
			baseGrid.BaseGridY = new double[baseGrid.Ny];
			for (int i = 0; i < baseGrid.Ny; i++)
				in >> baseGrid.BaseGridY[i];
			/***********************/

			/* ��������� ���������� */
			in >> baseGrid.L;
			baseGrid.CalculationArea = new int* [baseGrid.L];
			for (int i = 0; i < baseGrid.L; i++)
				baseGrid.CalculationArea[i] = new int[baseGrid.SizeOfCalculationAreaElemet];

			for (int i = 0; i < baseGrid.L; i++)
			{
				for (int j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
				{
					in >> baseGrid.CalculationArea[i][j];
					baseGrid.CalculationArea[i][j]--; // ���������� �������� � ���� 
				}
			}
			/***********************/

			/* �������� ������ */
			in >> baseGrid.P;
			baseGrid.BoundsArea = new int* [baseGrid.P];
			for (int i = 0; i < baseGrid.P; i++)
				baseGrid.BoundsArea[i] = new int[baseGrid.SizeOfBoundsAreaElement];

			for (int i = 0; i < baseGrid.P; i++)
			{
				for (int j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
				{
					in >> baseGrid.BoundsArea[i][j];
					baseGrid.BoundsArea[i][j]--;  // ���������� �������� � ����
				}
			}
			/***********************/

			/* ������� ��������� ������� ����� */
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
		// ������ ������ ���������� ����� � ����� 
		// ��� ����� ��������� �� ������� ��������� ������� ������� � �������� ����� ����� ����� 
		Grid_private::GetTotalNumberOfNodes(baseGrid, grid);

		// ��������� ������� ����� 
		Grid_private::GenerateBaseGrid(baseGrid, grid);

		// ���� ��������� ����� 
		Grid_private::DivisionIntoSubAreas(baseGrid, grid);
		// ���� �� � ����������� ������ 
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

		/* ���������� ������ ������� ��� ����� max ��������  */
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

		/* ������ �� ������ ��������� �������� �������� ��� ����� ������ ������� (��������� ����� 1)  */
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
		/* ������ ���������� �������� */
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

			/* �������� ����� �� �����������  */
			if (!PointInfo::IsFiFictitious(res.e[i].info))
			{
				/* ������� ��������� ������������� ���� ���� � ���������� ���� */
				res.isFictive = false; // �� ��������� 
				return res; // ������� �� �������
			}
		}
		/* ���������� ������ ��� ���������� � �������� �������� */
		/* ���� �� ����� ������ ������� �� ��������� */
		res.isFictive = true; // �� ��������� 

		/* ���� �� ���������� � ���������� �������������� � ������� */

		int Map[PointInfo::Comand::RightBoundValAreaInfo]; // ������ ������������ ��������� ���������� 
		ZeroInitArray(Map, PointInfo::Comand::RightBoundValAreaInfo);

		for (int i = 0; i < res.FinitElementSize; i++)
		{
			PointInfo::AreaInfo Area = PointInfo::GetAreaInfo(res.e[i].info);
			/* ���� �� ���������� � ����������� ��� ���� ��� ����  */
			for (int j = 0; j < Area.size; j++)
				Map[Area.Cond[j]]++;
		}
		/* � ����� ����� ���������� � ���������� ������������ ����� */
		/* ��� ����� ��� ������ ��� ����� 8 ��� ����� max*/
		res.AreaInfo = GetArgMax(Map, PointInfo::Comand::RightBoundValAreaInfo);

		/* ����������� � �������� */
		/* ��������� ��������� ������ */
		// {1, 2, 3, 4}, { 5,6,7,8 }, { 1,3,5,7 }, { 2,4,6,8 }, { 1,2,5,6 }, { 3,4,7,8 }
		res.Bound[0].LocalIdx[0] = 0; res.Bound[0].LocalIdx[1] = 1; res.Bound[0].LocalIdx[2] = 2; res.Bound[0].LocalIdx[3] = 3;
		res.Bound[1].LocalIdx[0] = 4; res.Bound[1].LocalIdx[1] = 5; res.Bound[1].LocalIdx[2] = 6; res.Bound[1].LocalIdx[3] = 7;
		res.Bound[2].LocalIdx[0] = 0; res.Bound[2].LocalIdx[1] = 2; res.Bound[2].LocalIdx[2] = 4; res.Bound[2].LocalIdx[3] = 6;
		res.Bound[3].LocalIdx[0] = 1; res.Bound[3].LocalIdx[1] = 3; res.Bound[3].LocalIdx[2] = 5; res.Bound[3].LocalIdx[3] = 7;
		res.Bound[4].LocalIdx[0] = 0; res.Bound[4].LocalIdx[1] = 1; res.Bound[4].LocalIdx[2] = 4; res.Bound[4].LocalIdx[3] = 5;
		res.Bound[5].LocalIdx[0] = 2; res.Bound[5].LocalIdx[1] = 3; res.Bound[5].LocalIdx[2] = 6; res.Bound[5].LocalIdx[3] = 7;


		/* ���� ���������� */
		for (int i = 0; i < res.BoundCount; i++)
		{
			for (int j = 0; j < PointInfo::Comand::RightBoundBoundCount; j++)
			{
				res.Bound[i].GlobalIdx[j] = res.GlobalIdx[res.Bound[i].LocalIdx[j]];
			}
		}

		/* �� ��� �������� ���� � ������������� ���������� � ������� */
		for (int i = 0; i < res.BoundCount; i++)
		{
			std::pair<int, int> MapB[PointInfo::Comand::RightBoundValBoundInfo];
			ZeroInitArrayPair(MapB, PointInfo::Comand::RightBoundValBoundInfo);

			for (int j = 0; j < PointInfo::Comand::RightBoundBoundCount; j++)
			{
				PointInfo::BoundInfo bnd = PointInfo::GetBoundInfo(res.e[res.Bound[i].LocalIdx[j]].info); // �������� ���������� � �������� �� ������� 

				/* ���� �� ��������� */
				for (int k = 0; k < bnd.size; k++)
				{
					MapB[bnd.Cond[k]].first++;
					MapB[bnd.Cond[k]].second = bnd.TypeCond[k];
				}
			}
			/* �������� ������������ �������� � ������������� �������� ��� ������� */
			std::pair<int, int> Cond = GetArgMaxPair(MapB, PointInfo::Comand::RightBoundValBoundInfo);
			res.Bound[i].BoundInfo = Cond.first;
			res.Bound[i].BoundType = Cond.second;
		}

		/* ����������� ��������� ������ */
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
				double step; // ��� �� �������
				double coef; // ����������� ���������� ���� 
				int num; // ���������� ���������� ���� �� num-1 � ����� ���� ��������� ������� 
			};

			/* ����������� ��� ��� �����  */
			/*
				@param
				int i - ����� ������� �� 0 �� 2
				int j - ����� �������� � �������
				double left - ����� ����� �������
				double right - ������ ������� �������
				ret: SettingForDivide -  ��������� � ������������ ����������� ������� �����
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

				//  ������� ����������� 
				if (std::abs(res.step) < eps)
					res.step = 0.0;

				res.num = num;
				res.coef = coef;
				return res;
			};

			/* ��������� ��������� �� X ��� Z( ��������������� ����� ��� ������������ ) � ������ ��������� */
			/*
				@param
				SettingForDivide &param - �������� ���������
				double left - ����� ������� �������
				double right - ������ ������� �������
				double *Line - ������������ ������
				int &idx - ������ � ������� �� ����� ������� ������� �������
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


			/* ������ ����� �������������� ��������� ������� �� ���� ��� ������  */
			/*
				@param
				int i - �����
				int axis - �������������� ���� 0 - x, 1 - z 2 - y
				ret int: �������� ������ � �����
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
				grid.Dim = total; // ������������� ��������� ����������� 
				grid.Grid = new Point[total];
				// ��������� ��� �������� ��������� 
				Point* Grid = grid.Grid;
				int Nx = baseGrid.Nx;
				int Nz = baseGrid.Nz;
				int Ny = baseGrid.Ny;
				BaseGrid::PointXZS** BaseGridXZ = baseGrid.BaseGridXZ;

				int GlobalNx = grid.Nx;
				int GlobalNz = grid.Nz;
				int GlobalNy = grid.Ny;
				double* LineX = new double[GlobalNx]; // ������ ��������� � ������ �� �
				double* LineZ = new double[GlobalNz]; // ������ ��������� � ������ �� Z
				double* LineY = new double[GlobalNy]; // ������ ��������� � ������ �� Y
				// ����������� ���� ��������� xz, � ����� �� ������������� �� y

				/* ��������� �� x � z */
				/* ����������� ��������� �������� ����� � ������ �� ������������ � ����� ( ������� ����� ) */

				/* ����������� ��������� �� x ( ������� ) */

				/* ����� �������� � ���������� � �������������� ������ ��� ������������ ��������� �� z */

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
					/* ������� ��������������� �������� x �� ���� ������� */

					int startIdx = Getlevel(i, 1) * GlobalNx;
					int endIdx = startIdx + GlobalNx;
					for (int k = startIdx, kk = 0; k < endIdx; k++, kk++)
						Grid[k].x = LineX[kk];
				}


				/* ����������� ��������� �� z ( ������� ) */
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

					/* ��������� ����������� ����� � ���������� ������ */
					int startIdx = Getlevel(i, 0); // ��������� ������ ��� ������� �� ������� 
					for (int k = 0; k < GlobalNz; k++)
					{
						// ������ ����� ����� �� �������� GlobalNx
						Grid[startIdx].z = LineZ[k];
						startIdx += GlobalNx;
					}
				}


				/****************************************************/

				/* ��������� ��������������� ������������ �����  */
				/*
					������ ��������:
					�������� � ���� Z �������������� ���������� ����� ����������� �� ���� ���
					� ����� ���� �� ���� �������� ������� �����

					����� �������� ����� � ������ ������� �� ������ ���������
					������������ ������ �������� ��  ������ ����������
					������� ���������� ������ � ���������� �����
				*/
				/* ���� �� ���� �������������� ������ */
				for (int i = 0; i < GlobalNx; i++)
				{
					int idx = 0;
					/* ���� �� ����������� ��� Z */
					for (int j = 0; j < Nz - 1; j++)
					{
						int startIdx = i + GlobalNx * Getlevel(j, 1);
						int endIdx = i + GlobalNx * Getlevel(j + 1, 1);
						double left = Grid[startIdx].x; // ����� ������� �� �
						double right = Grid[endIdx].x; // ������ ������� �� � 

						// ��������� ��������� ����������� ��������� �� ��� z
						SettingForDivide param = CalcSettingForDivide(1, j, left, right);
						GenerateDivide(param, left, right, LineZ, idx);
					}
					/* ��������� ���������� � �������� ������ */

					int startIdx = i; // ��������� �������
					for (int k = 0; k < GlobalNz; k++)
					{
						Grid[startIdx].x = LineZ[k];
						startIdx += GlobalNx;
					}
				}
				/****************************************************/

				/* ��������� ��������������� �������������� �����  */
				/* ���� �� ���� �������������� ������ */
				for (int i = 0; i < GlobalNz; i++)
				{
					int idx = 0;
					for (int j = 0; j < Nx - 1; j++)
					{
						int startIdx = Getlevel(j, 0) + i * GlobalNx;
						int endIdx = Getlevel(j + 1, 0) + i * GlobalNx;
						double left = Grid[startIdx].z;
						double right = Grid[endIdx].z;
						// ��������� ��������� ����������� ��������� �� ��� x
						SettingForDivide param = CalcSettingForDivide(0, j, left, right);
						GenerateDivide(param, left, right, LineX, idx);

					}
					/* ��������� ����������� � ���������� ����� */
					int  startIdx = i * GlobalNx;
					for (int k = 0; k < GlobalNx; k++)
					{
						Grid[startIdx].z = LineX[k];
						startIdx++;
					}
				}
				/****************************************************/



				/* ��������� �� y ���������� ������� � ������ ������� ���������� */
				// �������� ��������� ��������� ������� Y
				int idxY = 0;
				for (int i = 0; i < Ny - 1; i++)
				{
					double left = baseGrid.BaseGridY[i];
					double right = baseGrid.BaseGridY[i + 1];
					SettingForDivide param = CalcSettingForDivide(2, i, left, right);
					GenerateDivide(param, left, right, LineY, idxY);
				}

				// ��������� ������� Y0 � ������� ������� 
				int idx = 0;
				for (int j = 0; j < GlobalNz; j++)
				{
					for (int k = 0; k < GlobalNx; k++)
					{
						PointInfo::ClearInfo(Grid[idx].info); // ������������� ������� ����� ��� ���� ������� ����� ������ !!!
						Grid[idx].y = LineY[0];
						idx++;
					}
				}

				// ������������� ����� �� ���� ����� 
				// i �������� �� ���������� � ������� �� Y
				// Idx �� ���������� �� ��������� 


				for (int i = 1; i < GlobalNy; i++)
				{
					int ListIdx = 0;
					for (int j = 0; j < GlobalNz; j++)
					{
						for (int k = 0; k < GlobalNx; k++)
						{
							PointInfo::ClearInfo(Grid[idx].info); // ������������� ������� ����� ��� ���� ������� ����� ������ !!!
							Grid[idx].y = LineY[i];
							Grid[idx].x = Grid[ListIdx].x;
							Grid[idx].z = Grid[ListIdx].z;
							ListIdx++;
							idx++;
						}
					}
				}

				/* ������� ������ */
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
				int AreaNum = -1; // ����� ���������� (���������� ����� ������ ���������� ���������� ��)  
				int PlaneXZSize = 0; // ������ ������� PlaneXZ
				int* PlaneXZ; // ������ ����� �������������� (������)

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
				/* ������ ������� ����� � ������� ��� ������ ����� Nx*Nz
				   ������ ������� ����� �������������� ��������������� ����� �� �����������
				*/
				int* BoundGrid = new int[grid.Nx * grid.Nz];
				/* �������� ���������� � ����������� � �� ������� */
				/*
					@param
					int i - ����� ���������� (����������) � ������� ������������ �������
					ret BoundArea - �������������� ������ ����������
				*/
				auto GetBound = [&](int i) -> BoundArea
				{
					int GlobalNx = grid.Nx;
					int GlobalNz = grid.Nz;
					BoundArea Bound;

					Bound.AreaNum = baseGrid.CalculationArea[i][0]; // ��������� ����� ���������� 

					/* �������� ��� ������ ������������� ������ ������� */
					/* �������� ������ �������������� � ���������� ��������� */
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
					/* � ���� ���������� �� ��� ����� ����� ������ ������ � ��������� �������� */
					Bound.StartY = baseGrid.BaseGridY[baseGrid.CalculationArea[i][5]];
					Bound.EndY = baseGrid.BaseGridY[baseGrid.CalculationArea[i][6]];
					Bound.PlaneXZSize = idxBoundGrid;
					Bound.PlaneXZ = BoundGrid;

					return Bound;
				};

				/* �������� ����� �� ������� */
				/*
					@param
					const BoundArea& Bound - ������ ������
					int numPointGlobal - �������� ��� ������
					ret bool: true - ������� � ������� ���� false � ��������� ������

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

						if (numPointGlobal < arr[midd])       // ���� ������� ������ �������� � ������
							right = midd - 1;      // ������� ������ ������� ������
						else if (numPointGlobal > arr[midd])  // ���� ������� ������ �������� � ������
							left = midd + 1;    // ������� ����� ������� ������
						else                       // ����� (�������� �����)
							break;           // ������� ���������� ������ ������

						// ���� ������� ����������
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

				/* ��������� ����������� �� ����� �������� ������� */
				/*
					@param
					const BoundArea &Bound - �������� �������
					Point &point - ����� ��� ��������
					int numPointGlobal - ������ ����� � ���������� ���������
					ret bool: true - ����� ����������� �������� �������, false � ��������� ������
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

					numPointGlobal %= grid.Nx * grid.Nz; // �������� ����� �� ��������� 
					bool arg1 = BinarySerch(Bound, numPointGlobal);
					bool arg2 = false;
					if (le(point.z, Bound.EndY) && ge(point.z, Bound.StartY)) arg2 = true;

					return arg1 && arg2;

				};


				/* � ����� �� ���� �������� */
				int N = grid.Dim;
				for (int i = 0; i < baseGrid.L; i++)
				{
					BoundArea Bound = GetBound(i);

					/* �� ��� ����� ������� ��� ������������� �������� (������� ���������)*/
					/* ��� ����� ������� ����������� � ������������ �������� �������� */
					// int leftAreaBound = ;
					// int RightAreaBound = ;

					/* ���� ���� ����� �� �� ���� ��������� ���� */
					for (int j = 0; j < N; j++)
					{
						if (IsInArea(Bound, grid.Grid[j], j))
						{
							/* �� � �������� ������� ������������� ������ ��������� */
							PointInfo::SetFictitiousBit(grid.Grid[j].info, 1);
							PointInfo::SetAreaInfo(grid.Grid[j].info, Bound.AreaNum);
						}
					}
				}

				/* ������� ��������� ����������  */
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


			//� ����� �� ���� ��������� ����� ���������� �������������� �����-���� ��������� ������� 

		}

		void DivisionIntoSubBounds(const BaseGrid& baseGrid, GridS& grid)
		{

			struct Bound
			{
				int Size = 0;
				int* Plane = nullptr; // ������ ����� ������� 
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

			/* �������� ����� �� ������� */
				/*
					@param
					const BoundArea& Bound - ������ ������
					int numPointGlobal - �������� ��� ������
					ret bool: true - ������� � ������� ���� false � ��������� ������

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

					if (numPointGlobal < arr[midd])       // ���� ������� ������ �������� � ������
						right = midd - 1;      // ������� ������ ������� ������
					else if (numPointGlobal > arr[midd])  // ���� ������� ������ �������� � ������
						left = midd + 1;    // ������� ����� ������� ������
					else                       // ����� (�������� �����)
						break;           // ������� ���������� ������ ������

					// ���� ������� ����������
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
				/* ������� ������� �������� ����������� ��� �� ��� � ������� � ��� ���� �����������, ��� ���� �� ��������� ����������� */
				auto GetBound = [&](int i) -> Bound
				{
					int GlobalNx = grid.Nx;
					int GlobalNz = grid.Nz;
					int GlobalNy = grid.Ny;

					Bound bound;
					bound.BoundType = baseGrid.BoundsArea[i][1];
					bound.BoundFormula = baseGrid.BoundsArea[i][0];

					/* ��� ������������ ����������  ��������������� �����*/
					/* y - ���������� */
					if (baseGrid.BoundsArea[i][6] == baseGrid.BoundsArea[i][7])
					{
						/*  ��������� �������� ��� ��� y */
						int startBaseIdx = GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // ������ ���� ��� 

						/* ���������� ������� �������� */
						int leftStartX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1);
						int leftEndX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1);
						int rightStartX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1);
						int rightEndX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1);

						/* ������� ������  */
						int nX = (rightEndX - rightStartX) / GlobalNx;
						int nZ = (rightStartX - leftStartX);

						int idxBoundGrid = 0;
						for (int i = 0; i <= nX; i++)
						{
							int Idx = leftStartX + i * GlobalNx + startBaseIdx; // ��������� ������ 
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
					/* � - ���������� */
					else if (baseGrid.BoundsArea[i][2] == baseGrid.BoundsArea[i][3])
					{
						/* ��������� ������� ���� �� ������ OZ � ����� ������������� ���� �� ������� ( ��������� ������������� )  ��� �� ������� GlobalNx*/
						int StartPozitionZ = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // ��������� ������� ���  Z
						int EndPozitionZ = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][5], 1); // �������� ����� ��� ��� Z

						int StartPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // ��������� ������� ��� Y
						int EndPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][7], 2); // ������� ������� ��� Y

						int nZ = (EndPozitionZ - StartPozitionZ) / GlobalNx; //  ���������� ����� �� ��� Z
						int nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // ���������� ����� ��� Y

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
					/* z - ���������� */
					else if (baseGrid.BoundsArea[i][4] == baseGrid.BoundsArea[i][5])
					{
						/* ���������� ������� ���� �� ��� OX, � ����� ���������� ������������� ���� �� ������� (��������� �������������) ��� �� ������� + 1 */
						int StartPositionX = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // ��������� ������� �� ��� �
						int EndPositionX = Getlevel(baseGrid.BoundsArea[i][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[i][4], 1); // �������� ������� �� ��� �

						int StartPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][6], 2); // ��������� ������� ��� Y
						int EndPositionY = Getlevel(baseGrid.BoundsArea[i][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[i][7], 2); // ������� ������� ��� Y ����� ����� (���� �������� �� ������� �������������� �����������) 

						int nX = EndPositionX - StartPositionX; //  ���������� ����� �� ��� X
						int nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // ���������� ����� ��� Y 

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


				/* ���� �� ���� �������� */
				int N = grid.Dim;
				for (int i = 0; i < baseGrid.P; i++)
				{
					/* �������� ������ ����� �������  */
					Bound bound = GetBound(i);
					for (int j = 0; j < N; j++)
					{
						/* ����������� ������ �������� � ��������������� ����� ����� */
						if (IsInBound(bound, j))
						{
							/* �� � �������� ������� ������������� ������ ��������� */
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

	/* ������ � ��������� � ����������� ������������� �������������� ���������� � ����� � �������  */
	namespace PointInfo
	{
		/* ������� ���������� �������� ����������� */
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
				/* ���� �� ���� �������� */
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
						// ��������� ��������� 
						if ((Bound.Cond[i] ^ val) == 0) // ���� �������� 
							return;
					}

					/* ����� ������� */
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
						// ��������� ��������� 
						if ((Area.Cond[i] ^ val) == 0) // ���� �������� 
							return;
					}

					/* ����� ������� */
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