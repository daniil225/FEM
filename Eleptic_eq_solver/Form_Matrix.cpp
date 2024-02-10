#include "Form_Matrix.h"
#include <iostream>


namespace FormMatrix
{

	void GenerateSLAU(SLAU& slau, const GridS& grid, const ParamOfDE & Param)
	{
		/* Выделяем память под СЛАУ */
		try
		{
			PrivateFormMatrix::GeneratePortrait(slau, grid); //  Генерация портрета 

			int N = slau.N;
			int size = slau.size;
			if (slau.MainMemoryPool_ != nullptr)
				delete[] slau.MainMemoryPool_;

			slau.MainMemoryPool_ = new double[2 * N + 2 * size];
			double* MainMemoryPool_ = slau.MainMemoryPool_;
			int index = 0; // Индексация 

			/* Распределяем память массив double */
			slau.matrix.di = &MainMemoryPool_[index];
			index += N;
			slau.f = &MainMemoryPool_[index];
			index += N;
			slau.matrix.ggl = &MainMemoryPool_[index];
			index += size;
			slau.matrix.ggu = &MainMemoryPool_[index];
			/* Инициализация памяти нулями */
			for (int i = 0; i < 2 * N + 2 * size; i++)
				MainMemoryPool_[i] = 0;

			/* Генерация матрицы и ее сборка */
			int regionsNum = (grid.Nx - 1) * (grid.Nz - 1) * (grid.Ny - 1); // Количество конечных элементов 
			
			PrivateFormMatrix::Local localData;
			for (int i = 0; i < regionsNum; i++)
			{
				Element elem = GetFinitElement(grid, i); // Получаем конечный элемент 
				
				/* Элемент не фиктивный производим расчет  */
				if (elem.isFictive)
				{
					// Генерация локадльной матрицы 
					PrivateFormMatrix::GenerateLocal(localData, elem, Param);
					// Добавление в глобальную матрицу 
					PrivateFormMatrix::AddLocal(slau, localData);

				}
				else
				{
					/* Расставляем фиктивные элементы нужно на дигональ поставить 1 */
					for (int i = 0; i < Element::FinitElementSize; i++)
					{
						/* Элемент фиктивный  */
						if (!PointInfo::IsFiFictitious(elem.e[i].info))
						{
							/* Ставим на диагональ 1 */
							int idx = elem.GlobalIdx[i];
							PrivateFormMatrix::AddElement(slau, idx, idx, 1.0);
						}
					}
				}
			}

			///* Учет первых краевых */
			for (int i = 0; i < regionsNum; i++)
			{
				Element elem = GetFinitElement(grid, i); // Получаем конечный элемент 
				
				for (int j = 0; j < elem.BoundCount; j++)
				{
					/* если граница и тип КУ 1-ый */
					if (elem.Bound[j].IsBound && elem.Bound[j].BoundType == 1)
					{
						/* В цикле по вершинам */
						Func_3D f = GetBoundCond12(Param, elem.Bound[j].BoundInfo);
						for (int k = 0; k < 4; k++)
						{
							int i_global = elem.Bound[j].GlobalIdx[k]; // Глобальный номер узла он соответсвует строке в матрице
							//std::cout << i_global << "\n";
							/* В цикле по строке матрицы выставляем нулевые значения */
							for (int l = 0; l < N; l++)
								PrivateFormMatrix::SetElement(slau, i_global, l, 0.0);
							/* На диагонали 1-а */
							PrivateFormMatrix::SetElement(slau, i_global, i_global, 1.0);
							/* В вектор f значение на границе */
							slau.f[i_global] = f(elem.e[elem.Bound[j].LocalIdx[k]].x, elem.e[elem.Bound[j].LocalIdx[k]].y, elem.e[elem.Bound[j].LocalIdx[k]].z);
						}
					}
				}
			}

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

	namespace PrivateFormMatrix
	{
		void GeneratePortrait(SLAU& slau, const GridS& grid)
		{
			try
			{
				int Nx = grid.Nx;
				int Nz = grid.Nz;
				/* По номеру конечного элемента и локальному номеру на нем возвращает глобальный индекс конечного элемента
				* @param
				* int ielem - номер конечного элемента
				* int j - локальный номер
				*/
				auto IndexOfUnknown = [&](int idx, int j) -> int
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

					int res = shiftY + shiftXZ;

					switch (j)
					{
					case 0:
						res += 0;
						break;
					case 1:
						res += 1;
						break;
					case 2:
						res += Nx;
						break;
					case 3:
						res += (Nx + 1);
						break;
					case 4:
						res += (Nx * Nz);
						break;
					case 5:
						res += (Nx * Nz + 1);
						break;
					case 6:
						res += (Nx * Nz) + Nx;
						break;
					case 7:
						res += (Nx * Nz) + Nx + 1;
						break;
					default:
						break;
					}

					return res;
				};

				int globalN = grid.Dim;
				int regionsNum = (grid.Nx - 1) * (grid.Nz - 1) * (grid.Ny - 1);

				if (slau.matrix.ig != nullptr)
					delete[] slau.matrix.ig;

				slau.matrix.ig = new int[globalN + 1];
				int* ig = slau.matrix.ig;

				int* list[2]{};
				list[0] = new int[2 * globalN * (globalN - 2)];
				list[1] = new int[2 * globalN * (globalN - 2)];
				int* listbeg = new int[globalN];

				int listsize = -1;
				for (int i = 0; i < globalN; i++)
					listbeg[i] = -1;
				for (int ielem = 0; ielem < regionsNum; ielem++)
				{
					for (int i = 0; i < 8; i++)
					{
						int k = IndexOfUnknown(ielem, i);
						for (int j = i + 1; j < 8; j++)
						{
							int ind1 = k;
							int ind2 = IndexOfUnknown(ielem, j);
							if (ind2 < ind1) {
								ind1 = ind2;
								ind2 = k;
							}
							int iaddr = listbeg[ind2];
							if (iaddr == -1) {
								listsize++;
								listbeg[ind2] = listsize;
								list[0][listsize] = ind1;
								list[1][listsize] = -1;
							}
							else {
								while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
									iaddr = list[1][iaddr];

								if (list[0][iaddr] > ind1) {
									listsize++;
									list[0][listsize] = list[0][iaddr];
									list[1][listsize] = list[1][iaddr];
									list[0][iaddr] = ind1;
									list[1][iaddr] = listsize;
								}
								else {
									if (list[0][iaddr] < ind1) {
										listsize++;
										list[1][iaddr] = listsize;
										list[0][listsize] = ind1;
										list[1][listsize] = -1;
									}
								}
							}
						}
					}
				}

				if (slau.matrix.jg != nullptr)
					delete[] slau.matrix.jg;
				slau.matrix.jg = new int[listsize + 1];
				int* jg = slau.matrix.jg;

				ig[0] = 0;
				for (int i = 0; i < globalN; i++)
				{
					ig[i + 1] = ig[i];
					int iaddr = listbeg[i];
					
					while (iaddr != -1) {
						jg[ig[i + 1]] = list[0][iaddr];
						ig[i + 1]++;
						iaddr = list[1][iaddr];
					}
				}
				slau.N = globalN;
				slau.size = slau.matrix.ig[globalN];
				slau.matrix.N = globalN;
				slau.matrix.size = slau.size;

				/* Очистка памяти */
				delete[] listbeg;
				delete[] list[0];
				delete[] list[1];
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

		void AddLocal(SLAU& slau, const Local& local)
		{
			// В цикле по элементам Локальной матрицы добавляем все элементы в СЛАУ 
			for (int i = 0; i < Element::FinitElementSize; i++)
			{
				int i_global = local.LocalToGlobal[i];
				for (int j = 0; j < Element::FinitElementSize; j++)
				{
					int j_global = local.LocalToGlobal[j];
					double a = local.LocalMatrix[i][j];
					//std::cout << "Element add to " << i_global << " " << j_global << "\n";
					AddElement(slau, i_global, j_global, a);
				}
				/* Добавляем вектор f */
				slau.f[i_global] += local.LocalF[i];
			}
		}

		void AddElement(SLAU& slau, int i, int j, double a)
		{
			double* di = slau.matrix.di;
			double* ggu = slau.matrix.ggu;
			double* ggl = slau.matrix.ggl;
			int* ig = slau.matrix.ig;
			int* jg = slau.matrix.jg;
			int N = slau.N;

			if (i == j)
			{
				di[i] += a;
				return;
			}

			int ind;
			if (i < j)
			{
				int igj0 = ig[j];
				int igj1 = ig[j + 1];
				for (ind = igj0; ind < igj1; ind++)
				{
					if (jg[ind] == i)
					{
						ggu[ind] += a;
						return;
					}
				}
				
			}
			else
			{
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (ind = igi0; ind < igi1; ind++)
				{
					if (jg[ind] == j)
					{
						ggl[ind] += a;
						return;
					}
				}
			}
		}

		void SetElement(SLAU& slau, int i, int j, double a)
		{
			double* di = slau.matrix.di;
			double* ggu = slau.matrix.ggu;
			double* ggl = slau.matrix.ggl;
			int* ig = slau.matrix.ig;
			int* jg = slau.matrix.jg;
			int N = slau.N;

			if (i == j)
			{
				di[i] = a;
				return;
			}

			int ind;
			if (i < j)
			{
				int igj0 = ig[j];
				int igj1 = ig[j + 1];
				for (ind = igj0; ind < igj1; ind++)
				{
					if (jg[ind] == i)
					{
						ggu[ind] = a;
						return;
					}
				}
			}
			else
			{
				int igi0 = ig[i];
				int igi1 = ig[i + 1];
				for (ind = igi0; ind < igi1; ind++)
				{
					if (jg[ind] == j)
					{
						ggl[ind] = a;
						return;
					}
				}
			}
		}

		void GenerateLocal(Local& local, const Element& element, const ParamOfDE& Param)
		{

			auto mu = [](int i) -> int { return i % 4; };
			auto nu = [](int i) -> int { return std::floor(i / 4); };
			// Перенесем глобальные номера 
			for (int i = 0; i < Element::FinitElementSize; i++)
				local.LocalToGlobal[i] = element.GlobalIdx[i];

			LocalPartMatrix localPartMatr = GenerateLocalPartMatrix(element);

			/* Сборка локальной матрицы из кусочков локальных  */
			
			/* Матрица массы и жесткости расчитанные с lambda = gamma = 1 */
			double hy = element.e[5].y - element.e[0].y;
			//std::cout << hy << "\n";
			for (int i = 0; i < Element::FinitElementSize; i++)
			{
				for (int j = 0; j < Element::FinitElementSize; j++)
				{
					local.LocalStiffnessMatrix[i][j] = (hy * localPartMatr.C_xz[mu(i)][mu(j)] * localPartMatr.M_y[nu(i)][nu(j)])/6.0 + (localPartMatr.M_xz[mu(i)][mu(j)] * localPartMatr.C_y[nu(i)][nu(j)])/hy;
					local.LocalMassMatrix[i][j] = (hy * localPartMatr.M_xz[mu(i)][mu(j)] * localPartMatr.M_y[nu(i)][nu(j)]) / 6.0;
					//std::cout << local.LocalMassMatrix[i][j] << " ";
				}
				//std::cout << "\n";
			}

			/* В курсаче lambda = const  и gamma = const */
			Func_3D lambda = GetLambda(Param, element.AreaInfo);
			Func_3D gamma = GetGamma(Param, element.AreaInfo);
			Func_3D F = GetF(Param, element.AreaInfo);

			for (int i = 0; i < Element::FinitElementSize; i++)
			{
				double lambda_ = lambda(element.e[i].x, element.e[i].y, element.e[i].z);
				double gamma_ = gamma(element.e[i].x, element.e[i].y, element.e[i].z);
				for (int j = 0; j < Element::FinitElementSize; j++)
				{
					local.LocalMatrix[i][j] = lambda_ * local.LocalStiffnessMatrix[i][j] + gamma_ * local.LocalMassMatrix[i][j];
					//std::cout << local.LocalMatrix[i][j] << " ";
				}
				//std::cout << "\n";
				
			}
			//std::cout << "\n";
			/* Формаируем вектор парвой части  */
			for (int i = 0; i < Element::FinitElementSize; i++)
			{
				double fi = F(element.e[i].x, element.e[i].y, element.e[i].z);
				local.LocalF[i] = 0;
				for (int j = 0; j < Element::FinitElementSize; j++)
				{
					local.LocalF[i] += fi * local.LocalMassMatrix[i][j];
				}
			}

			/* Учет 2-х и 3-х КУ */
			double BaseMassMatrix[4][4] =
			{
				{4,2,2,1},
				{2,4,1,2},
				{2,1,4,2},
				{1,2,2,4}
			};
			for (int j = 0; j < element.BoundCount; j++)
			{
				/* если граница и тип КУ 2-ые */
				if (element.Bound[j].IsBound && element.Bound[j].BoundType == 2)
				{
					/* В цикле по вершинам */
					Func_3D f = GetBoundCond12(Param, element.Bound[j].BoundInfo);
					
					/* По полученной границе сформируем Матрицу массы она состоит из локальных матриц */
					if (element.Bound[j].LocalIdx[0] == 2 && element.Bound[j].LocalIdx[1] == 3 && element.Bound[j].LocalIdx[2] == 6 && element.Bound[j].LocalIdx[3] == 7)
					{
						double hx = element.e[3].x - element.e[2].x;
						double hy = element.e[6].y - element.e[2].y;
						
						/* По всем 4-ем строкам */
						for (int i = 0; i < 4; i++)
						{
							int i_local = element.Bound[j].LocalIdx[i];
							double fi = 0.0;
							/* По 4-ем столбцам */
							for (int k = 0; k < 4; k++)
							{
								double x = element.e[element.Bound[j].LocalIdx[k]].x;
								double y = element.e[element.Bound[j].LocalIdx[k]].y;
								double z = element.e[element.Bound[j].LocalIdx[k]].z;
								fi += f(x,y,z)*hx*hy*BaseMassMatrix[i][k]/36.0;
							}
							local.LocalF[i_local] += fi;
						}
					}
					else if (element.Bound[j].LocalIdx[0] == 4 && element.Bound[j].LocalIdx[1] == 5 && element.Bound[j].LocalIdx[2] == 6 && element.Bound[j].LocalIdx[3] == 7)
					{
						double hx = element.e[5].x - element.e[4].x;
						double hz = element.e[6].z - element.e[4].z;

						/* По всем 4-ем строкам */
						for (int i = 0; i < 4; i++)
						{
							int i_local = element.Bound[j].LocalIdx[i];
							double fi = 0.0;
							/* По 4-ем столбцам */
							for (int k = 0; k < 4; k++)
							{
								double x = element.e[element.Bound[j].LocalIdx[k]].x;
								double y = element.e[element.Bound[j].LocalIdx[k]].y;
								double z = element.e[element.Bound[j].LocalIdx[k]].z;
								fi += f(x, y, z) * hx * hz * BaseMassMatrix[i][k] / 36.0;
							}
							local.LocalF[i_local] += fi;
						}
					}
					else if (element.Bound[j].LocalIdx[0] == 0 && element.Bound[j].LocalIdx[1] == 2 && element.Bound[j].LocalIdx[2] == 4 && element.Bound[j].LocalIdx[3] == 6)
					{
						double hy = element.e[4].y - element.e[0].y;
						double hz = element.e[2].z - element.e[0].z;

						/* По всем 4-ем строкам */
						for (int i = 0; i < 4; i++)
						{
							int i_local = element.Bound[j].LocalIdx[i];
							double fi = 0.0;
							/* По 4-ем столбцам */
							for (int k = 0; k < 4; k++)
							{
								double x = element.e[element.Bound[j].LocalIdx[k]].x;
								double y = element.e[element.Bound[j].LocalIdx[k]].y;
								double z = element.e[element.Bound[j].LocalIdx[k]].z;
								fi += f(x, y, z) * hy * hz * BaseMassMatrix[i][k] / 36.0;
							}
							local.LocalF[i_local] += fi;
						}
					}
				}
			}
		}

		LocalPartMatrix GenerateLocalPartMatrix(const Element& element)
		{
			LocalPartMatrix Matr;
			
			const Point* e = element.e;
			Matr.a0 = (e[1].x - e[0].x) * (e[2].z - e[0].z) - (e[1].z - e[0].z) * (e[2].x - e[0].x);
			Matr.a1 = (e[1].x - e[0].x) * (e[3].z - e[2].z) - (e[1].z - e[0].z) * (e[3].x - e[2].x);
			Matr.a2 = (e[2].z - e[0].z) * (e[3].x - e[1].x) - (e[2].x - e[0].x) * (e[3].z - e[1].z);
			Matr.sign_a0 = (std::abs(Matr.a0) - Matr.a0) < 1e-10 ? 1 : -1;

			Matr.b1 = e[2].x - e[0].x;
			Matr.b2 = e[1].x - e[0].x;
			Matr.b3 = e[2].z - e[0].z;
			Matr.b4 = e[1].z - e[0].z;
			Matr.b5 = e[0].x - e[1].x - e[2].x + e[3].x;
			Matr.b6 = e[0].z - e[1].z - e[2].z + e[3].z;

			Matr.J = [&](double Ksi, double Eth) -> double { return Matr.a0 + Matr.a1 * Ksi + Matr.a2 * Eth; };
			/* Собираем матрицу массы xz*/
			double a0 = Matr.a0 / 36.0;
			double a1 = Matr.a1 / 72.0;
			double a2 = Matr.a2 / 72.0;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Matr.M_xz[i][j] = a0 * Matr.M_xz0[i][j] + a1 * Matr.M_xz1[i][j] + a2 * Matr.M_xz2[i][j];
					//std::cout << Matr.M_xz[i][j] << " ";
				}
				
			}

			/* Генерация матрицы жескости */
			/* Точки Гауса  */
			double t[3] = { 0, 0.77459666924148337,-0.77459666924148337 };
			double tau[3] = { 8.0 / 9.0, 5.0 / 9.0, 5.0 / 9.0 };

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					/* Подинтегральная функция */
					std::function<double(double, double)> f = [&](double Ksi, double Eth) -> double
					{
						double arg1 = Matr.b6 * Ksi + Matr.b3;
						double arg2 = Matr.b6 * Eth + Matr.b4;
						double arg3 = Matr.b5 * Eth + Matr.b2;
						double arg4 = Matr.b5 * Ksi + Matr.b1;

						double f1 = ((Matr.dPhi_dKsi[i](Eth)*arg1 - Matr.dPhi_dEth[i](Ksi)*arg2) * (Matr.dPhi_dKsi[j](Eth)*arg1 - Matr.dPhi_dEth[j](Ksi)*arg2)) / Matr.J(Ksi, Eth);
						double f2 = ((Matr.dPhi_dEth[i](Ksi)*arg3 - Matr.dPhi_dKsi[i](Eth)*arg4) * (Matr.dPhi_dEth[j](Ksi)*arg3 - Matr.dPhi_dKsi[j](Eth)*arg4)) / Matr.J(Ksi, Eth);

						double f = Matr.sign_a0 * (f1 + f2);
						return f;
					};

					double Cij = 0.0;
					/* Квадратура Гауса */
					for (int k = 0; k < 3; k++)
					{
						for (int l = 0; l < 3; l++)
							Cij += tau[k] * tau[l] * f(( 1 + t[k] )/2.0,(1+t[l])/2.0);
					}

					Matr.C_xz[i][j] = Cij / 4.0;
				}
			}

			return Matr;
		}
	};
};