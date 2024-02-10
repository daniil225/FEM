#include "Form_Matrix.h"
#include "Grid.h"
#include <iostream>
#include <iomanip>
using namespace FormMatrix;
using namespace Grid;


double PrintDenseSLAU(SLAU& slau, int i, int j);

auto f = [](double x, double y, double z) -> double { return x*x + y*y*sin(z); };


void GenTrueq(GridS& grid, double *x)
{
	auto f = [](double x, double y, double z) -> double { return x*x + y * y * sin(z); };

	int idx = 0;
	for (int i = 0; i < grid.Dim; i++)
	{
		x[i] = f(grid.Grid[i].x, grid.Grid[i].y, grid.Grid[i].z);
	}
}



int main()
{

	setlocale(LC_ALL, "RU_ru");


	BaseGrid base;
	GridS grid;
	SLAU slau;

	Load("CalcArea2.txt", base);
	DivideGrid(base, 2);
	GenerateGrid(base, grid);
	
	
	int count = (grid.Nx - 1) * (grid.Ny - 1) * (grid.Nz - 1);
	std::cout << count << "\n";
	//Element el = GetFinitElement(grid, 0);
	//PrintElement(el);



	ParamOfDE test1 = InitTest1();
	GenerateSLAU(slau, grid, test1);



	SLAUSolvers::IterationSolvers::SetSolveMode(slau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_LUsq_FACT);
	slau.maxiter = 50;
	slau.eps = 1e-15;
	SLAUSolvers::IterationSolvers::SolveSLAU(slau, true);




	
	/* Ввиду недостатка времени костыльный вариант расчета нормы */
	/* Пусть в частном случае мы находимся на параллелипипеде  */
	double norm = 0.0;
	/* Линейные базисные функции */
	std::function<double(double, double, double)> X[2] =
	{
		[](double x, double xp, double xp1) -> double {
			if (x >= xp && x <= xp1) return (xp1 - x) / (xp1 - xp);
			else return 0.0;
		},
		[](double x, double xp, double xp1) -> double { 
			
			if (x >= xp && x <= xp1) return (x - xp) / (xp1 - xp);
			else return 0.0;
		}
	};
	
	
	
	/* Расчет функции решения ДУ  */
	std::function<double(double, double, double)> UH = [&](double x, double y, double z) -> double
		{
			double res = 0.0;
			for (int i = 0; i < count; i++)
			{
				Element el = GetFinitElement(grid, i);
				double xp = el.e[0].x;
				double xp1 = el.e[1].x;
				double ys = el.e[0].z;
				double ys1 = el.e[2].z;
				double zr = el.e[0].y;
				double zr1 = el.e[4].y;

				std::function<double(double, double, double)> Psii[8] =
				{
					[&](double x, double y, double z) -> double { return X[0](x, xp, xp1) * X[0](z, ys, ys1) * X[0](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[1](x, xp, xp1) * X[0](z, ys, ys1) * X[0](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[0](x, xp, xp1) * X[1](z, ys, ys1) * X[0](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[1](x, xp, xp1) * X[1](z, ys, ys1) * X[0](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[0](x, xp, xp1) * X[0](z, ys, ys1) * X[1](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[1](x, xp, xp1) * X[0](z, ys, ys1) * X[1](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[0](x, xp, xp1) * X[1](z, ys, ys1) * X[1](y, zr, zr1); },
					[&](double x, double y, double z) -> double { return X[1](x, xp, xp1) * X[1](z, ys, ys1) * X[1](y, zr, zr1); }
				};

				std::function<double(double, double, double)> uh = [&](double x, double y, double z) -> double
					{
						double res = 0.0;
						for (int i = 0; i < 8; i++)
						{
							res += slau.x[el.GlobalIdx[i]] * Psii[i](x, y, z);
						}
						return res;
					};

				res += uh(x, y, z);
			}
			return res;
		};

	/* Расчет нормы */
	for (int i = 0; i < count; i++)
	{
		Element el = GetFinitElement(grid, i);
		double xp = el.e[0].x;
		double xp1 = el.e[1].x;
		double ys = el.e[0].z;
		double ys1 = el.e[2].z;
		double zr = el.e[0].y;
		double zr1 = el.e[4].y;

		double hx = xp1 - xp;
		double hy = ys1 - ys;
		double hz = zr1 - zr;
		double func = UH(xp + hx / 2.0, ys + hy / 2, zr + hz / 2) - f(xp + hx / 2.0, ys + hy / 2, zr + hz / 2);
		double integ = hx * hy * hz * func * func;
		norm += integ;
	}


	std::cout << "|| u - uh ||L2 = " << std::setprecision(10) << std::sqrt(norm);
	


	GenTrueq(grid, slau.f);
	//SLAUSolvers::IterationSolvers::MultA(slau.matrix, slau.x, slau.x);
	for (int i = 0; i < slau.N; i++)
	{
		slau.x[i] -= slau.f[i];
	}

	//std::cout << "\n";
	//for (int i = 0; i < slau.N; i++)
	//{
	//	for (int j = 0; j < slau.N; j++)
	//	{
	//		std::cout << std::setw(6) << std::setprecision(2) << PrintDenseSLAU(slau, i, j) << " ";
	//	}

	//	std::cout << std::setw(7) << " = " << slau.f[i];
	//	std::cout << "\n";
	//}


	SLAUSolvers::IterationSolvers::SaveX(slau, "res.txt", "\n");
	

	SLAUSolvers::IterationSolvers::DeAllocateMemory(slau);
	DeallocateBaseGrid(base);
	DeallocateGrid(grid);
	

}

double PrintDenseSLAU(SLAU& slau, int i, int j)
{
	double* di = slau.matrix.di;
	double* ggu = slau.matrix.ggu;
	double* ggl = slau.matrix.ggl;
	int* ig = slau.matrix.ig;
	int* jg = slau.matrix.jg;
	int N = slau.N;

	if (i == j)
	{
		return di[i];
	}

	
	if (i < j)
	{
		int ind;
		int igj0 = ig[j];
		int igj1 = ig[j + 1];
		for (ind = igj0; ind < igj1; ind++)
			if (jg[ind] == i) return ggu[ind];
		
	}
	else
	{
		int ind;
		int igi0 = ig[i];
		int igi1 = ig[i + 1];
		for (ind = igi0; ind < igi1; ind++)
			if (jg[ind] == j) return ggl[ind];
		
	}

	return 0;
}