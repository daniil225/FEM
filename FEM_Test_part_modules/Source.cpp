#include <functional>
#include <iostream>
#include <cmath>


struct Point
{
	double x;
	double y;
};

struct Quadrilateral
{
	Point v1;
	Point v2;
	Point v3;
	Point v4;
};



/* Локальные функции определеяемые на Конечном элементе */
class LocalBasisFunction
{
private:
	struct QuadParametrs
	{
		double a0, a1, a2;
		double b1, b2, b3, b4, b5, b6;
	};

	double eps = 1e-7;
	Quadrilateral Quad; // Структура определяющая 4-х угольник

	QuadParametrs QuadParam;

	inline double W1(double dzeta) const noexcept { return 1 - dzeta; }
	inline double W2(double dzeta) const noexcept { return dzeta; }

	bool Check_Ksi_Eth_Param(const std::pair<double, double> KsiEth) const noexcept
	{
		if (KsiEth.first >= 0 && KsiEth.first <= 1.0 && KsiEth.second >= 0 && KsiEth.second <= 1.0) return true;
		else return false;
	}

	std::pair<double, double> Ksi_Eth_Calc(double x, double y) const noexcept
	{

		double ksi = -2.0;
		double eth = -2.0;
		/* Добавлено чисто для удобства чтения кода и сокращения записи формул */
		//double a0 = QuadParam.a0;
		double a1 = QuadParam.a1;
		double a2 = QuadParam.a2;
		double b1 = QuadParam.b1;
		double b2 = QuadParam.b2;
		double b3 = QuadParam.b3;
		double b4 = QuadParam.b4;
		double b5 = QuadParam.b5;
		double b6 = QuadParam.b6;

		Point v1 = Quad.v1;
		//Point v2 = Quad.v2;
		//Point v3 = Quad.v3;
		//Point v4 = Quad.v4;

		std::function<double(double, double)> w = [&](double x, double y) -> double { return b6* (x - v1.x) - b5 * (y - v1.y); };

		double w_xy = w(x, y);
		if (std::abs(QuadParam.a1) < eps && std::abs(QuadParam.a2) < eps) // а1 = а2 = 0 соответствует параллелограмму 
		{
			double b2b3_b1b4 = b2 * b3 - b1 * b4;
			ksi = (b3*(x - v1.x) - b1*(y-v1.y)) / b2b3_b1b4;
			eth = (b2 * (y - v1.y) - b4 * (x - v1.x)) / b2b3_b1b4;
		}
		else if (std::abs(QuadParam.a1) < eps && std::abs(QuadParam.a2) > eps) // параллельны только стороны между 1 и 2 и между 3 и 4 вершинами
		{
			ksi = (a2*(x-v1.x) + b1*w_xy) / (a2*b2 - b5*w_xy);
			eth = -w_xy / a2;
		}
		else if (std::abs(QuadParam.a1) > eps && std::abs(QuadParam.a2) < eps) // 1 и 3 и между 2 и 4
		{
			ksi = w_xy / a1;
			eth = (a1 * (y - v1.y) - b4 * w_xy) / (a1 * b3 + b6 * w_xy);
		}
		else
		{
			// Расчитаем коэффтциенты квадратоного уравнения 
			double a = b5 * a2;
			double b = a2 * b2 + a1 * b1 + b5 * w(x, y);
			double c = a1 * (v1.x - x) + b2 * w(x, y);

			double D = b * b - 4.0 * a * c;

			// Решаем квадратное уравнение
			double eth1 = -1.0;
			double eth2 = -1.0;
			double ksi1 = -1.0;
			double ksi2 = -1.0;
			if (D > eps)
			{
				double D_sqrt = std::sqrt(D);
				eth1 = (-b - D_sqrt) / (2.0 * a);
				eth2 = (-b + D_sqrt) / (2.0 * a);

				// Считаем ksi
				ksi1 = (a2 * eth1 / a1) + w_xy/ a1;
				ksi2 = (a2 * eth2 / a1) + w_xy / a1;
			}
			if (D < eps)
			{
				eth1 = -b / (2.0 * a);
				ksi1 = (a2 * eth1 / a1) + w_xy / a1;
				if (eth1 >= 0 && eth1 <= 1 && ksi1 >= 0 && ksi1 <= 1) { ksi = ksi1; eth = eth1; }
			}

			if (eth1 >= 0 && eth1 <= 1 && ksi1 >= 0 && ksi1 <= 1) { ksi = ksi1; eth = eth1; }
			else if (eth2 >= 0 && eth2 <= 1 && ksi2 >= 0 && ksi2 <= 1) { ksi = ksi2; eth = eth2; }
		}

		return std::make_pair(std::abs(ksi), std::abs(eth));
	}


public:

	LocalBasisFunction() = default;
	LocalBasisFunction(const Quadrilateral& quad_): Quad(quad_)
	{
		Point &v1 = Quad.v1;
		Point &v2 = Quad.v2;
		Point &v3 = Quad.v3;
		Point &v4 = Quad.v4;

		QuadParam.a0 = (v2.x - v1.x) * (v3.y - v1.y) - (v2.y - v1.y) * (v3.y - v1.x);
		QuadParam.a1 = (v2.x - v1.x) * (v4.y - v3.y) - (v2.y - v1.y) * (v4.x - v3.x);
		QuadParam.a2 = (v3.y - v1.y) * (v4.x - v2.x) - (v3.x - v1.x) * (v4.y - v2.y);

		QuadParam.b1 = v3.x - v1.x;
		QuadParam.b2 = v2.x - v1.x;
		QuadParam.b3 = v3.y - v1.y;
		QuadParam.b4 = v2.y - v1.y;
		QuadParam.b5 = v1.x - v2.x - v3.x + v4.x;
		QuadParam.b6 = v1.y - v2.y - v3.y + v4.y;
	}
	
	
	double Psi1(double x, double y) const noexcept 
	{
		/* Пройесс расчета значений выражения */
		std::pair<double, double> ksi_eth = Ksi_Eth_Calc(x, y);

		std::cout << "\nKsi = " << ksi_eth.first << "  Eth = " << ksi_eth.second << "\n";
		if (Check_Ksi_Eth_Param(ksi_eth)) return W1(ksi_eth.first) * W1(ksi_eth.second);
		else return -1;
	}

	~LocalBasisFunction() {}
};

int main()
{
	Quadrilateral quad;
	quad.v1.x = 1; quad.v1.y = 1;
	quad.v2.x = 3; quad.v2.y = 0.8;
	quad.v3.x = 1.3; quad.v3.y = 1.8;
	quad.v4.x = 2.7; quad.v4.y = 2.1;
	LocalBasisFunction Func(quad);

	double x = 2.6955018;
	double y = 2.1003985;

	std::cout << "x = " << x << " y = " << y << " Psi1(x,y) = " << Func.Psi1(x, y) << "\n";


	return 0;
}