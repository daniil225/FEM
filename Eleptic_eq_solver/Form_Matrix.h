#pragma once

/* ����� ��� ���� ��������������� ����� ��������� ��� ���������� ������ ��� �������� ����� �������� � ��� ���������� */

/* ���������� ��� ������ � ��������� ������� ������� ��������� ���� */
#include "SlauSolver.h"
#include "Grid.h"
#include "ParamOfDE.h"


namespace FormMatrix
{
	using namespace SLAUSolvers::IterationSolvers;
	using namespace Grid;


	/*	��������� ����������������� ����
	* @param
	* SLAU& slau - ����
	* const GridS &grid - ����������������� �����
	* const ParamOfDE & Param - ��������� ������������ ������� �� 
	* ret: void
	* 
	* ��������������� ������� �������� � ����������� ����������������� ������� ��� !!! �� ������������ - ������������� ����� �������� �������������� ��������  
	*/
	void GenerateSLAU(SLAU& slau, const GridS& grid, const ParamOfDE & Param);

	namespace PrivateFormMatrix
	{

		/* ��������� ����� ������ */
		struct LocalPartMatrix
		{
			/* ������������ ��� ���������� �������� � ��������� �������������� ������� */
			double sign_a0;
			double a0;
			double a1;
			double a2;

			/* ������������ ���������� � ������� ��������� ������ ��� ���������� �� ���������*/
			double b1;
			double b2;
			double b3;
			double b4;
			double b5;
			double b6;

			/* ����������� �������� ������� */
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


			std::function<double(double, double)> J; // �������

			/* ��������� ����� ������� ����� */
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

			// ��������� ������� ����� �� ������ M = sign(a0)*(a0/36*M_xz0 + a1/72*M_xz1 + a2/72*M_xz2)
			// a0, a1, a2 ������� �� ������� 
			double M_xz[4][4];

			/* ��������� ������� ������������ ��������� ����������� ���������  */
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


		/* ��������� ������� � ������ F */
		struct Local
		{
			double LocalMassMatrix[Element::FinitElementSize][Element::FinitElementSize]; // ������� ����� ��� ��������� ��������� �������� �� ��� ����� ��������� ��� ����������� ������ 
			double LocalStiffnessMatrix[Element::FinitElementSize][Element::FinitElementSize]; // ������� ��������� 

			double LocalMatrix[Element::FinitElementSize][Element::FinitElementSize]; // �������������� ������� ������� ����� ��������� � ���������� 
			double LocalF[Element::FinitElementSize]; // ��������� ������ F � ��� ����� ���� 2-� � 3-� �� + ���� ������ ����� � ��
			int LocalToGlobal[Element::FinitElementSize]; // ����������� ��������� ���������� � ���������� 
		};
		
		/*	��������� �������� ������� ����
		* @param
		* SLAU& slau - ����
		* const GridS &grid - ����������������� �����
		* ret: void
		*/
		void GeneratePortrait(SLAU& slau, const GridS& grid);

		/* ���������� ��������� ������� � ������� � ����������
		* @param
		* SLAU& slau - ���������� ������� � ����� ������ ����� 
		* const Local& local - ��������� ������� 
		*  ret: void 
		* ��������� �������� ��������� ������� � ������ f � ���������� 
		*/
		void AddLocal(SLAU& slau, const Local& local);
		
		/* ��������� �������� �������� � � ���������� ������� � ������� Aij
		* @param
		* SLAU& slau - ������� ���� ������� 
		* int i - ������ 
		* int j - ������
		* double a - ��������� �������� 
		* ret: void 
		* ��������� ������� ������� � ������� ���� 
		*/
		void AddElement(SLAU& slau, int i, int j, double a);

		/* ������������� ������� � � ���������� ������� � �������� Aij
		* @param
		* SLAU& slau - ������� ���� ������� 
		* int i - ������ 
		* int j - ������
		* double a - ��������� �������� 
		* ret: void 
		* ��������� ������� ���������� � ������� ���� 
		*/
		void SetElement(SLAU& slau, int i, int j, double a);

		/* ��������� ���������� ��������
		* Local& local - ��������� ������� 
		* const Element &element - �������� ������� �������� � ���� ��� ��������� � �������� + ��� ���������� �� �������
		* const ParamOfDE& Param - ��������� �� 
		* ret: void 
		* ��������� ��������������� ��������� ������� � ����� ������ ����� �� �������� ��������� + �� 2 � 3 ��� ���� 
		*/
		void GenerateLocal(Local& local, const Element &element,const ParamOfDE& Param);

		/* ���� 1-� �� */
		void AddMainCond(SLAU& slau, const Element& element, const ParamOfDE& Param);

		/* ��������� ��������� ������ ��������� ������ ��� ��������� ��������� ��������
		* @param
		* const Element& element - �������� ������� 
		* ret: LocalPartMatrix - ��������� ����� ������ ��� ������ ����� ��������� ������ 
		*/
		LocalPartMatrix GenerateLocalPartMatrix(const Element& element);
	}
};