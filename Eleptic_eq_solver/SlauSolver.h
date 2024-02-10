#pragma once
/* ��������� ����������� ������� ��� ������ */
#include <string>
namespace SLAUSolvers
{

	enum class MatrixType
	{
		NONE, // �� ������������� ��� - default
		DENSE,
		SPARSE,
		SPARSE_SYMETRIC,
		PROFILE,
		PROFILE_SYMETRIC
	};

	/* ������ ������ ��������� */
	namespace ForwardSolvers
	{

	};

	/* ������ ������������ ��������� �������� � ���� ��������� ������:
		��� ��� ������������ � �� ������������
		��� ��� ������������ � �� ������������
		��� � ������������ �������������������
		��� � ������������ �������������������
		��� � ������������� LLT, LUsq, LU
		��� � ������������� LLT, LUsq, LU

		������ �������� ������� ����������� ������� ����������  
	*/
	namespace IterationSolvers
	{
		using namespace std;
		struct InitDateFile
		{
			// ��� ������� 
			MatrixType type = MatrixType::NONE;

			string gg = "";
			string ggl = "";
			string ggu = "";
			string di = "";
			string kuslau = "";
			string f = ""; // ������ ������ ����� 
			string ig = "";
			string jg = "";
		};

		/* ��������� ������ */
		struct Sparse_matrix
		{
			int N = -1; // ������ ������� ���� 
			int size = -1; // ������ ������  ggl, ggu, jg
			double* ggu = nullptr; // size
			double* ggl = nullptr; // size
			double* di = nullptr; // N
			int* ig = nullptr; // N+1
			int* jg = nullptr; // size
		};

		struct Sparse_matrix_symetric
		{
			int N = -1; // ������ ������� ���� 
			int size = -1; // ������ ������  ggl, ggu, jg
			double* gg = nullptr;
			double* di = nullptr; // N
			int* ig = nullptr; // N+1
			int* jg = nullptr; // size
		};

		/* ��� ��������� ��������� ��������� ��� ���������� ����� ���� ���� ��� �� �������� ������� */
		using Fact_matrix_symetric = Sparse_matrix_symetric;
		using Fact_matrix = Sparse_matrix;
		
		/******************/

		/* ��������� ������� ���������� */
		enum class Solvemode
		{
			NONE, // �� ������������� ������ - default
			/* ��� ��� ������������ */
			LOS_SYMETRIC_CLASSIC,
			LOS_SYMETRIC_DIAG_FACT,
			LOS_SYMETRIC_LLT_FACT,

			/* ��� ��� �� ����������� */
			LOS_NOSYMETRIC_CLASSIC,
			LOS_NOSYMETRIC_DIAG_FACT,
			LOS_NOSYMETRIC_LUsq_FACT,
			LOS_NOSYMETRIC_LU_FACT,

			/* ��� ��� ������������ */
			MSG_SYMETRIC_CLASSIC,
			MSG_SYMETRIC_DIAG_FACT,
			MSG_SYMETRIC_LLT_FACT,
			/* ��� ��� �� ������������ ?*/
		

			/* ����� ���� ��� ������ ������� ���� � ����������� �������� */

		};



		/* ��������� ���� */
		struct SLAU
		{
			/*  ��������� ���� ��������� (������ � �������) !!! ��� ������� �� ��� ������   !!! */
			double *MainMemoryPool_ = nullptr; //  ��� ������ ��� ���� �������� ������� ���� � ������� ������ �����
			double *FactMemoryPool_ = nullptr; // ��� ������ ��� �������� ������ ������������ ( ��� �� � �������� ������ �� �������� ��������� ����� ���� ������ ������ ������������)
			double *SuportMemoryPool_ = nullptr; // ��� ������ ��� ��������������� �������� ( ��������� ���������� ������� ��� ������ �������  Allocate_Memory_Support() )
			
			/* ����� ���������� �����  */

			/* ���������������� ����  */
			int N = -1; // ����������� ����
			int size = -1; // ���������� ��������� ��� �������� gg, ggl, ggu, jg
			double eps = 1e-7; // �������� �������
			int maxiter = -1; // ������������ ���������� �������� 
			Sparse_matrix matrix; // �� ������������ �������
			Sparse_matrix_symetric matrix_s; // ������������ �������
			Fact_matrix_symetric Fmatr_s; // ��������������� ������� ������������
			Fact_matrix Fmatr; // ��������������� ������� �� ������������

			double* f = nullptr; // ������ ������ ����� 
			double* x = nullptr; // ������ ���������� ����������� ���������������� ������ ��� ������ ���������� ������� � ����� ��������� ���������� �����������

			/* ����� ����������������� ����� */

			/* ������������� ��� �� ���� � ������ �� ������� */
			Solvemode mode = Solvemode::NONE;

			/* ������� �������� ������ ��� �������� ��������� */
			friend ostream& operator <<(ostream& cnt, const SLAU& slau);
			
		};

		/*****************/

		
		/*******************/

		/* ������� ������ � ������� (��� ���� ���������� ��������� ������ �� ����� ����) */

		/* ��������� ������ ��� ���� */
		/*@param
			SLAU & - ���������� ���� 
			const int N - ����������� ���� 
			const int size - ����������� �������� gg, jg, ggu, ggl
			const MatrixType -  ��� ���������������� ������� 
		*/
		void Allocate_Memory_SLAU(SLAU &slau, const int N, const int size, const MatrixType type);


		/* ����� ������� ������ � ������� , ����� ���� ����������������.*/
		/* ��������� ������ ��� ������� ������������ */
		/*@param
			SLAU & - ���������� ���� 
		*/
		void Allocate_Memory_Fact(SLAU &slau);

		/* ��������� ������ ��� ��������������� ������� ��� ���������� �������.  ������������� ����� ���� �������� ���������� ��� ����� ������������ ������ ������ N */
		/*@param
			SLAU & - ���������� ���� 
		*/
		void Allocate_Memory_Support(SLAU& slau);

		/* ������� ���� ���������� ������ */
		/*@param
			SLAU & - ���������� ����
		*/
		void DeAllocateMemory(SLAU &slau);

		/********************************/
		
		

		/* ��������� ���� �� ������ */
		/* @param
		* InitDateFile & - ����� ������������� 
		* SLAU & - ��������� ����
		*/
		void Load(InitDateFile &initfile, SLAU &slau, bool printProgress = false);


		/* ������������� ��� ��� ������� ���� � �������� ������ ��� ���������� ���������������� ������ ���� ����� ��� ������ �� ��� ��� ���� ���������������� ���������� 
			� ���������� ����� �������� ������ 
		*/
		/* @param
			SLAU& - ���� 
			Solvermode - ����� ������� 
		*/
		void SetSolveMode(SLAU& slau, Solvemode mode);

		/* DEBUG functions */
		/* Ceck ������� ��������� ��� �� ������ ���� �������������������� ��� ��������� ������� � ��������� �������� �� �� ����������� ������  */

		/* ������� ������ � ��������� � ��������� */
		
		/* ��������� ������� �� ������ ������������ 
			@param
			const Sparse_matrix_symetric& matr - ������� ������������ � ����������� ������� ���������� �������
			const double* x - ������ �� ������� ������� 
			double *res - ������ ���������� 
		*/
		void MultA(const Sparse_matrix_symetric& matr, const double* x, double *res);

		/* ��������� ������� �� ������ ��� �� ������������ ������� 
			@param
			const Sparse_matrix& matr - ������� � ����������� ������� ���������� �������
			const double* x - ������ �� ������� ������� 
			double *res - ������ ���������� 
		*/
		void MultA(const Sparse_matrix& matr, const double* x, double* res);

		/* ���������� �������� ��� ��������� ����: y = a*x1 + b*x2 , ��� a,b - ��������� ����� � x1, x2 - �������, y = ������ ���������� */
		/* 
			@param
			const double a - ��������� �� ������� ���������� ������ ������
			const double* x1 - ������ ������
			const double b - ��������� �� ������� ���������� ������
			const double *x2 - ������ ������ 
			double *y - ������ ��������� 
			const int N -  ������ �������� 
		*/
		void ActionVec(const double a, const double* x1, const double b, const double *x2, double *y ,const int N);

		/* ��������� ������������ �������� */
		/* 
			@param
			const double* x1 - ������ ������
			const double* x2 - ������ ������
			const int N - ����� �������� 
			@ret
			double -  ��������� ������������ ��������
		*/
		double ScalarMult(const double* x1, const double* x2, const int N);

		/* ����� ������� */
		/* 
			@param
			const double* x - ������
			const int N - ������ ������� 
			@ret
			double - ����� ������� � E��������� ����� 
		*/
		double Norma(const double* x, const int N);


		/* ����������� �������� */
		/*
			@param
			const double* from - ������ ��������
			double* to - ���� �������� 
			const int N - ������ ������� 
		*/
		void CopyVec(const double* from, double* to, const int N);


		/* ������ ��� ����� */
		/*  
			@param
			Sparse_matrix_symetric& matr - ������� ��� ������� �� �����������  
			const double* b - ������ ������ ����� 
			double *res - ������ � ������� ��������� ��������� 
		*/
		void normal(Fact_matrix_symetric& matr,const double* b, double *res);

		/*
			@param
			Sparse_matrix& matr - ������� ��� ������� �����������
			const double* b - ������ ������ �����
			double *res - ������ � ������� ��������� ���������
		*/
		void normal(Fact_matrix& matr,const double* b, double *res);

		/* �������� ��� ����� */
		/*  
			@param
			Sparse_matrix_symetric& matr - ������� ��� ������� �� �����������  
			const double* x - ������ ������ ����� 
			double *res - ������ � ������� ��������� ��������� 
		*/
		void reverse(Fact_matrix_symetric& matr, const double* x, double* res);

		/*
			@param
			Sparse_matrix& matr - ������� ��� ������� �� �����������
			const double* x - ������ ������ �����
			double *res - ������ � ������� ��������� ���������
		*/
		void reverse(Fact_matrix& matr, const double* x, double* res);

		/********************************/

		/* ������������ */

		/* ������������ ������������ � ������������ � ��������� ����� ���������� ������������ � �� ������ ����� �������� ����������� ������ 
			��������� ��� ��� ������������ ������� ��� � ��� �� ������������ 
		*/
		/*
			@param
			SLAU& slau - ����
		*/
		void DiagFactor(SLAU& slau);

		/* LLT ������������ */
		/* ����������� �������� ����������� ���������� 
			��������� ��� ������������ ������ � �������������� ���������� �� ��������� � �� �������� �������� �������� 
		*/
		/*
			@param
			SLAU& slau - ����
		*/
		void LLTFactor(SLAU& slau);

		/* LUsq ������������ */
		/* ������������ ���������� LUsq  
			�������� ��� �� ������������ ������ � �������������� ���������� �� ��������� ( ����� ����� �������� ���������� � ������������  )
		*/
		/*
			@param
			SLAU& slau - ����
		*/
		void LUsqFactor(SLAU& slau);

		/* LU ������������ */
		/*  ������������ ���������� LU
			�������� ��� ������������� ������ � ���������� �������� ��������
		*/
		/*
			@param
			SLAU& slau - ����
		*/
		void LUFactor(SLAU& slau);
		/********************************/


		/* ������� �������� ��� ������� ���� �������� */

		/* MSG - ����������� ���������� 
			������������ - ��� ������ ����������� � ������������ ������������ 

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������) 
			bool printIteration = false -  ������������� �� �������� � �������� �� ������� 
			@ret
			double - ������� ����� ������� ������� 
		*/
		double MSG_Symetric_Classic(SLAU &slau, bool printIteration = false);

		/* MSG - � ������������ ������������������� ��� ������������ ������ 
			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������) 
			bool printIteration = false -  ������������� �� �������� � �������� �� ������� 
			@ret
			double - ������� ����� ������� ������� 
		*/
		double MSG_Symetric_DiagFact(SLAU& slau, bool printIteration = false);


		/* MSG - � �������� ������������� �� ������ ����������
			���������� ��� ����������� ������������ ������������ ������ 

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������) 
			bool printIteration = false -  ������������� �� �������� � �������� �� ������� 
			@ret
			double - ������� ����� ������� ������� 
		
		*/
		double MSG_Symetric_LLTFact(SLAU& slau, bool printIteration = false);
		
		/* LOS - ����������� ����� ��� ������������ ������
			���������� ��� ����������� ������������ ������������ ������

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������)
			bool printIteration = false -  ������������� �� �������� � �������� �� �������
			@ret
			double - ������� ����� ������� �������

		*/
		double LOS_Symetric_Classic(SLAU& slau, bool printIteration = false);

		/* LOS - � ������������ ������������������� ��� ������������ ������
			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������)
			bool printIteration = false -  ������������� �� �������� � �������� �� �������
			@ret
			double - ������� ����� ������� �������
		*/
		double LOS_Symetric_DiagFact(SLAU& slau, bool printIteration = false);


		/* LOS - � �������� ������������� �� ������ ����������
			���������� ��� ����������� ������������ ������������ ������ 

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������) 
			bool printIteration = false -  ������������� �� �������� � �������� �� ������� 
			@ret
			double - ������� ����� ������� ������� 
		*/
		double LOS_Symetric_LLTFact(SLAU& slau, bool printIteration = false);

		
		/********************************/
		/* �������� ��� �� ������������ ������ */

		/* LOS - ����������� ����� ��� �� ������������ ������
			���������� ��� �� ����������� ������������ ������������ ������

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������)
			bool printIteration = false -  ������������� �� �������� � �������� �� �������
			@ret
			double - ������� ����� ������� �������

		*/
		double LOS_Classic(SLAU& slau, bool printIteration = false);

		/* LOS - � ������������ ������������������� ��� ������������ ������
			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������)
			bool printIteration = false -  ������������� �� �������� � �������� �� �������
			@ret
			double - ������� ����� ������� �������
		*/
		double LOS_DiagFact(SLAU& slau, bool printIteration = false);

		/* LOS - � �������� ������������� �� ������ LUsq
			���������� ��� �� ����������� ������������ ������������ ������

			@param
			SLAU &slau - ��������� ���� � ������������� ����������� (�������� ������������ �� ����������)
			bool printIteration = false -  ������������� �� �������� � �������� �� �������
			@ret
			double - ������� ����� ������� �������
		*/
		double LOS_LUsqFact(SLAU& slau, bool printIteration = false);


		/* ������� �������� */

		double SolveSLAU(SLAU& slau, bool printIteration = false);

		/********************************/

		/* ������ � �������� X */
		/* �������� ������� X  �� ����� */
		/* 
			@param
			SLAU& slau - ���� 
			const string filename - ���� �� �������� ������ 
		*/
		void LoadX(SLAU& slau, const string filename);
		
		/* ��������� ������ X � ���� */
		/*
			@param
			SLAU& slau - ����
			const string filename - ���� ���� ��������� ��������� 
		*/
		void SaveX(SLAU& slau, const string filename, const string delimetr = " ");
		
		/* ���������� ������� �������  */
		void PrintX(SLAU& slau, const char *fmt = "%.7f ");

		/********************************/
	};
};