#pragma once

#include <string>

/* GridBuilder v1.0.0 */

/* �������� ���� ��� ��� ������� � ��������� ������ ���������� � ����� 
	������ ��� ����� ������������� 
	1) ���� ���� ���������� �� ��������� � ��������� ���������� 
	2) ����������� ������������ ������� ���������� ������� RAII
	3) ���������� ������������� ��������� ����� ����� �� ������ ������ �� ��������� ������������ 
	4) ������� ���������� ������������� ������� ����������� ��� ������������ ��������� �������������� ����� ���������� ����� ��������� ��� ������� � ������� ������ 
	5) ���� ���� �������������� ������� ��������� ��� ��������� 
*/


/* ������ � ����� �� �����. ��������� �����.  */
/* �������������� ��������� ������������ ��� �������� ����� � ���� ������������ ��������� 
   ������ ������ � ��������������� ���������� � ���������� � ��������� (x;z)
*/

/* ������ ����� 
* Nx Nz
* x1,1 z1,1 x2,1 z1,2 ... xn,1 z1,n
* x1,2 z2,1 x2,2 z2,2 ... x2,n z2,n
* .....
* Ny
* y1, y2, ..., yk
* 
* L - ����� ����� - ���������� ����������� � ����� L ������� ����� �� 7 ���� 
* 1 ����� - ����� ������ ������������ ��������� �� � ����������
* 2 ����� - ������ ������������ �������� ����������� ������ ������� 
* 3 ����� - ������ ������������ �������� ���������� ����� ������� 
  �� ������ ������������ � ������������ � �������� ��� k ���� �� ������ � �������� ����� �����: k, (Nx + k), (2*Nx + k), ... , ( (Ny-1)*Nx + k )
* 4 ����� - ����������� �������������� ������� �����
* 5 ����� - ���������� �������������� ������� ������ 
* 6 ����� - ������� �� y ������
* 7 ����� - ������� �� y �����
* ������ 3 ������� �������� ��������� ����� 
* P - ����� ����� - ���������� ����������� ��������� ������� � ����� P ������� ����� �� 8 ����
* 1 ����� - ����� ������ ������������ ��������� ��������� ������� 
* 2 ����� - ��� �������� ������� 0 - ������������ �� �������� ������� 
* 3 ����� - ������ ������������ �������� ����������� ������ ������� 
* 4 ����� - ������ ������������ �������� ���������� ����� ������� 
* 5 ����� - ����������� �������������� ������� �����
* 6 ����� - ���������� �������������� ������� ������ 
* 7 ����� - ������� �� y ������
* 8 ����� - ������� �� y �����
* !!! ������������� ����������� ���� �� ������������ ����� ������ ���� ������������� 
* N1 coef1 N2 coef 2 ...
* ...
* ...
* ������ ����� ������ ��������� ������� ������� �� �
* ������ ������ ������ ��������� �� z
* ������ ����� ������ ��������� �� �
*/

namespace Grid
{

	/* ������ � ��������� � ����������� ������������� �������������� ���������� � ����� � �������  */
	namespace PointInfo
	{
		/*
		��������� ������������ ���������� � ����� � �������
		������ ���������
		*/
		struct Info
		{
			uint64_t BaseInfo : 8; /* ������� ���������� */
			uint64_t AreaInfo : 32; /* ���������� �� ������� ������������ ��������� �� */
			uint64_t BoundInfo : 16; /*������� ������� ������ ����� ������ */
			uint64_t TypeBoundCond : 8; /* ��� �������� ������� */
		};

		/* ��������� ���������� � ���� ��  */
		struct BoundInfo
		{
			uint8_t size = 0; // ������� ����� ������ 
			uint8_t Cond[4]{ 0,0,0,0 }; // ����� ������ ������������ ��
			uint8_t TypeCond[4]{ 0,0,0,0 }; // ��� �������� ������� 
		};

		/* ��������� ���������� �� ������� */
		struct AreaInfo
		{
			uint8_t size = 0;
			uint8_t Cond[8]{ 0,0,0,0,0,0,0,0 };
		};

		struct Comand
		{
			typedef const uint32_t ComandType;
			/* Clear comand */
			static ComandType BaseInfoClear = 0x00;
			static ComandType AreaInfoClear = 0x00000000;
			static ComandType BoundInfoClear = 0x0000;
			static ComandType TypeBoundCondClear = 0x00;

			/* Base Info comand */
			static ComandType SetZeroFiFictitious = 0xFE;
			static ComandType GetFiFictitious = 0x01;

			/* Area Comand */
			static ComandType SetZeroAreaCount = 0xE1;
			static ComandType GetAreaCount = 0x1E;
			static ComandType ShiftAreaCountBits = 1;
			static ComandType RightBoundAreaCount = 8;
			static ComandType RightBoundValAreaInfo = 15;
			static ComandType BaseAreaShift = 4;
			static ComandType GetAreaInfoBaseComand = 0x0000000F;

			/* Bound Comand */
			static ComandType SetZeroBountCount = 0x1F;
			static ComandType ShiftBoundCountBits = 5;
			static ComandType GetBoundCount = 0xE0;
			static ComandType RightBoundBoundCount = 4;
			static ComandType RightBoundValBoundInfo = 15;
			static ComandType RightBoundValBoundType = 3;
			static ComandType BaseBoundShift = 4;
			static ComandType BaseTypeBoundShift = 2;
			static ComandType GetBoundInfoBaseComand = 0x000F;
			static ComandType GetTypeBoundBaseComand = 0x03;
		};

		/* ������� ���������� �������� ����������� */
		/* Private Section */

		/* ������� ��������� */
		void ClearInfo(Info& info);

		/* ��������� ���������/�� ��������� ���� 
			@param
			Info& info - �������������� ��������� 
			uint8_t val true - �� ��������� false - ���������
			ret void
		*/
		void SetFictitiousBit(Info& info, uint8_t val);
		
		/* ��������� ���������� ���������� ����� �������� ������� ����������� ����� 
			@param
			Info& info - �������������� ��������� 
			uint8_t val - �������� [0;8]
			ret void
		*/
		void SetAreaCountBits(Info& info, uint32_t val);

		/* ��������� ���������� ��������� ����� ������ � ������� ��������� ������ �����
			@param
			Info& info - �������������� ���������
			uint8_t val - �������� [0;4]
			ret void
		*/
		void SetBoundCountBits(Info& info, uint32_t val);

		/* �������� ���������� ��������� �������� � ������� ��������� �����
			@param
			Info& info - �������������� ���������
			ret uint32_t - ����������� � ��� � ������� �������� ���������  �������� [0;8]
			0 - �� � ����� 
			1 - � ����� 
			...
		*/
		uint32_t GetAreaCount(const Info& info);

		/* �������� ���������� ��������� �������� � ������� ��������� ������� 
			@param
			Info& info - �������������� ���������
			ret uint32_t - ����������� � ��� � ������� �������� ���������(���������)  �������� [0;4]
			0 - �� � ����� 
			1 - � ����� 
			2 - 
			3 - 
			4 - 
		*/
		uint32_t GetBoundCount(const Info& info);
		/*********************************************/

		/* Public section */

		/*
			@param
			Info& info 
			uint32_t val - ����� ������ �������� �������� �� ������� �������� [0,15]
			uint32_t boundType - ��� �� ��������  [0,3]
		*/
		void SetBoundInfo(Info& info, uint32_t val, uint32_t boundType);
		/*
			@param
			Info& info
			uint32_t val - ����� ������ �������� �������� � ������� ������������ ��������� ��  �������� [0,15]
		*/
		void SetAreaInfo(Info& info, uint32_t val);

		/* 
			@param
			const Info& info
			ret bool - false - ��������� true - �� ��������� 
		*/
		bool IsFiFictitious(const Info& info);

		/* 
			@param
			const Info& info
			ret bool - false - �� ������� true - ������� 
		*/
		bool IsBound(const Info& info);

		/*
			@param
			const Info& info
			ret BoundInfo -��������� ���������� ���������� � ����� �� � ������ ������� �� ������
		*/
		BoundInfo GetBoundInfo(const Info& info);

		/*
			@param
			const Info& info
			ret AreaInfo -��������� ���������� ���������� � ����� �������� ������� ����������� ����� � ��� ��������� ��� ���������� �������� �� 
		*/
		AreaInfo GetAreaInfo(const Info& info);
		/********************************************/

		/* Debug functions */
		void PrintInfo(const Info& info);

		void PrintBoundInfo(const BoundInfo& Bound);

		void PrintAreaInfo(const AreaInfo& Area);
	};
	/***********************************************************************************/

	/*  ��������� ������ ��� ������������� �������  */
	const double eps = 1e-13; // �������� ���� 

	/* ����������  ������� ����� �� �� ������ ������������ �������� ����� ��� ������� */
	struct BaseGrid
	{
		int CountOfDivision = 1; // ���������� ������� ������� ��������� �����; 1 - ������������ ���� ��� ������� �� ���� 2 ���� ��� �������� ������� ��������� � ���� �.�.�
		const int SizeOfCalculationAreaElemet = 7;
		const int SizeOfBoundsAreaElement = 8;
		const int SizeOfDivideParam = 3;

		/* ��������� ��� ������� �����  */
		struct PointXZS
		{
			double x;
			double z;
		};

		struct DivideParamS
		{
			int num; // ���������� ���������� �� ������� ����� ��������� ������� 
			double coef; // ����������� ���������� ��� ������ 
		};
		int Nx = 0; // ���������� ����� ����� ���������������� ����������� 
		int Ny = 0; // ���������� ����� ����� ������������� ����������� 
		int Nz = 0; // ���������� ����� ����� ��� z
		int L = 0; // ���������� ����������� 
		int P = 0; // ���������� ����� ������ 
		PointXZS **BaseGridXZ = nullptr; //  ������� ����� � ���������� XZ 
		double *BaseGridY = nullptr; //  ������� ����� �� Y
		/*
			| ����� ������ | ������� �� x(start) | ������� �� x(end) | | ������� �� z(start) | ������� �� z(end) | ������� �� y(start) | ������� �� y(end) |
		*/
		int** CalculationArea = nullptr; // ������ ��������� ��������
		/*
			| ����� ������ | | ��� �� | ������� �� x(start) | ������� �� x(end) | | ������� �� z(start) | ������� �� z(end) | ������� �� y(start) | ������� �� y(end) |
		*/
		int** BoundsArea = nullptr; // �������

		/*
			3 ������� 
			DivideParam[0] - ������ ��� ��������� �� ��� x ������ Nx-1
			DivideParam[1] - ������ ��� ��������� �� ��� z ������ Nz-1
			DivideParam[2] - ������ ��� ��������� �� ��� y ������ Ny-1
		*/
		DivideParamS **DivideParam = nullptr; // ������ ��� ��������� 
	};



	// �������������� ����� ��� ���� ��������� �������� 
	struct Point
	{
		// ���������� � ������� ���������� ����� ������ + ���������� �� ���� ��������� ��� ��� + ���������� � ��� ��������� �� ���� ��� ��� 
		// + ���������� � �������� �������� + ��� �� + ����� � ����������� 
		// ��� ��� ���������� �������� � �������������� ����� ����� ������� ������ ����
		
		PointInfo::Info info;

		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
	};


	// �������� ������� 
	struct Element
	{
		struct BoundS
		{
			bool IsBound = false;
			int LocalIdx[4]{ -1,-1,-1,-1 };
			int GlobalIdx[4]{ -1, -1, -1, -1 };
			int BoundType = -1;
			int BoundInfo = -1;
		};
		// ��������� ���� �� �������� � ��� ����� ���������� ����������, �� ��� ���������� ��������� ��������� ����� ���������� ����� ������� 
		// �� ����� ��������� ����� ����� �������� ��� ����������

		static const int FinitElementSize = 8; // ������ �������� 
		static const int BoundCount = 6;
		/* ���������� � �������� �������� */
		bool isFictive = false; // �� ��������� ��������� 
		int AreaInfo = -1; // ����� ������� �������� �������� �� �� �������� ������� 
		/* ������� */
		/* ����� �� �������� �������� ����� ���� 6 ������ (8-� ��������) 
			�� ��������� ��������� ���������� 
			{1,2,3,4}, {5,6,7,8}, {1,3,5,7}, {2,4,6,8}, {1,2,5,6}, {3,4,7,8} - � ��������� � 1
			�� ������ ������� ����� ���� ������ ���� �������. ��� ���� ������� ������� ������������ 2 ������� 
			1 ����� - ����� ������ �������� ��
			2 ����� - ��� ��
			������� 6 �������� ��� ������ ������� �������� 

			��������� � ��������� �������:
			bool IsBound = false; - ������� �� ���. �� ��������� ��� 
			int LocalIdx[4]; - ��������� ��������� 
			int GlobalIdx[4]; - �������������� ���������� ��������� ��� �������� ����������� � ���������� ������� 
			int BoundType = -1; - ��� ��
			int BoundInfo = -1; - ����� ����� 
		*/
		BoundS Bound[BoundCount]; // ��� ���������� � �������� �������� 

		Point e[FinitElementSize];
		int GlobalIdx[FinitElementSize];
	};

	struct GridS
	{
		int Dim = 0; // ����������� ����� � ���� � �� �� �����  
		int Nx = 0; // �������� ���������� ����� �� ��� �
		int Nz = 0; // �������� ���������� ����� �� ��� Z
		int Ny = 0; // �������� ���������� ����� �� ��� � ��� �������� ���������� ����� ��������� ����� 
		Point* Grid = nullptr; // ������ ����� ������������ ��� ��������� �������� ��������� 
	};

	/* �������� ������ �� ����� 
		@param
		std::string filename - ��� ����� ��������� ��������� �����  
		BaseGrid &baseGrid - ������� ����� �� ������� ����� �������������� ��������� ��������� ������� 
		@ret: void
	*/
	void Load(std::string filename, BaseGrid &baseGrid);

	/*  ��������� ��������� �������
		@param
		const BaseGrid &baseGrid - ������� ����� �������� ��������� ������� 7
		Grid &grid - ��������������� ��������� ������� 
		ret: void
	*/
	void GenerateGrid(BaseGrid &baseGrid, GridS &grid);


	/* ������������� �������, ���� ��� ��� ���� ������
		@param
		const BaseGrid &baseGrid - ������� ����� �������� ��������� ������� 7
		Grid &grid - ��������������� ��������� �������
		ret: void
		!!! ���������� ������� ������, ���� ��� ����� �� �������� ����� ������ 
	*/
	void ReGenerateGrid(BaseGrid& baseGrid, GridS& grid);


	/* ���������� �������� ������� �� ��� ������ 
		@param 
		const  GridS& grid - ����� 
		int idx
		ret Element - ��������� �� 8-� ����� � ��������� ����������. ������ ����� �������� �������������� ���������� � ��� ����� ������� ����������� 
	*/
	Element GetFinitElement(const GridS& grid,int idx);

	/* ��������� ����� �� �������� ����������� 
	* !!! �������� ������ �� ��������� ������� 
		@param
		BaseGrid& baseGrid - ������� ����� (� ���� ���� ����� �������� ���������� ������������� 
		int coef
		ret void
		��������� ������������� ����� � �������� ��������� ��� 
	*/
	void DivideGrid(BaseGrid& baseGrid, int coef);


	/* ������� ������ ��������� �������� �����
		@param
		GridS& grid - �����
		ret void
	*/
	void DeallocateGrid(GridS& grid);


	/* ����������� ������ */
	void DeallocateBaseGrid(BaseGrid& baseGrid);

	/* Debug function */
	void PrintListGrid(const GridS& grid, int list);

	void PrintElement(const Element& element);

	/********************************************************************************************************/

	/* ��������� ������ � ��� ���������� ����������� � �������� �� ������ ������ */

	namespace Grid_private
	{
		/* ����������� ����� ����� ����� ������� ��������� � �� ����� */
		/*
			@param
			const BaseGrid &baseGrid - 
			Grid &grid - 
			ret void 
			��������� � grid ���� Nx,Ny,Nz ���������������� ������ ���������� 
		*/
		void GetTotalNumberOfNodes(const BaseGrid &baseGrid, GridS &grid);
		
		/* ��������� ���� ��������� ������� ��� ������ ��������� ��������� � �������������� � ����� ���� ������� � ��������� ������� */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			��������� � grid ���� Grid ������� ���� �������� ����� 
		*/
		void GenerateBaseGrid(const BaseGrid& baseGrid, GridS& grid);

		/* ������� ������ �������������� � ������� � ������ ����� ��������� ��� ��� ���� */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			��������� � grid ���� Points.Info  ������� ���������� � ��� ����������� �� ����� �����-���� �������  + ����� ������ ��� ������ ������� 
		*/
		void DivisionIntoSubAreas(const BaseGrid& baseGrid, GridS& grid);

		/* ������� ������ ���� �� � ��������� ����� �������� �� ������� ��������� */
		/*
			@param
			const BaseGrid &baseGrid -
			Grid &grid -
			ret void
			��������� � grid ���� Points.Info  ������� ���������� � ��� ����������� �� ����� �����-���� ������� + ����� ������ ��� ������ �������
		*/
		void DivisionIntoSubBounds(const BaseGrid& baseGrid, GridS& grid);


	};

};