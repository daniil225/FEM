#pragma once

#include <functional>
#include <tuple>
#include <vector>

/* ��������� �� ����������� � ����� 
	��� ������ ������� ����������� ���� ��������� 
*/

typedef std::function < double(double, double, double)> Func_3D;
typedef std::vector<Func_3D> ArrayFunc_3D;
typedef std::pair<Func_3D, Func_3D> PairFunc_3D;
typedef std::vector<PairFunc_3D> ArrayPairFunc_3D;
typedef std::vector<std::pair<int,int>> VecPairInt;
/* ����������� ���������� � ����*/
struct ParamOfDE
{
	ArrayFunc_3D Lambda; // �������� ��������� Lambda
	ArrayFunc_3D Gamma; // �������� ��������� �����
	ArrayFunc_3D F; // ������ ����� 
	
	/* ���������� ������� �� ������ ������ ����������� ��� �� */

	// ��� ��������������� ��������� � ��� �������������� ������ ������ � ������ ������ � �������� 1-� �� 2-� �� 3-� �� 
	// ������ ��� ����� ������ �������� � ������ ������ � �������������� ������� � ��� �� ������� � ������� ��� ���������� ���� 
	VecPairInt Map;  

	ArrayFunc_3D FirstBoundCond;
	ArrayFunc_3D SecondBoundCond;
	ArrayPairFunc_3D ThirdBoundCond;
};

/* �������� ������� ������ 
* @param
* const ParamOfDE &param - ��������� ��
* int num - ����� ������
*/
Func_3D GetLambda(const ParamOfDE &param,int num);

/* �������� ������� �����
* @param
* const ParamOfDE &param - ��������� ��
* int num - ����� ������
*/
Func_3D GetGamma(const ParamOfDE& param, int num);

/* �������� ������� ������ ����� 
* @param
* const ParamOfDE &param - ��������� ��
* int num - ����� ������
*/
Func_3D GetF(const ParamOfDE& param, int num);

/* �������� ������� �� 1 � 2 ��� ����� 
* @param
* const ParamOfDE &param - ��������� ��
* int num - ����� ������
*/
Func_3D GetBoundCond12(const ParamOfDE& param, int num);


/* �������� ���� ������� 3-� ��
* @param
* const ParamOfDE &param - ��������� ��
* int num - ����� ������
*/
PairFunc_3D GetBoundCond3(const ParamOfDE& param, int num);

/* ������������� ������� �����
* @param
* void
* ret ParamOfDE - ��������� ��
*/
ParamOfDE InitTest1();
