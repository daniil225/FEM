#ifndef UNITTESTS_H_
#define UNITTESTS_H_

/* 
    Система предназначена для автоматического тестирования всей системы 
    
    Реализованные модули тестирования:
    [-] UnitTests_Grid3D_StreightQuadPrismatic: Тестирование получения корректного построения расчетной области 
    [-] UnitTests_BasisLinearFunctionOnStreightQuadPrismatic: Тестирование процедуры корректного построения базисных функций + операция получения значения функции 
    [-] UnitTests_PointInfo: Тестирование модуля установки параметров точки 
    [-] UnitTests_SlauSolvers: Процедура тестирования решателей для матрицы  
    [-] UnitTests_Matrix: Тестирование корректности работы с матрицыми в различных форматах хранения
    [-] UnitTests_GenerateFEMLocalMatrix: Тестирование процедуры генерации локальных конечно-элементных матрицы для различных базисов и геметрии
    [-] UnitTests_ElepticEquationSolver: Тестирование операции тестирования Элептического уравнения методом конечных элементов. Проверка собранной МКЭ Слау + сверка результата
    [-] UnitTests_HyperbolicEquationSolver: Тестирование операции решения гиперболического уравнения с различными схемами по времени. В работе используется 4-х слойная не явная

    !!! Нужно будет подумать каким сделать классы. Логично занести базу тестов и прогонть все это непосредственно 
        Проблема с тем, что объекты могуть быть различной природы и соответсвенно нужно брать и делать, что-то похожее на 
        шаблонный класс. Думаю лучше сделать что-то около статического. Хотя можно брать и внутри класса создавать объекты тестов 
        для соответсвующих типов объектов (их не так уж и много) и в соответствие с типом уже применять подходящий тест для прогона. 
        Это все конечно в планах. В будущем это позволит делать прогоны при написании каких-либо других программ с походим функционалом
        В любом случае например если это матрицы и решатель, то их удобно проверять набором базовых тестов

    !!! Полагаю в отдельном файле создать базу тестов и уже их ставить на прогон. Таким образом можно будет увеличивать базу тестов
 */

/* Интерфейс для тестирования */
class UnitTests
{
    private:
    protected:
    public:
};


class UnitTests_Grid3D_StreightQuadPrismatic
{
    private:
    protected:
    public:

};

class UnitTests_BasisLinearFunctionOnStreightQuadPrismatic
{

    private:
    protected:
    public:
};

class UnitTests_PointInfo
{

    private:
    protected:
    public:
};

class UnitTests_SlauSolvers
{

    private:
    protected:
    public:
};

class UnitTests_Matrix
{
    private:
    protected:
    public:

};

class UnitTests_GenerateFEMLocalMatrix
{
    private:
    protected:
    public:

};

class UnitTests_ElepticEquationSolver
{
    private:
    protected:
    public:

};

class UnitTests_ElepticEquationSolver
{
    private:
    protected:
    public:

};

class UnitTests_HyperbolicEquationSolver
{
    private:
    protected:
    public:

};


#endif