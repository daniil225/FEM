#ifndef GRID3D_INTEGRATE_STREIGHTQUADPRISMATIC_H_
#define GRID3D_INTEGRATE_STREIGHTQUADPRISMATIC_H_


#include "GridIntegrateI.h"
#include <vector>
#include <tuple>

using namespace std;
/*
    Структура для построения области инетгрирования из прямых призм. Интегрирование можно производить по всей 
    расчетной области. Даже для области построенной при МКЭ аппрокимации
    # Поля:
    + int32_t CountOfDivision = 1 - Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д
     + struct PointXZS - Структура для базовой точки:
        * double x = 0.0;
        * double z = 0.0;

    + struct DivideParamS Параметры структуры:
        * int32_t num = 0 - количество интервалов на которое нужно разделить отрезок
        * double coef = 0 - Коэффициент растяжения или сжатия

    + int32_t Nx = 0 - количество узлов вдоль горизонатального направления
    + int32_t Ny = 0 - количество узлов вдоль вертикального направления
    + int32_t Nz = 0 - количество узлов вдоль оси z
    + vector<vector<PointXZS>> BaseGridXZ - Базовая сетка в плосткости XZ
    + vector<double> BaseGridY - Базовая сетка по Y
*/
struct BaseGrid3D_Integrate_StreightQuadPrismatic
{
    int32_t CountOfDivision = 1; // Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д

     /*
        Структура для базовой точки
    */
    struct PointXZS
    {
        double x = 0.0;
        double z = 0.0;

        /* Копирование/присваивание */
        PointXZS() = default;
        PointXZS(const PointXZS &) = default;
        PointXZS &operator=(const PointXZS &) = default;
        PointXZS(PointXZS &&) = default;
        PointXZS &operator=(PointXZS &&) = default;

        ~PointXZS() = default;
    };

    /*
        Параметры структуры:
        + int32_t num = 0 - количество интервалов на которое нужно разделить отрезок
        + double coef = 0 - Коэффициент растяжения или сжатия
    */
    struct DivideParamS
    {
        int32_t num = 0; // количество интервалов на которое нужно разделить отрезок
        double coef = 0; // Коэффициент растяжения или сжатия

        /* Копирование/присваиваение */
        DivideParamS() = default;
        DivideParamS(const DivideParamS &) = default;
        DivideParamS &operator=(const DivideParamS &) = default;
        DivideParamS(DivideParamS &&) = default;
        DivideParamS &operator=(DivideParamS &&) = default;

        ~DivideParamS() = default;
    };
    int32_t Nx = 0; // количество узлов вдоль горизонатального направления
    int32_t Ny = 0; // количество узлов вдоль вертикального направления
    int32_t Nz = 0; // количество узлов вдоль оси z

    vector<vector<PointXZS>> BaseGridXZ; //  Базовая сетка в плосткости XZ
    vector<double> BaseGridY;            //  Базовая сетка по Y


     /*
        3 массива
        DivideParam[0] - массив для разбиения по оси x размер Nx-1
        DivideParam[1] - массив для разбиения по оси z размер Nz-1
        DivideParam[2] - массив для разбиения по оси y размер Ny-1
    */
    vector<vector<DivideParamS>> DivideParam; // Массив для разбиения

    bool isReadyToUse = false;

    /************************************************************************/
    /* Constructors and operator= */
    BaseGrid3D_Integrate_StreightQuadPrismatic() = default;
    BaseGrid3D_Integrate_StreightQuadPrismatic(BaseGrid3D_Integrate_StreightQuadPrismatic &&) = default;            // Для возврата из функции с захватом параметра
    BaseGrid3D_Integrate_StreightQuadPrismatic &operator=(BaseGrid3D_Integrate_StreightQuadPrismatic &&) = default; // Аналогично присваивание с захватом ресурсов

   BaseGrid3D_Integrate_StreightQuadPrismatic(const BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_)
    {
        this->BaseGridXZ = baseGrid_.BaseGridXZ;
        this->BaseGridY = baseGrid_.BaseGridY;
        this->CountOfDivision = baseGrid_.CountOfDivision;
        this->DivideParam = baseGrid_.DivideParam;
        this->isReadyToUse = baseGrid_.isReadyToUse;
        this->Nx = baseGrid_.Nx;
        this->Ny = baseGrid_.Ny;
        this->Nz = baseGrid_.Nz;
    }

    BaseGrid3D_Integrate_StreightQuadPrismatic &operator=(const BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_)
    {
        this->BaseGridXZ = baseGrid_.BaseGridXZ;
        this->BaseGridY = baseGrid_.BaseGridY;
        this->CountOfDivision = baseGrid_.CountOfDivision;
        this->DivideParam = baseGrid_.DivideParam;
        this->isReadyToUse = baseGrid_.isReadyToUse;
        this->Nx = baseGrid_.Nx;
        this->Ny = baseGrid_.Ny;
        this->Nz = baseGrid_.Nz;
        return *this;
    }

    BaseGrid3D_Integrate_StreightQuadPrismatic(BaseGrid3D_Integrate_StreightQuadPrismatic &) = delete;
    BaseGrid3D_Integrate_StreightQuadPrismatic &operator=(BaseGrid3D_Integrate_StreightQuadPrismatic &) = delete;
    /************************************************************************/

    ~BaseGrid3D_Integrate_StreightQuadPrismatic() = default;
};

/* Точка в пространстве для сетки интегрирования */
struct IntegratePoint3D
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/* 4-х угольная прямая призма набор из 8-и точек*/
struct StreightQuadPrismatic { IntegratePoint3D points[8]; };

class Grid3D_Integrate_StreightQuadPrismatic : public GridIntegrateI<BaseGrid3D_Integrate_StreightQuadPrismatic>
{
private:

  /*         все величины получаются после генерации сетки       */
    int32_t Dim = 0;      // Размерность сетки(количество узлов) и СЛАУ в то же время
    int32_t ElemCount = 0; // Количество полученных элементов 
    int32_t GlobalNx = 0; // Сумарное количество узлов по оси Х
    int32_t GlobalNy = 0; // Сумарное количество узлов по оси У
    int32_t GlobalNz = 0; // Сумарное количество узлов по оси Z

     /* Базовая сетка */
     BaseGrid3D_Integrate_StreightQuadPrismatic baseGrid;
    
     // Массив точек получающийся при генерации
    vector<IntegratePoint3D> Grid;

    /**************************************************************/

    /* Private Method */
    /*
        @param void
        @return void
        @note Dim - будет равняться размерности и матрицы СЛАУ и МКЭ сетки Расчитать общее число узлов получающееся в сетке
    */
    void GetTotalNumberOfNodes() noexcept;

    /*
        @param void
        @return void
        @note Генерация всей расчетной области без учетка фиктивных элементов и принадлежности к какой либо границе и расчетной области

    */
    void GenerateBaseGrid() noexcept;



public:
    Grid3D_Integrate_StreightQuadPrismatic() = default;
    Grid3D_Integrate_StreightQuadPrismatic(const BaseGrid3D_Integrate_StreightQuadPrismatic &_BaseGrid) : GridIntegrateI(_BaseGrid) {}

    /*
    @param void
    @return GridStatus
    @result Генерация сетки
    @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus GenerateGrid() noexcept;

    /*
        @param const int coef - Коэффициент дробления
        @return GridStatus
        @result Дробление сетки в заданное количество раз
        @warning Производит исключительно установку новых параметров дробление. Для их применения нужно вызвать метод ReGenerateGrid()
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus DivideGrid(const int coef) noexcept;

    /*
        @param void
        @return GridStatus
        @result Перегенерация сетки при изменении ее параметров
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus ReGenerateGrid() noexcept;

    /* Гетеры и сеттеры */
    /*
        @param void
        @return BaseGrid3DStreightQuadPrismatic\\
        @result -
    */
    inline BaseGrid3D_Integrate_StreightQuadPrismatic GetBaseGrid() const noexcept { return baseGrid; }

    /*
        @param int idx - Индекс центральной точки в глобальной нумерации
        @return Finit_Element_StreightQuadPrismatic - Структура содержащая всю необходимую информацию о конечном элементе
        @result Получить конечный элемент с полным его описанием (Границы, Область)
    */
    StreightQuadPrismatic GetElement(const int32_t idx) const noexcept;

    /*
        @param const BaseGrid3DStreightQuadPrismatic& baseGrid_ - Базовая сетка области
        @return GridStatus\\
        @result Установка параметров базовой сетки
        @warning Не изменяет уже построенной сетки. Для применения изменений нужно вызвать метод ReGenerateGrid()
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] GridStatus SetBaseGrid(const BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_) noexcept;

    /*
        @param int idx - индекс точки в МКЭ сетке
        @return Point& - точка в МКЭ Области
        @result -
    */
    inline IntegratePoint3D &operator[](const int32_t idx) noexcept { return Grid[static_cast<uint64_t>(idx)]; }

    

    ~Grid3D_Integrate_StreightQuadPrismatic() = default;
};

#endif