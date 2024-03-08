#ifndef GRID_INTEGRATE_I_H_
#define GRID_INTEGRATE_I_H_

#include "GridStatus.h"

/*
    Общий интерфес сеток для областей инетгрирования. Подразумевает построение исключительно равномерных сеток в 1D, 2D, 3D
    Определяем общий набор функциональных возможностей необходимых для выполнения операции интегрирования

    # Шаблонные параметры:
    BaseGridIntegrate - Структура содержащая базовую разметку для расчетной области 
        для 1D - 2 точки 
        для 2D - 4 точки интегрирование по 4-х уголнику 
        для 3D - 8 точек интегрирование по 8-и угольнику 
*/
template <typename BaseGridIntegrate>
class GridIntegrateI
{

protected:
    BaseGridIntegrate BaseGrid;
    double DivideCoef = 1.0; // Коэффициент дробления 

public:
    GridIntegrateI() = default;
    GridIntegrateI(const BaseGridIntegrate &_BaseGrid) : BaseGrid(_BaseGrid) {}

    /*
        @param void
        @return GridStatus
        @result Генерация сетки
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus GenerateGrid() noexcept = 0;

    /*
        @param const int coef - Коэффициент дробления
        @return GridStatus
        @result Дробление сетки в заданное количество раз
        @warning Производит исключительно установку новых параметров дробление. Для их применения нужно вызвать метод ReGenerateGrid()
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus DivideGrid(const int coef) noexcept = 0;

    /*
        @param void
        @return GridStatus
        @result Перегенерация сетки при изменении ее параметров
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus ReGenerateGrid() noexcept = 0;

    /*
        @param void
        @return BaseGridIntegrate 
    */
    virtual BaseGridIntegrate GetBaseGrid() const { return BaseGrid; }

    /*
        @param const BaseGridIntegrate& _BaseGrid
        @return void
    */
    virtual void SetBaseGrid(const BaseGridIntegrate& _BaseGrid) { BaseGrid = _BaseGrid; }

    virtual ~GridIntegrateI() = default;
};

#endif