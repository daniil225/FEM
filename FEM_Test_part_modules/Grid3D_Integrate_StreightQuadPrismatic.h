#ifndef GRID3D_INTEGRATE_STREIGHTQUADPRISMATIC_H_
#define GRID3D_INTEGRATE_STREIGHTQUADPRISMATIC_H_


#include "GridIntegrateI.h"
#include <vector>
#include <tuple>

using namespace std;
/*
    Структура содержит описание расчетной области в данном случае это 8-и угольник
    # Поля:
    + double Xstart, Xend;
    + double Ystart, Yend;
    + double Zstart, Zend;
*/
struct BaseGrid3D_Integrate_StreightQuadPrismatic
{
    
};

/* Размерность области */
typedef tuple<uint64_t, uint64_t, uint64_t> Grid3D_Size;

class Grid3D_Integrate_StreightQuadPrismatic : public GridIntegrateI<BaseGrid3D_Integrate_StreightQuadPrismatic>
{
private:
    /* Размерность по координатам */
    uint64_t _Xsize = 0;
    uint64_t _Ysize = 0;
    uint64_t _Zsize = 0;

    /* Массивы по координатам */
    vector<double> _XGrid;
    vector<double> _YGrid;
    vector<double> _ZGrid;

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

    /* Get/Set methods */

    inline uint64_t Xsize() const { return _Xsize; }
    inline uint64_t Ysize() const { return _Ysize; }
    inline uint64_t Zsize() const { return _Zsize; }

    /* Возврат значений в соответсвующих массивах/ модификация запрещена */
    inline double XGrid(uint64_t idx) const {return _XGrid[idx]; }
    inline double YGrid(uint64_t idx) const {return _YGrid[idx]; }
    inline double ZGrid(uint64_t idx) const {return _ZGrid[idx]; }
    

    ~Grid3D_Integrate_StreightQuadPrismatic() = default;
};

#endif