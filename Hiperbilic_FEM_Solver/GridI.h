#ifndef GRID_H_
#define GRID_H_

#include <string>
#include <fstream>
#include "GridStatus.h"

using namespace std;

template<class BaseGridXD, class ElementXD>
class GridI
{

private:


protected:
    double eps = 1e-7; //  машинный ноль 
    ofstream fout; // Файловый поток на запись
    ifstream fin; // Файловый поток на чтение 

public:

    GridI() = default;
    
    /* Загрузка базовой сетки из файла и по ней строится уже все */
    [[nodiscard]]
    virtual GridStatus Load(const string& filename) noexcept = 0;

    /* Генерация сетки */
    [[nodiscard]]
    virtual GridStatus GenerateGrid() noexcept = 0;

    /* Дробление сетки в заданное количество раз */
    [[nodiscard]]
    virtual GridStatus DivideGrid(const int coef) noexcept = 0;

    /* Получить базовую сетку */
    virtual BaseGridXD GetBaseGrid() const noexcept = 0;

    /* Получить шаблонный элемент который необходим */
    virtual ElementXD GetElement(int idx) const noexcept = 0;

    /* Перегенерация сетки при изменении ее параметров */
    [[nodiscard]]
    virtual GridStatus ReGenerateGrid() noexcept = 0;
    virtual ~GridI() = default;

};

#endif