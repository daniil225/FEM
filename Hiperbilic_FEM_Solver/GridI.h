#ifndef GRID_H_
#define GRID_H_

#include <string>
#include <fstream>
#include "GridStatus.h"

using namespace std;

template <class BaseGridXD, class ElementXD>
class GridI
{

private:
protected:
    double eps = 1e-10; //  машинный ноль
    ofstream fout;      // Файловый поток на запись
    ifstream fin;       // Файловый поток на чтение

public:
    GridI() = default;

    /*
    @param const string &filename - Текстовый файл с разметкой
    @return GridStatus\\
    @result Загрузка базовой сетки из файла
    @warning Валидация входных данных не предусмотрена
    @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus Load(const string &filename) noexcept = 0;

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
    @return BaseGridXD
    @result -
    */
    inline virtual BaseGridXD GetBaseGrid() const noexcept = 0;

    /*
    @param int idx - Индекс центральной точки в глобальной нумерации
    @return ElementXD - Структура содержащая всю необходимую информацию о конечном элементе
    @result Получить конечный элемент с полным его описанием (Границы, Область)
    */
    virtual ElementXD GetElement(int idx) const noexcept = 0;

    /*
    @param void
    @return GridStatus
    @result Перегенерация сетки при изменении ее параметров
    @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus ReGenerateGrid() noexcept = 0;
    virtual ~GridI() = default;
};

#endif