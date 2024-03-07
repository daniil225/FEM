#include "Grid3D_StreightQuadPrismatic.h"

#include <iostream>
#include <cmath>
#include <functional>

/* C lib */
#include <stdlib.h>

/* Private section */
void Grid3D_StreightQuadPrismatic::GetTotalNumberOfNodes() noexcept
{
    for (int32_t i = 0; i < baseGrid.Nx - 1; i++)
        GlobalNx += baseGrid.DivideParam[0][static_cast<uint64_t>(i)].num;

    for (int32_t i = 0; i < baseGrid.Nz - 1; i++)
        GlobalNz += baseGrid.DivideParam[1][static_cast<uint64_t>(i)].num;

    for (int32_t i = 0; i < baseGrid.Ny - 1; i++)
        GlobalNy += baseGrid.DivideParam[2][static_cast<uint64_t>(i)].num;

    FEMCount = GlobalNx * GlobalNy * GlobalNz;
    GlobalNx++;
    GlobalNy++;
    GlobalNz++;
    Dim = GlobalNx * GlobalNy * GlobalNz;
}

void Grid3D_StreightQuadPrismatic::GenerateBaseGrid() noexcept
{

    /* Вспомогальельная структура для определения параметров разбиения */

    struct SettingForDivide
    {
        double step = 0; // Шаг на отрезке
        double coef = 0; // Коэффициент увеличения шага
        int32_t num = 0; // Количество интервалов идем то num-1 и потом явно вставляем элемент

        /* Копирование и присваивание */
        SettingForDivide() = default;
        SettingForDivide(const SettingForDivide &) = default;
        SettingForDivide &operator=(const SettingForDivide &) = default;
        SettingForDivide(SettingForDivide &&) = default;
        SettingForDivide &operator=(SettingForDivide &&) = default;
    };
    /*******************************************************************/

    /* Вспомогательные функции для генерации сетки опеределены в виде лямбда функций */

    /*
        @param:
            int32_t i - Номер массива от 0 до 2
            int32_t j - Номер элемента в массиве
            double left - левая грани отрезка
            double right - правая граница отрезка
        @return: SettingForDivide -  структура с вычесленными параметрами деления сетки
        @result: Расчитываем шаг для сетки
    */
    std::function<SettingForDivide(int32_t, int32_t, double, double)> CalcSettingForDivide =
        [&](int32_t i, int32_t j, double left, double right) -> SettingForDivide
    {
        SettingForDivide res;
        int32_t num = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].num;
        double coef = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].coef;

        if (std::abs(coef - 1.0) > eps)
        {
            double coefStep = 1.0 + (coef * (std::pow(coef, num - 1) - 1.0)) / (coef - 1.0);

            res.step = (right - left) / coefStep;
        }
        else
        {
            res.step = (right - left) / num;
        }

        //  Убираем погрешность
        if (std::abs(res.step) < eps)
            res.step = 0.0;

        res.num = num;
        res.coef = coef;
        return res;
    };

    /*
        @param:
            SettingForDivide &param - параметр разбиения
            double left - левая граница отрезка
            double right - правая граница отрезка
            vector<double>& Line - генерируемый массив
            int &idx - индекс в массиве на какую позицию ставить элемент
        @return: void
        @result: Генерация разбиения по X или Z( гороизонтальная линия или вертикальная ) с учетом разбиения
    */
    std::function<void(SettingForDivide &, double, double, vector<double> &, int &idx)> GenerateDivide =
        [](SettingForDivide &param, double left, double right, vector<double> &Line, int &idx) -> void
    {
        int32_t num = param.num;
        double coef = param.coef;
        double step = param.step;

        Line[static_cast<uint64_t>(idx)] = left;
        idx++;
        double ak = left;
        for (int32_t k = 0; k < num - 1; k++)
        {
            ak = ak + step * std::pow(coef, k);
            Line[static_cast<uint64_t>(idx)] = ak;
            idx++;
        }
        Line[static_cast<uint64_t>(idx)] = right;
    };

    try
    {
        /* Подготовка переменных для генерации сетки */
        int32_t total = GlobalNx * GlobalNy * GlobalNz;
        Grid.resize(static_cast<uint64_t>(total));

        /* Псевдонимы поддержка локализации */
        int32_t Nx = baseGrid.Nx;
        int32_t Ny = baseGrid.Ny;
        int32_t Nz = baseGrid.Nz;
        vector<vector<BaseGrid3DStreightQuadPrismatic::PointXZS>> &BaseGridXZ = baseGrid.BaseGridXZ;
        vector<double> &BaseGridY = baseGrid.BaseGridY;

        /* Глобальные размеры есть */

        /* разбиения покоординатные */
        vector<double> LineX(static_cast<uint64_t>(GlobalNx)); // Массив элементов в строке по Х
        vector<double> LineZ(static_cast<uint64_t>(GlobalNz)); // Массив элементов в строке по Z
        vector<double> LineY(static_cast<uint64_t>(GlobalNy)); // Массив элементов в строке по Y
                                                               /**********************************************************************/

        /* Блок генерации разбиения */

        /*  Описание основной идеи генерации разбиения
            1) Сгенерируем одну плоскость xz, а потом ее растиражируем по y разбиение по x и z
            2) Расстановка элементов основных линий с учетом их расположения в сетке ( Опорная сетка ):
                Расстановка элементов по x ( Опорных )
                Нужно значения х расставить в соответсвующие строки они соответсвуют разбиению по z
        */

        for (int32_t i = 0; i < Nz; i++)
        {
            int32_t idx = 0;

            for (int32_t j = 0; j < Nx - 1; j++)
            {
                double left = BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].x;
                double right = BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j + 1)].x;
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Заносим соответствующие значения x на свои позиции */

            int32_t startIdx = Getlevel(i, 1) * GlobalNx;
            int32_t endIdx = startIdx + GlobalNx;
            for (int32_t k = startIdx, kk = 0; k < endIdx; k++, kk++)
                Grid[static_cast<uint64_t>(k)].x = LineX[static_cast<uint64_t>(kk)];
        }

        /* Расстановка элементов по z ( Опорных ) */
        for (int32_t i = 0; i < Nx; i++)
        {
            int32_t idx = 0;
            for (int32_t j = 0; j < Nz - 1; j++)
            {
                double left = BaseGridXZ[static_cast<uint64_t>(j)][static_cast<uint64_t>(i)].z;
                double right = BaseGridXZ[static_cast<uint64_t>(j + 1)][static_cast<uint64_t>(i)].z;
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineZ, idx);
            }

            /* Процедура расстановки узлов в глобальный массив */
            int32_t startIdx = Getlevel(i, 0); // Стартовый индекс для прохода по массиву
            for (int32_t k = 0; k < GlobalNz; k++)
            {
                // Скачки будут ровно на величину GlobalNx
                Grid[static_cast<uint64_t>(startIdx)].z = LineZ[static_cast<uint64_t>(k)];
                startIdx += GlobalNx;
            }
        }

        /************************************************************************/
        /*
                    Кратко Алгоритм:
                    Работаем с осью Z соответственно индексация будет происходить по этой оси
                    в цикле идем по всем столбцам массива сетки

                    Нужно получить левую и правую границу на каждом интервале
                    Сформировать массив отрезков по  данной координате
                    Занести полученный массив в Глобальную сетку
                */
        /* Цикл по всем горизонтальным линиям */
        for (int32_t i = 0; i < GlobalNx; i++)
        {
            int32_t idx = 0;
            /* Цикл по интеравалам оси Z */
            for (int32_t j = 0; j < Nz - 1; j++)
            {
                int32_t startIdx = i + GlobalNx * Getlevel(j, 1);
                int32_t endIdx = i + GlobalNx * Getlevel(j + 1, 1);
                double left = Grid[static_cast<uint64_t>(startIdx)].x; // Левая граница по х
                double right = Grid[static_cast<uint64_t>(endIdx)].x;  // Правая граница по х

                // Разбиение интервала подчиняется разбиению по оси z
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineZ, idx);
            }
            /* Занесение результата в Итоговый массив */

            int32_t startIdx = i; // Стартовая позиция
            for (int32_t k = 0; k < GlobalNz; k++)
            {
                Grid[static_cast<uint64_t>(startIdx)].x = LineZ[static_cast<uint64_t>(k)];
                startIdx += GlobalNx;
            }
        }

        /*
            Генерация вспомогательных горизонтальных линий
            Цикл по всем горизонтальным линиям
        */
        for (int32_t i = 0; i < GlobalNz; i++)
        {
            int32_t idx = 0;
            for (int32_t j = 0; j < Nx - 1; j++)
            {
                int32_t startIdx = Getlevel(j, 0) + i * GlobalNx;
                int32_t endIdx = Getlevel(j + 1, 0) + i * GlobalNx;
                double left = Grid[static_cast<uint64_t>(startIdx)].z;
                double right = Grid[static_cast<uint64_t>(endIdx)].z;
                // Разбиение интервала подчиняется разбиению по оси x
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Занесение результатов в Глобальную сетку */
            int32_t startIdx = i * GlobalNx;
            for (int32_t k = 0; k < GlobalNx; k++)
            {
                Grid[static_cast<uint64_t>(startIdx)].z = LineX[static_cast<uint64_t>(k)];
                startIdx++;
            }
        }
        /************************************************************************/

        /*
            Разбиение по y теражируем сечение с учетом сечения плоскостью
            Создадим разбиение элементов массива Y
        */

        int32_t idxY = 0;
        for (int32_t i = 0; i < Ny - 1; i++)
        {
            double left = BaseGridY[static_cast<uint64_t>(i)];
            double right = BaseGridY[static_cast<uint64_t>(i + 1)];
            SettingForDivide param = CalcSettingForDivide(2, i, left, right);
            GenerateDivide(param, left, right, LineY, idxY);
        }

        // Вставляем элемент Y0 в базовом сечении
        int32_t idx = 0;
        for (int32_t j = 0; j < GlobalNz; j++)
        {
            for (int32_t k = 0; k < GlobalNx; k++)
            {
                // // Инициализация битовых полей без этой команды будет Ошибка !!! ( В актуальной версии исправлен этот недочет )
                // PointInfo::ClearInfo(Grid[idx].info);
                Grid[static_cast<uint64_t>(idx)].y = LineY[0];
                idx++;
            }
        }

        // Тиражирование сетки по всем узлам
        // i отвечает за индексацию в массиве по Y
        // Idx за индексацию по элементам

        for (int32_t i = 1; i < GlobalNy; i++)
        {
            int32_t ListIdx = 0;
            for (int32_t j = 0; j < GlobalNz; j++)
            {
                for (int32_t k = 0; k < GlobalNx; k++)
                {
                    // Инициализация битовых полей без этой команды будет Ошибка !!! ( В актуальной версии исправлен этот недочет )
                    // PointInfo::ClearInfo(Grid[idx].info);
                    Grid[static_cast<uint64_t>(idx)].y = LineY[static_cast<uint64_t>(i)];
                    Grid[static_cast<uint64_t>(idx)].x = Grid[static_cast<uint64_t>(ListIdx)].x;
                    Grid[static_cast<uint64_t>(idx)].z = Grid[static_cast<uint64_t>(ListIdx)].z;
                    ListIdx++;
                    idx++;
                }
            }
        }
        /************************************************************************/
    }
    catch (const std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        // exit(State::MEMORY_ALLOC_ERROR);
    }
    catch (...)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, "Unknown error in Grid3D_StreightQuadPrismatic::GenerateBaseGrid\n");
    }
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubAreas() noexcept
{

    /* Вспомогальельная структура для определения параметров границы */
    struct BoundArea
    {
        int32_t AreaNum = -1;    // номер подобласти (определяет набор формул отвечающих параметрам ДУ)
        int32_t PlaneXZSize = 0; // Размер массива PlaneXZ
        vector<int32_t> PlaneXZ; // Массив точек многоугольника (номера)

        double StartY = 0;
        double EndY = 0;
    };

    /************************************************************************/
    try
    {
        /*
            Размер массива берем с запасом его размер равен Nx*Nz
            контур номеров точек многоугольника соттветствующих каким то координатам
        */
        vector<int32_t> BoundGrid(static_cast<uint64_t>(GlobalNx * GlobalNz));

        /* Вспомогательные фунции */

        /*
            @param:
                int32_t i - Номер подобласти (Порядковый) в массиве определяется порядок
            @return: BoundArea - сформированный массив подобласти
            @result: Получает подобласть в соответсвии с ее номером
        */
        std::function<BoundArea(int32_t)> GetBound = [&](int32_t i) -> BoundArea
        {
            BoundArea Bound;

            Bound.AreaNum = baseGrid.CalculationArea[static_cast<uint64_t>(i)][0]; // Выставили номер подобласти

            /* Вычислим все номера принадлежащие данной области */
            /* Инедексы границ многоугольника в глобальной нумерации */
            /* XZ */

            int32_t leftStartX = Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][3], 1);
            // int32_t leftEndX = Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][1], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][4], 1);
            int32_t rightStartX = Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][3], 1);
            int32_t rightEndX = Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.CalculationArea[static_cast<uint64_t>(i)][4], 1);

            int32_t nX = (rightEndX - rightStartX) / GlobalNx;
            int32_t nZ = (rightStartX - leftStartX);

            int32_t idxBoundGrid = 0;
            for (int32_t k = 0; k <= nX; k++)
            {
                int32_t Idx = leftStartX + k * GlobalNx;
                for (int32_t j = 0; j <= nZ; j++)
                {
                    BoundGrid[static_cast<uint64_t>(idxBoundGrid)] = Idx;
                    idxBoundGrid++;
                    Idx++;
                }
            }

            /* Y */
            /* В силу линейности по оси можно взять только первое и последнее значение */
            Bound.StartY = baseGrid.BaseGridY[static_cast<uint64_t>(baseGrid.CalculationArea[static_cast<uint64_t>(i)][5])];
            Bound.EndY = baseGrid.BaseGridY[static_cast<uint64_t>(baseGrid.CalculationArea[static_cast<uint64_t>(i)][6])];
            Bound.PlaneXZSize = idxBoundGrid;
            Bound.PlaneXZ = BoundGrid;

            return Bound;
        };

        /*
            @param
                const BoundArea& Bound - Массив границ
                int32_t numPointGlobal - Значение для поиска
            @return: bool - true - элемент в массиве есть false в противном случае
            @details: Бинарный поиск по массиву
        */
        std::function<bool(const BoundArea &, int32_t)> BinarySerch =
            [](const BoundArea &Bound, int32_t numPointGlobal) -> bool
        {
            const vector<int32_t> &arr = Bound.PlaneXZ;
            int32_t left = 0;
            int32_t right = Bound.PlaneXZSize - 1;
            int32_t midd = 0;
            while (1)
            {
                midd = (left + right) / 2;

                if (numPointGlobal < arr[static_cast<uint64_t>(midd)])      // если искомое меньше значения в ячейке
                    right = midd - 1;                                       // смещаем правую границу поиска
                else if (numPointGlobal > arr[static_cast<uint64_t>(midd)]) // если искомое больше значения в ячейке
                    left = midd + 1;                                        // смещаем левую границу поиска
                else                                                        // иначе (значения равны)
                    break;                                                  // функция возвращает индекс ячейки

                // если границы сомкнулись
                if (left > right)
                {
                    midd = -1;
                    break;
                }
            }
            if (midd == -1)
                return false;
            else
                return true;
        };

        /*
            @param
                const BoundArea &Bound - Заданная граница
                Point &point - Точка для проверки
                int numPointGlobal - индекс точки в глобальной нумерации
            @return: bool - true - Точка принадлежит заданной области, false в противном случае
            @result: Проверяет принадлежит ли точка Заданной области
        */
        auto IsInArea = [&](const BoundArea &Bound, Point &point, int numPointGlobal) -> bool
        {
            std::function<bool(double, double)> le =
                [&](double x1, double x2) -> bool
            {
                if (x1 < x2)
                    return true;
                if (std::abs(x1 - x2) < eps)
                    return true;
                return false;
            };

            std::function<bool(double, double)> ge =
                [&](double x1, double x2) -> bool
            {
                return le(x2, x1);
            };

            numPointGlobal %= GlobalNx * GlobalNz; // Проекция точки на плоскость
            bool arg1 = BinarySerch(Bound, numPointGlobal);
            bool arg2 = false;
            if (le(point.z, Bound.EndY) && ge(point.z, Bound.StartY))
                arg2 = true;

            return arg1 && arg2;
        };
        /************************************************************************/

        /* В цикле по всем областям */
        for (int32_t i = 0; i < baseGrid.L; i++)
        {
            BoundArea Bound = GetBound(i);

            /* По той части области где распологаются элементы (включая фиктивный)*/
            /* Для этого получим минимальный и максимальный значений индексов */
            // int leftAreaBound = ;
            // int RightAreaBound = ;

            /* Этот цикл тупой он по всем элементам идет */
            for (int32_t j = 0; j < Dim; j++)
            {
                if (IsInArea(Bound, Grid[static_cast<uint64_t>(j)], j))
                {
                    /* Мы в заданной области устанавливаем нужные параметры */
                    InfoManeger::SetFictitious(Grid[static_cast<uint64_t>(j)].info, 1);
                    InfoManeger::SetAreaInfo(Grid[static_cast<uint64_t>(j)].info, static_cast<uint32_t>(Bound.AreaNum));
                }
            }
        }

        /************************************************************************/
    }
    catch (const std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        // exit(State::MEMORY_ALLOC_ERROR);
    }
    catch (...)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, "Unknown error in Grid3D_StreightQuadPrismatic::DivisionIntoSubAreas\n");
    }
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubBounds() noexcept
{
    /* Вспомогальельная структура для определения параметров границы */
    struct Bound
    {
        int32_t Size = 0;
        vector<int32_t> Plane; // Массив точек границы
        int32_t BoundType = 0;
        int32_t BoundFormula = -1;
    };

    /* Вспомогательные функции */

    std::function<std::pair<int32_t, int32_t>(int32_t, int32_t, int32_t)> GetMaxPair =
        [](int32_t a, int32_t b, int32_t c) -> std::pair<int, int>
    {
        int32_t arr[3]{a, b, c};
        std::sort(arr, arr + 3);
        return std::make_pair(arr[1], arr[2]);
    };

    /*
        @param:
            const BoundArea& Bound - Массив границ
            int32_t numPointGlobal - Значение для поиска
        @return: bool - true - элемент в массиве есть false в противном случае
        @details: Бинарный поиск
    */
    std::function<bool(const Bound &, int32_t)> BinarySerch =
        [](const Bound &bound, int32_t numPointGlobal) -> bool
    {
        const vector<int32_t> &arr = bound.Plane;
        int32_t left = 0;
        int32_t right = bound.Size - 1;
        int32_t midd = 0;
        while (1)
        {
            midd = (left + right) / 2;

            if (numPointGlobal < arr[static_cast<uint64_t>(midd)])      // если искомое меньше значения в ячейке
                right = midd - 1;                                       // смещаем правую границу поиска
            else if (numPointGlobal > arr[static_cast<uint64_t>(midd)]) // если искомое больше значения в ячейке
                left = midd + 1;                                        // смещаем левую границу поиска
            else                                                        // иначе (значения равны)
                break;                                                  // функция возвращает индекс ячейки

            // если границы сомкнулись
            if (left > right)
            {
                midd = -1;
                break;
            }
        }
        if (midd == -1)
            return false;
        else
            return true;
    };

    std::function<bool(const Bound &, int32_t)> IsInBound =
        [&](const Bound &bound, int32_t numPointGlobal) -> bool
    { return BinarySerch(bound, numPointGlobal); };

    try
    {
        std::pair<int, int> max2 = GetMaxPair(GlobalNx, GlobalNz, GlobalNy);
        vector<int32_t> MemoryPool(static_cast<uint64_t>(max2.first * max2.second));
        /* Краевые условия задаются практически так же как и области с тем лишь исключением, что одна из координат фиксируется */
        std::function<Bound(int32_t)> GetBound = [&](int32_t i) -> Bound
        {
            Bound bound;
            bound.BoundType = baseGrid.BoundsArea[static_cast<uint64_t>(i)][1];
            bound.BoundFormula = baseGrid.BoundsArea[static_cast<uint64_t>(i)][0];

            /* Это определяется равенством  соответствующих точек*/
            /* y - фиксирован */
            if (baseGrid.BoundsArea[static_cast<uint64_t>(i)][6] == baseGrid.BoundsArea[static_cast<uint64_t>(i)][7])
            {
                /*  Определим смещение для оси y */
                int32_t startBaseIdx = GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][6], 2); // Узнаем край для

                /* Определяем границы областей */
                int32_t leftStartX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][4], 1);
                // int32_t leftEndX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][5], 1);
                int32_t rightStartX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][4], 1);
                int32_t rightEndX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][5], 1);

                /* Сколько шагать  */
                int32_t nX = (rightEndX - rightStartX) / GlobalNx;
                int32_t nZ = (rightStartX - leftStartX);

                int32_t idxBoundGrid = 0;
                for (int32_t k = 0; k <= nX; k++)
                {
                    int32_t Idx = leftStartX + k * GlobalNx + startBaseIdx; // Смещенный индекс
                    for (int32_t j = 0; j <= nZ; j++)
                    {
                        MemoryPool[static_cast<uint64_t>(idxBoundGrid)] = Idx;
                        idxBoundGrid++;
                        Idx++;
                    }
                }

                bound.Plane = MemoryPool;
                bound.Size = idxBoundGrid;
            }
            /* х - фиксирован */
            else if (baseGrid.BoundsArea[static_cast<uint64_t>(i)][2] == baseGrid.BoundsArea[static_cast<uint64_t>(i)][3])
            {
                /* Определим базовые узлы на прямой OZ а потом растиражируем узлы по границе ( Дробление фиксированное )  шаг по массиву GlobalNx*/
                int32_t StartPozitionZ = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][4], 1); // Стартовая позиция для  Z
                int32_t EndPozitionZ = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][5], 1);   // Конечная точка для оси Z

                int32_t StartPositionY = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][6], 2); // Стартовая позиция для Y
                int32_t EndPositionY = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][7], 2);   // Конечня позиция для Y

                int32_t nZ = (EndPozitionZ - StartPozitionZ) / GlobalNx;              //  Количество узлов по оси Z
                int32_t nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // Количество шагов для Y

                int32_t idxBoundGrid = 0;
                for (int32_t k = 0; k <= nY; k++)
                {
                    int32_t Idx = StartPositionY + k * GlobalNx * GlobalNz;
                    for (int32_t j = 0; j <= nZ; j++)
                    {
                        MemoryPool[static_cast<uint64_t>(idxBoundGrid)] = Idx;
                        idxBoundGrid++;
                        Idx += GlobalNx;
                    }
                }
                bound.Plane = MemoryPool;
                bound.Size = idxBoundGrid;
            }
            /* z - фиксирован */
            else if (baseGrid.BoundsArea[static_cast<uint64_t>(i)][4] == baseGrid.BoundsArea[static_cast<uint64_t>(i)][5])
            {
                /* Определяем базовые узлы по оси OX, а потом аналогично растиражируем узлы по границе (Дробление фиксированное) шаг по массиву + 1 */
                int32_t StartPositionX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][4], 1); // Стартовая позиция по оси Х
                int32_t EndPositionX = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][3], 0) + GlobalNx * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][4], 1);   // Конечная позиция по оси Х

                int32_t StartPositionY = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][6], 2); // Стартовая позиция для Y
                int32_t EndPositionY = Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][2], 0) + GlobalNx * GlobalNz * Getlevel(baseGrid.BoundsArea[static_cast<uint64_t>(i)][7], 2);   // Конечня позиция для Y правй конец (если смотреть со стороны положительного направления)

                int32_t nX = EndPositionX - StartPositionX;                           //  Количество узлов по оси X
                int32_t nY = (EndPositionY - StartPositionY) / (GlobalNx * GlobalNz); // Количество шагов для Y

                int32_t idxBoundGrid = 0;
                for (int32_t k = 0; k <= nY; k++)
                {
                    int32_t Idx = StartPositionY + StartPositionX + k * GlobalNx * GlobalNz;
                    for (int32_t j = 0; j <= nX; j++)
                    {
                        MemoryPool[static_cast<uint64_t>(idxBoundGrid)] = Idx;
                        idxBoundGrid++;
                        Idx++;
                    }
                }
                bound.Plane = MemoryPool;
                bound.Size = idxBoundGrid;
            }

            return bound;
        };

        /* Цикл по всем границам */
        for (int32_t i = 0; i < baseGrid.P; i++)
        {
            /* Получить список точек границы  */
            Bound bound = GetBound(i);
            for (int32_t j = 0; j < Dim; j++)
            {
                /* Расставляем нужные значения в соответствующие точки сетки */
                if (IsInBound(bound, j))
                {
                    /* Мы в заданной области устанавливаем нужные параметры */
                    InfoManeger::SetBoundInfo(Grid[static_cast<uint64_t>(j)].info, static_cast<uint32_t>(bound.BoundFormula), static_cast<uint32_t>(bound.BoundType + 1));
                }
            }
        }
    }
    catch (const std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        // exit(State::MEMORY_ALLOC_ERROR);
    }
    catch (...)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, "Unknown error in Grid3D_StreightQuadPrismatic::DivisionIntoSubBounds\n");
    }

    /************************************************************************/
}

int32_t Grid3D_StreightQuadPrismatic::Getlevel(int32_t i, int32_t axis) const noexcept
{
    int32_t res = 0;
    for (int32_t k = 0; k < i; k++)
        res += baseGrid.DivideParam[static_cast<uint64_t>(axis)][static_cast<uint64_t>(k)].num;
    return res;
}

/* Public section */

Grid3D_StreightQuadPrismatic::Grid3D_StreightQuadPrismatic(const string &filename)
{
    // Если загрузка не выполнена успешно то генерируем исключение с возвратом статуса
    if (Load(filename).GetState() != State::OK)
        throw Status;
}

Grid3D_StreightQuadPrismatic::Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)
{
    // Если структура не валидная то выбрасываем ошибку
    if (!baseGrid_.isReadyToUse)
    {
        Status.SetStatus(State::LOAD_ERROR, "Incorrect Base Grid. Call from Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)");
        throw Status;
    }
    else
        baseGrid = baseGrid_;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::Load(const string &filename) noexcept
{
    try
    {
        /* code */
        fin.open(filename);
        if (!fin.is_open())
        {
            Status.SetStatus(State::FILE_OPEN_ERROR, "Ошибка открытия файла " + filename + " в Grid2D_Quad::Grid2D_Quad(const string &filename)\n");
            fout.close();
            return Status;
        }
        /* Базовая сетка по XZ */
        fin >> baseGrid.Nx >> baseGrid.Nz;
        baseGrid.BaseGridXZ = vector(static_cast<uint64_t>(baseGrid.Nz), vector<BaseGrid3DStreightQuadPrismatic::PointXZS>(static_cast<uint64_t>(baseGrid.Nx)));
        for (int32_t i = 0; i < baseGrid.Nz; i++)
        {
            for (int32_t j = 0; j < baseGrid.Nx; j++)
            {
                fin >> baseGrid.BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].x >> baseGrid.BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].z;
            }
        }
        /***********************/

        /* Базовая сетка по Y */
        fin >> baseGrid.Ny;
        baseGrid.BaseGridY = vector<double>(static_cast<uint64_t>(baseGrid.Ny));
        for (int32_t i = 0; i < baseGrid.Ny; i++)
            fin >> baseGrid.BaseGridY[static_cast<uint64_t>(i)];
        /***********************/

        /* Расчетные подобласти */
        fin >> baseGrid.L;
        baseGrid.CalculationArea = vector(static_cast<uint64_t>(baseGrid.L), vector<int32_t>(static_cast<uint64_t>(baseGrid.SizeOfCalculationAreaElemet)));

        for (int32_t i = 0; i < baseGrid.L; i++)
        {
            for (int32_t j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
            {
                fin >> baseGrid.CalculationArea[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)];
                baseGrid.CalculationArea[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)]--; // приведение нумераии с нуля
            }
        }
        /***********************/

        /* Описание Границ */
        fin >> baseGrid.P;
        baseGrid.BoundsArea = vector(static_cast<uint64_t>(baseGrid.P), vector<int32_t>(static_cast<uint64_t>(baseGrid.SizeOfBoundsAreaElement)));
        for (int32_t i = 0; i < baseGrid.P; i++)
        {
            for (int32_t j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
            {
                fin >> baseGrid.BoundsArea[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)];
                baseGrid.BoundsArea[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)]--; // приведение нумераии с нуля
            }
        }
        /***********************/

        /* Правила дробления базовой сетки */
        baseGrid.DivideParam.resize(static_cast<uint64_t>(baseGrid.SizeOfDivideParam));
        baseGrid.DivideParam[0].resize(static_cast<uint64_t>(baseGrid.Nx - 1));
        baseGrid.DivideParam[1].resize(static_cast<uint64_t>(baseGrid.Nz - 1));
        baseGrid.DivideParam[2].resize(static_cast<uint64_t>(baseGrid.Ny - 1));

        for (int32_t i = 0; i < baseGrid.Nx - 1; i++)
            fin >> baseGrid.DivideParam[0][static_cast<uint64_t>(i)].num >> baseGrid.DivideParam[0][static_cast<uint64_t>(i)].coef;

        for (int32_t i = 0; i < baseGrid.Nz - 1; i++)
            fin >> baseGrid.DivideParam[1][static_cast<uint64_t>(i)].num >> baseGrid.DivideParam[1][static_cast<uint64_t>(i)].coef;

        for (int32_t i = 0; i < baseGrid.Ny - 1; i++)
            fin >> baseGrid.DivideParam[2][static_cast<uint64_t>(i)].num >> baseGrid.DivideParam[2][static_cast<uint64_t>(i)].coef;
        /***********************/
        baseGrid.isReadyToUse = true; // Структуру можно использовать она инициализированна
        fin.close();
    }
    catch (std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        return Status;
    }
    catch (std::exception &e)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, e.what());
        return Status;
    }
    catch (...)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, "Unknown Error in  Grid3D_StreightQuadPrismatic::Load");
        return Status;
    }

    return Status;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::GenerateGrid() noexcept
{
    // Расчет общего количества узлов в сетке
    // Для этого пройдемся по массиву разбиения каждого отрезка и вычислим общее число узлов
    GetTotalNumberOfNodes();

    // Генерация базовой сетки
    if (Status.GetState() == State::OK)
        GenerateBaseGrid();
    else
        return Status;
    // Учет фиктивных узлов
    if (Status.GetState() == State::OK)
        DivisionIntoSubAreas();
    else
        return Status;

    // Учет КУ и расстановка границ
    if (Status.GetState() == State::OK)
        DivisionIntoSubBounds();
    else
        return Status;

    return Status;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::DivideGrid(const int32_t coef) noexcept
{
    /* Для этого нужно сделать преобразование базовой сетки на заданный коэффициент и пересчитать всю сетку*/

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Nx - 1); i++)
    {
        baseGrid.DivideParam[0][i].num *= coef;
        baseGrid.DivideParam[0][i].coef = pow(baseGrid.DivideParam[0][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Nz - 1); i++)
    {
        baseGrid.DivideParam[1][i].num *= coef;
        baseGrid.DivideParam[1][i].coef = pow(baseGrid.DivideParam[1][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Ny - 1); i++)
    {
        baseGrid.DivideParam[2][i].num *= coef;
        baseGrid.DivideParam[2][i].coef = pow(baseGrid.DivideParam[2][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    return Status;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::ReGenerateGrid() noexcept
{
    /* Предварительно выставляем по дефолту все параметры класса */
    Dim = GlobalNx = GlobalNy = GlobalNz = 0;
    Grid.clear();

    /* Перегенерация сетки */
    Status = GenerateGrid();

    return Status;
}

Finit_Element_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const int32_t idx) const noexcept
{
    /* Оперделим список локальных функций вспомогательных */

    /*
        @param int32_t *arr
        @param const int32_t size
        @return Занууляет массив
    */
    auto ZeroInitArray = [](int32_t *arr, const int32_t size)
    {
        for (int32_t i = 0; i < size; i++)
            arr[static_cast<uint64_t>(i)] = 0;
    };

    /*
        @param std::pair<int32_t, int32_t>* arr
        @param const int32_t size
        @return Занууляет массив пар
    */
    auto ZeroInitArrayPair = [](std::pair<int32_t, int32_t> *arr, const int32_t size)
    {
        for (int32_t i = 0; i < size; i++)
        {
            arr[static_cast<uint64_t>(i)].first = 0;
            arr[static_cast<uint64_t>(i)].second = 0;
        }
    };

    /*
        @param int32_t* arr
        @param const int32_t size
        @return int32_t - Возвращает индекс массива где лежит max значение
    */
    auto GetArgMax = [](int32_t *arr, const int32_t size) -> int32_t
    {
        int32_t res = -1;
        for (int32_t i = 0; i < size; i++)
        {
            if (arr[static_cast<uint64_t>(i)] == 8)
            {
                res = i;
                break;
            }
        }

        return res;
    };

    /*
        @param std::pair<int32_t, int32_t> *arr
        @param const int32_t size
        @return std::pair<int32_t, int32_t> - хз что тут давно писал но работает
    */
    auto GetArgMaxPair = [](std::pair<int32_t, int32_t> *arr, const int32_t size) -> std::pair<int32_t, int32_t>
    {
        std::pair<int32_t, int32_t> res = std::make_pair(-1, -1);
        for (int32_t i = 0; i < size; i++)
        {
            if (arr[static_cast<uint64_t>(i)].first == 4)
            {
                res.first = i;
                res.second = arr[i].second;
                break;
            }
        }

        return res;
    };

    int32_t Nx = GlobalNx;
    int32_t Ny = GlobalNy;
    int32_t Nz = GlobalNz;

    /*
        @param int32_t idx - Индекс конечного элемента
        @return int32_t\\
        @return  Вернет по номеру конечного элемента значение его праой нижней границе (локальный номер 1)
    */
    auto K = [Nx, Ny, Nz](int32_t _idx) -> int32_t
    {
        int32_t shiftXZ = 0;
        int32_t projidx = _idx % (((Nx - 1) * (Nz - 1)));
        if (projidx < Nx - 1)
        {
            shiftXZ = projidx;
        }
        else
        {
            int32_t level = static_cast<int32_t>(std::floor(projidx / (Nx - 1)));
            int32_t start = level * Nx;
            int32_t shift = projidx - (Nx - 1) * level;
            shiftXZ = start + shift;
        }

        int32_t levelY = static_cast<int32_t>(std::floor((_idx) / ((Nx - 1) * (Nz - 1))));
        int32_t shiftY = levelY * Nx * Nz;
        return shiftY + shiftXZ;
    };

    Finit_Element_StreightQuadPrismatic FinitElement;
    /* Расчет глобальных индексов */
    FinitElement.GlobalIdx[0] = K(idx);
    FinitElement.GlobalIdx[1] = FinitElement.GlobalIdx[0] + 1;
    FinitElement.GlobalIdx[2] = FinitElement.GlobalIdx[0] + Nx;
    FinitElement.GlobalIdx[3] = FinitElement.GlobalIdx[1] + Nx;

    FinitElement.GlobalIdx[4] = FinitElement.GlobalIdx[0] + Nx * Nz;
    FinitElement.GlobalIdx[5] = FinitElement.GlobalIdx[4] + 1;
    FinitElement.GlobalIdx[6] = FinitElement.GlobalIdx[4] + Nx;
    FinitElement.GlobalIdx[7] = FinitElement.GlobalIdx[5] + Nx;

    for (int32_t i = 0; i < FinitElement.FinitElementSize; i++)
    {
        FinitElement.e[static_cast<uint64_t>(i)] = Grid[static_cast<uint64_t>(FinitElement.GlobalIdx[static_cast<uint64_t>(i)])];

        /* Проверим сразу на фиктивность  */
        if (!InfoManeger::IsFiFictitious(FinitElement.e[static_cast<uint64_t>(i)].info))
        {
            /* Элемент фиктивный устанавливаем этот факт и прекращаем цикл */
            FinitElement.isFictive = false; // Да фиктивный
            return FinitElement;            // Возврат из функции
        }
    }
    /* Опредеелим теперь всю информацию о конечном элементе */
    /* Если мы здесь значит элемент не фиктивный */
    FinitElement.isFictive = true; // НЕ ФИКТИВНЫЙ

    /* Идем по подобласти и определяем принадлежность к области */
    
    int32_t Map[Comand::RightBoundValAreaInfo]; // Массив всевозможных элементов подобласти
    ZeroInitArray(Map, Comand::RightBoundValAreaInfo);

    for (int32_t i = 0; i < FinitElement.FinitElementSize; i++)
    {
        AreaInfo Area = InfoManeger::GetAreaInfo(FinitElement.e[static_cast<uint64_t>(i)].info);
        /* Идем по структурке и расставляем все типы что есть  */
        for (int32_t j = 0; j < Area.size; j++)
            Map[Area.Cond[static_cast<uint64_t>(j)]]++;
    }
    /* В карте лежит информация о количестве всевозможных типах */
    /* Нам нужна тот индекс где лежит 8 это будет max*/
    FinitElement.AreaInfo = GetArgMax(Map, Comand::RightBoundValAreaInfo);

    /* Разбираемся с краевыми */
    /* Установим локальные номера */
    // {1, 2, 3, 4}, { 5,6,7,8 }, { 1,3,5,7 }, { 2,4,6,8 }, { 1,2,5,6 }, { 3,4,7,8 }
    FinitElement.Bound[0].LocalIdx[0] = 0;
    FinitElement.Bound[0].LocalIdx[1] = 1;
    FinitElement.Bound[0].LocalIdx[2] = 2;
    FinitElement.Bound[0].LocalIdx[3] = 3;

    FinitElement.Bound[1].LocalIdx[0] = 4;
    FinitElement.Bound[1].LocalIdx[1] = 5;
    FinitElement.Bound[1].LocalIdx[2] = 6;
    FinitElement.Bound[1].LocalIdx[3] = 7;

    FinitElement.Bound[2].LocalIdx[0] = 0;
    FinitElement.Bound[2].LocalIdx[1] = 2;
    FinitElement.Bound[2].LocalIdx[2] = 4;
    FinitElement.Bound[2].LocalIdx[3] = 6;

    FinitElement.Bound[3].LocalIdx[0] = 1;
    FinitElement.Bound[3].LocalIdx[1] = 3;
    FinitElement.Bound[3].LocalIdx[2] = 5;
    FinitElement.Bound[3].LocalIdx[3] = 7;

    FinitElement.Bound[4].LocalIdx[0] = 0;
    FinitElement.Bound[4].LocalIdx[1] = 1;
    FinitElement.Bound[4].LocalIdx[2] = 4;
    FinitElement.Bound[4].LocalIdx[3] = 5;

    FinitElement.Bound[5].LocalIdx[0] = 2;
    FinitElement.Bound[5].LocalIdx[1] = 3;
    FinitElement.Bound[5].LocalIdx[2] = 6;
    FinitElement.Bound[5].LocalIdx[3] = 7;

    /* Даем глобальные */
    for (int32_t i = 0; i <  FinitElement.BoundCount; i++)
    {
        for (int32_t j = 0; j < static_cast<int32_t>(Comand::RightBoundBoundCount); j++)
        {
            FinitElement.Bound[static_cast<uint64_t>(i)].GlobalIdx[static_cast<uint64_t>(j)] =  FinitElement.GlobalIdx[ FinitElement.Bound[static_cast<uint64_t>(i)].LocalIdx[static_cast<uint64_t>(j)]];
        }
    }

    /* По все границам идем и устанавливаем информацию о границе */
    for (int32_t i = 0; i < FinitElement.BoundCount; i++)
    {
        std::pair<int32_t, int32_t> MapB[Comand::RightBoundValBoundInfo];
        ZeroInitArrayPair(MapB, Comand::RightBoundValBoundInfo);

        for (int32_t j = 0; j < static_cast<int32_t>(Comand::RightBoundBoundCount); j++)
        {
            BoundInfo bnd = InfoManeger::GetBoundInfo(FinitElement.e[FinitElement.Bound[static_cast<uint64_t>(i)].LocalIdx[static_cast<uint64_t>(j)]].info); // Получить информацию о границах на границе

            /* Идем по структуре */
            for (int32_t k = 0; k < bnd.size; k++)
            {
                MapB[bnd.Cond[k]].first++;
                MapB[bnd.Cond[k]].second = bnd.TypeCond[k];
            }
        }
        /* Получаем максимальный аргумент и устанавливаем значения для границы */
        std::pair<int32_t, int32_t> Cond = GetArgMaxPair(MapB, Comand::RightBoundValBoundInfo);
        FinitElement.Bound[i].BoundInfo = Cond.first;
        FinitElement.Bound[i].BoundType = Cond.second;
    }

    /* Расстановка фиктивных границ */
    for (int32_t i = 0; i < FinitElement.BoundCount; i++)
    {
        if (FinitElement.Bound[i].BoundInfo != -1)
            FinitElement.Bound[i].IsBound = true;
    }

    return FinitElement;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::SetBaseGrid(const BaseGrid3DStreightQuadPrismatic &baseGrid_) noexcept
{
    // Если структура не валидная то выбрасываем ошибку
    if (!baseGrid_.isReadyToUse)
        Status.SetStatus(State::LOAD_ERROR, "Incorrect Base Grid. Call from Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)");
    else
        baseGrid = baseGrid_;
    return Status;
}

Finit_Element_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const double x, const double y, const double z) const noexcept
{
    Finit_Element_StreightQuadPrismatic FEMElement;
    /* Определим конечный элемент которому принадлежит точка и после запустип процедуру формирования КЭ по индексу */
    
    int32_t Nx = GlobalNx;
    int32_t Ny = GlobalNy;
    int32_t Nz = GlobalNz;
    auto K = [Nx, Ny, Nz](int32_t _idx) -> int32_t
    {
        int32_t shiftXZ = 0;
        int32_t projidx = _idx % (((Nx - 1) * (Nz - 1)));
        if (projidx < Nx - 1)
        {
            shiftXZ = projidx;
        }
        else
        {
            int32_t level = static_cast<int32_t>(std::floor(projidx / (Nx - 1)));
            int32_t start = level * Nx;
            int32_t shift = projidx - (Nx - 1) * level;
            shiftXZ = start + shift;
        }

        int32_t levelY = static_cast<int32_t>(std::floor((_idx) / ((Nx - 1) * (Nz - 1))));
        int32_t shiftY = levelY * Nx * Nz;
        return shiftY + shiftXZ;
    };

    
    /* Расчет глобальных индексов */
    int32_t idx = 0;
    for(idx = 0; idx < FEMCount; idx++)
    {
        Finit_Element_StreightQuadPrismatic FinitElement;
        FinitElement.GlobalIdx[0] = K(idx);
        FinitElement.GlobalIdx[1] = FinitElement.GlobalIdx[0] + 1;
        FinitElement.GlobalIdx[2] = FinitElement.GlobalIdx[0] + Nx;
        FinitElement.GlobalIdx[3] = FinitElement.GlobalIdx[1] + Nx;

        FinitElement.GlobalIdx[4] = FinitElement.GlobalIdx[0] + Nx * Nz;
        FinitElement.GlobalIdx[5] = FinitElement.GlobalIdx[4] + 1;
        FinitElement.GlobalIdx[6] = FinitElement.GlobalIdx[4] + Nx;
        FinitElement.GlobalIdx[7] = FinitElement.GlobalIdx[5] + Nx;

        for (int32_t i = 0; i < FinitElement.FinitElementSize; i++)
        {
            FinitElement.e[static_cast<uint64_t>(i)] = Grid[static_cast<uint64_t>(FinitElement.GlobalIdx[static_cast<uint64_t>(i)])];
        }

        double XStart = FinitElement.e[0].x;
        double XEnd = FinitElement.e[1].x;

        double ZStart = FinitElement.e[0].z;
        double ZEnd = FinitElement.e[2].z;

        double YStart = FinitElement.e[0].y;
        double YEnd = FinitElement.e[4].y;

        std::cout << "idx = " << idx << "\n";
        std::cout << "XStart = " << XStart << " XEnd = " << XEnd << "\n";
        std::cout << "ZStart = " << ZStart << " ZEnd = " << ZEnd << "\n";
        std::cout << "YStart = " << YStart << " YEnd = " << YEnd << "\n";
        std::cout << "\n";
        if((x <= XEnd && x >= XStart) && (y <= YEnd && y >= YStart) && (z <= ZEnd && z >= ZStart)) break;

    }
    std::cout << "idx = " << idx << "\n";
    return GetElement(idx);
}

/* Friend Functions + вспомогательные функции для некоторых нужд */
ostream &operator<<(ostream &os, Grid3D_Size &Grid3D_Size_param)
{
    os << "Grid3D_Size parametrs:\n";
    os << "Dim = " << Grid3D_Size_param.Dim << "\n";
    os << "FEMCount = " << Grid3D_Size_param.FEMCount << "\n";
    os << "GlobalNx = " << Grid3D_Size_param.GlobalNx << "\n";
    os << "GlobalNy = " << Grid3D_Size_param.GlobalNy << "\n";
    os << "GlobalNz = " << Grid3D_Size_param.GlobalNz << "\n";
    os << "\n";

    return os;
}

void Grid3D_StreightQuadPrismatic::PrintGridSlice(int32_t level) const
{
    uint64_t idx = static_cast<uint64_t>(level * GlobalNx * GlobalNz);
    cout << "Start idx: " << idx << "  End idx: " << idx + static_cast<uint64_t>(GlobalNx * GlobalNz - 1) << " Step Row: " << GlobalNx << "\n";
    cout << "Format = (x,z)\n";
    cout << "y = " << Grid[idx].y << "\n";
    for (int32_t i = 0; i < GlobalNz; i++)
    {
        for (int32_t j = 0; j < GlobalNx; j++)
        {
            // cout << fixed << std::setprecision(2) << "(" << grid.Grid[idx].x << ";" << grid.Grid[idx].z << ";" << grid.Grid[idx].y << ") ";
            printf("(%.2f;%.2f)", Grid[idx].x, Grid[idx].z);
            idx++;
        }
        cout << "\n";
    }
}