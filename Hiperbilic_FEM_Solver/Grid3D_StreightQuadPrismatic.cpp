#include "Grid3D_StreightQuadPrismatic.h"

#include <iostream>
#include <cmath>
#include <functional>

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

void Grid3D_StreightQuadPrismatic::GenerateBaseGrid(GridStatus &status) noexcept
{
    struct SettingForDivide
    {
        double step = 0; // Шаг на отрезке
        double coef = 0; // Коэффициент увеличения шага
        int num = 0;     // Количество интервалов идем то num-1 и потом явно вставляем элемент

        /* Копирование и присваивание */
        SettingForDivide() = default;
        SettingForDivide(const SettingForDivide &) = default;
        SettingForDivide &operator=(const SettingForDivide &) = default;
        SettingForDivide(SettingForDivide&&) = default;
        SettingForDivide &operator=(SettingForDivide &&) = default;
    };

    /* Расчитываем шаг для сетки  */
    /*
        @param
        int32_t i - Номер массива от 0 до 2
        int32_t j - Номер элемента в массиве
        double left - левая грани отрезка
        double right - правая граница отрезка
        @return: SettingForDivide -  структура с вычесленными параметрами деления сетки
    */
    std::function<SettingForDivide(int32_t, int32_t, double, double)>  CalcSettingForDivide = 
    [&](int32_t i, int32_t j, double left, double right) -> SettingForDivide
    {
        SettingForDivide res;
        int32_t num = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].num;
        double coef = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].coef;

        if (std::abs(coef-1.0) > eps)
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
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubAreas(GridStatus &status) noexcept
{
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubBounds(GridStatus &status) noexcept
{
}

int Grid3D_StreightQuadPrismatic::Getlevel(int32_t i, int32_t axis) const noexcept
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
        GenerateBaseGrid(Status);
    else
        return Status;
    // Учет фиктивных узлов
    if (Status.GetState() == State::OK)
        DivisionIntoSubAreas(Status);
    else
        return Status;

    // Учет КУ и расстановка границ
    if (Status.GetState() == State::OK)
        DivisionIntoSubBounds(Status);
    else
        return Status;

    return Status;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::DivideGrid(const int coef) noexcept
{
    return Status;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::ReGenerateGrid() noexcept
{
    return Status;
}

FEM_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const int32_t idx) const noexcept
{
    FEM_StreightQuadPrismatic FEMElement;
    return FEMElement;
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

FEM_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const double x, const double y, const double z) const noexcept
{
    FEM_StreightQuadPrismatic FEMElement;
    return FEMElement;
}