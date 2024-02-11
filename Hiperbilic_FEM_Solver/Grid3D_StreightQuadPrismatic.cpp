#include <iostream>
#include "Grid3D_StreightQuadPrismatic.h"

/* Private section */
void Grid3D_StreightQuadPrismatic::GetTotalNumberOfNodes() noexcept
{
    for (int i = 0; i < baseGrid.Nx - 1; i++)
        GlobalNx += baseGrid.DivideParam[0][i].num;

    for (int i = 0; i < baseGrid.Nz - 1; i++)
        GlobalNz += baseGrid.DivideParam[1][i].num;

    for (int i = 0; i < baseGrid.Ny - 1; i++)
        GlobalNy += baseGrid.DivideParam[2][i].num;

    FEMCount = GlobalNx * GlobalNy * GlobalNz;
    GlobalNx++;
    GlobalNy++;
    GlobalNz++;
    Dim = GlobalNx * GlobalNy * GlobalNz;
}

void Grid3D_StreightQuadPrismatic::GenerateBaseGrid(GridStatus &status) noexcept
{
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubAreas(GridStatus &status) noexcept
{
}

void Grid3D_StreightQuadPrismatic::DivisionIntoSubBounds(GridStatus &status) noexcept
{
}

int Grid3D_StreightQuadPrismatic::Getlevel(int i, int axis) const noexcept
{
    int res = 0;
    for (int k = 0; k < i; k++)
        res += baseGrid.DivideParam[axis][k].num;
    return res;
}

/* Public section */

Grid3D_StreightQuadPrismatic::Grid3D_StreightQuadPrismatic(const string &filename)
{
    // Если загрузка не выполнена успешно то генерируем исключение с возвратом статуса 
    if(Load(filename).GetState() != State::OK)
        throw Status;

}

Grid3D_StreightQuadPrismatic::Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)
{
    // Если структура не валидная то выбрасываем ошибку
    if(!baseGrid_.isReadyToUse)
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
        baseGrid.BaseGridXZ = vector(baseGrid.Nz, vector<BaseGrid3DStreightQuadPrismatic::PointXZS>(baseGrid.Nx));
        for (int i = 0; i < baseGrid.Nz; i++)
        {
            for (int j = 0; j < baseGrid.Nx; j++)
            {
                fin >> baseGrid.BaseGridXZ[i][j].x >> baseGrid.BaseGridXZ[i][j].z;
            }
        }
        /***********************/

        /* Базовая сетка по Y */
        fin >> baseGrid.Ny;
        baseGrid.BaseGridY = vector<double>(baseGrid.Ny);
        for (int i = 0; i < baseGrid.Ny; i++)
            fin >> baseGrid.BaseGridY[i];
        /***********************/

        /* Расчетные подобласти */
        fin >> baseGrid.L;
        baseGrid.CalculationArea = vector(baseGrid.L, vector<int>(baseGrid.SizeOfCalculationAreaElemet));

        for (int i = 0; i < baseGrid.L; i++)
        {
            for (int j = 0; j < baseGrid.SizeOfCalculationAreaElemet; j++)
            {
                fin >> baseGrid.CalculationArea[i][j];
                baseGrid.CalculationArea[i][j]--; // приведение нумераии с нуля
            }
        }
        /***********************/

        /* Описание Границ */
        fin >> baseGrid.P;
        baseGrid.BoundsArea = vector(baseGrid.P, vector<int>(baseGrid.SizeOfBoundsAreaElement));
        for (int i = 0; i < baseGrid.P; i++)
        {
            for (int j = 0; j < baseGrid.SizeOfBoundsAreaElement; j++)
            {
                fin >> baseGrid.BoundsArea[i][j];
                baseGrid.BoundsArea[i][j]--; // приведение нумераии с нуля
            }
        }
        /***********************/

        /* Правила дробления базовой сетки */
        baseGrid.DivideParam.resize(baseGrid.SizeOfDivideParam);
        baseGrid.DivideParam[0].resize(baseGrid.Nx - 1);
        baseGrid.DivideParam[1].resize(baseGrid.Nz - 1);
        baseGrid.DivideParam[2].resize(baseGrid.Ny - 1);

        for (int i = 0; i < baseGrid.Nx - 1; i++)
            fin >> baseGrid.DivideParam[0][i].num >> baseGrid.DivideParam[0][i].coef;

        for (int i = 0; i < baseGrid.Nz - 1; i++)
            fin >> baseGrid.DivideParam[1][i].num >> baseGrid.DivideParam[1][i].coef;

        for (int i = 0; i < baseGrid.Ny - 1; i++)
            fin >> baseGrid.DivideParam[2][i].num >> baseGrid.DivideParam[2][i].coef;
        /***********************/
        baseGrid.isReadyToUse = true; // Структуру можно использовать она инициализированна
        fin.close();
    }
    catch (std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        return Status;
    }
    catch(std::exception &e)
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

BaseGrid3DStreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetBaseGrid() const noexcept
{
    return baseGrid;
}

FEM_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const int32_t idx) const noexcept
{
    FEM_StreightQuadPrismatic FEMElement;
    return FEMElement;
}

[[nodiscard]] GridStatus Grid3D_StreightQuadPrismatic::SetBaseGrid(const BaseGrid3DStreightQuadPrismatic &baseGrid_) noexcept
{
    return Status;
}

Point &Grid3D_StreightQuadPrismatic::operator[](const int32_t idx) noexcept
{
    return Grid[static_cast<uint64_t>(idx)];
}

FEM_StreightQuadPrismatic Grid3D_StreightQuadPrismatic::GetElement(const double x, const double y, const double z) const noexcept
{
    FEM_StreightQuadPrismatic FEMElement;
    return FEMElement;
}