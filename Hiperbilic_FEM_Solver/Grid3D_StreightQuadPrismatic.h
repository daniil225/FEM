#ifndef GRID3D_STREIGHTQUADPRISMATIC_H_
#define GRID3D_STREIGHTQUADPRISMATIC_H_

#include "GridI.h"
#include "PointInfo.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;
/*
    Данная Библиотека является ООП оберткой над структрурами и представляет более удобный интерфейс для упрвления
    сетками в 3D области. Сетка состоит из прямых призм с четырех угольным основанием
 */

/********* Входные данные  ***********/

/* Чтение и сетки из файла. Разбиение сетки.  */
/* Представленная структура используется для хранения сетки в виде шестигранных элементов
   Прямые призмы с четырехугольным основанием с основанием в плоскости (x;z)
*/

/* Формат файла
* Nx Nz
* x1,1 z1,1 x2,1 z1,2 ... xn,1 z1,n
* x1,2 z2,1 x2,2 z2,2 ... x2,n z2,n
* .....
* Ny
* y1, y2, ..., yk
*
* L - целое число - количество подобластей и далее L наборов чисел по 7 штук
* 1 число - номер формул определяющих параметры ДУ в подобласти
* 2 число - первая вертикальная ломанная определеяет правую границу
* 3 число - вторая вертикальная ломанная определяет левую границу
  их номера определяются в соответствии с правилом для k узла их номера в массивах точек будут: k, (Nx + k), (2*Nx + k), ... , ( (Ny-1)*Nx + k )
* 4 число - опеределяет горизонтальную границу снизу
* 5 число - определяет горизонтальную границу сверху
* 6 число - граница по y начало
* 7 число - граница по y конец
* Дальше 3 сторики задающие разбиение сетки
* P - целое число - количество подоблостей граничных условий и далее P наборов чисел по 8 штук
* 1 число - номер формул определяющих параметры граничных условий
* 2 число - Тип краевого условия 0 - соответсвует не заданным краевым
* 3 число - первая вертикальная ломанная определеяет правую границу
* 4 число - вторая вертикальная ломанная определяет левую границу
* 5 число - опеределяет горизонтальную границу снизу
* 6 число - определяет горизонтальную границу сверху
* 7 число - граница по y начало
* 8 число - граница по y конец
* !!! Отличительная особенность одна из координатных линий должна быть фиксированной
* N1 coef1 N2 coef 2 ...
* ...
* ...
* Первый набор задает разбиение каждого отрезка по х
* Второй набора задает разбиение по z
* Третий набор задает разбиение по у
*/
/***********************************************************************************************/

/*
    Структура данных для представления области
    Струкутура  базовой сетки на ее основе генерируется основная сетка для расчета
    Описываемая область не должна содержать внутренних узлов при ее генерации. Исключительно граничные. Далее алгоритм будет сам разбивать
    область и генерировать внутренние узлы сетки

    # Описание полей:
    + int32_t CountOfDivision = 1 - Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д
    + const int32_t SizeOfCalculationAreaElemet = 7
    + const int32_t SizeOfBoundsAreaElement = 8
    + const int32_t SizeOfDivideParam = 3
    + struct PointXZS - Структура для базовой точки:
        * double x = 0.0;
        * double z = 0.0;

    + struct DivideParamS Параметры структуры:
        * int32_t num = 0 - количество интервалов на которое нужно разделить отрезок
        * double coef = 0 - Коэффициент растяжения или сжатия

    + int32_t Nx = 0 - количество узлов вдоль горизонатального направления
    + int32_t Ny = 0 - количество узлов вдоль вертикального направления
    + int32_t Nz = 0 - количество узлов вдоль оси z
    + int32_t L = 0 - Количество подоблостей
    + int32_t P = 0 - Количество видов границ
    + vector<vector<PointXZS>> BaseGridXZ - Базовая сетка в плосткости XZ
    + vector<double> BaseGridY - Базовая сетка по Y
    + vector<vector<int32_t>> CalculationArea - Массив расчетных областей
        Формат: | номер формул | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
    + vector<vector<int32_t>> BoundsArea - Границы
        Формат:  | номер формул | | тип КУ | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
    + vector<vector<DivideParamS>> DivideParam - Массив для разбиения
        Формат:
            DivideParam[0] - массив для разбиения по оси x размер Nx-1
            DivideParam[1] - массив для разбиения по оси z размер Nz-1
            DivideParam[2] - массив для разбиения по оси y размер Ny-1
    + bool isReadyToUse = false - Готовность структуры к использованию. Поведение переменной и проверка валидности структуры целиком лежит на программисте
 */
struct BaseGrid3DStreightQuadPrismatic
{
    int32_t CountOfDivision = 1; // Количество делений базовой настройки сетки; 1 - соответсвует тому что делений не было 2 тому что исходная область разделена в двое и.т.д
    const int32_t SizeOfCalculationAreaElemet = 7;
    const int32_t SizeOfBoundsAreaElement = 8;
    const int32_t SizeOfDivideParam = 3;

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
    int32_t L = 0;  // Количество подоблостей
    int32_t P = 0;  // Количество видов границ

    vector<vector<PointXZS>> BaseGridXZ; //  Базовая сетка в плосткости XZ
    vector<double> BaseGridY;            //  Базовая сетка по Y
    /*
    | номер формул | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
    */
    vector<vector<int32_t>> CalculationArea; // Массив расчетных областей

    /*
        | номер формул | | тип КУ | граница по x(start) | граница по x(end) | | граница по z(start) | граница по z(end) | граница по y(start) | граница по y(end) |
    */
    vector<vector<int32_t>> BoundsArea; // Границы

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
    BaseGrid3DStreightQuadPrismatic() = default;
    BaseGrid3DStreightQuadPrismatic(BaseGrid3DStreightQuadPrismatic &&) = default;            // Для возврата из функции с захватом параметра
    BaseGrid3DStreightQuadPrismatic &operator=(BaseGrid3DStreightQuadPrismatic &&) = default; // Аналогично присваивание с захватом ресурсов

    BaseGrid3DStreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)
    {
        this->BaseGridXZ = baseGrid_.BaseGridXZ;
        this->BaseGridY = baseGrid_.BaseGridY;
        this->BoundsArea = baseGrid_.BoundsArea;
        this->CalculationArea = baseGrid_.CalculationArea;
        this->CountOfDivision = baseGrid_.CountOfDivision;
        this->DivideParam = baseGrid_.DivideParam;
        this->isReadyToUse = baseGrid_.isReadyToUse;
        this->L = baseGrid_.L;
        this->P = baseGrid_.P;
        this->Nx = baseGrid_.Nx;
        this->Ny = baseGrid_.Ny;
        this->Nz = baseGrid_.Nz;
    }

    BaseGrid3DStreightQuadPrismatic &operator=(const BaseGrid3DStreightQuadPrismatic &baseGrid_)
    {
        this->BaseGridXZ = baseGrid_.BaseGridXZ;
        this->BaseGridY = baseGrid_.BaseGridY;
        this->BoundsArea = baseGrid_.BoundsArea;
        this->CalculationArea = baseGrid_.CalculationArea;
        this->CountOfDivision = baseGrid_.CountOfDivision;
        this->DivideParam = baseGrid_.DivideParam;
        this->isReadyToUse = baseGrid_.isReadyToUse;
        this->L = baseGrid_.L;
        this->P = baseGrid_.P;
        this->Nx = baseGrid_.Nx;
        this->Ny = baseGrid_.Ny;
        this->Nz = baseGrid_.Nz;
        return *this;
    }

    BaseGrid3DStreightQuadPrismatic(BaseGrid3DStreightQuadPrismatic &) = delete;
    BaseGrid3DStreightQuadPrismatic &operator=(BaseGrid3DStreightQuadPrismatic &) = delete;
    /************************************************************************/

    ~BaseGrid3DStreightQuadPrismatic() = default;
};

/*
    Геометрическая точка или узел конечного элемента
    Информация о границе конкретный набор формул + информация об узле фиктивный или нет + информация о том граничный ли узел или нет
    + информация о конечном элементе + тип КУ + точка в пространтве
    + Вся эта информация хранится в соответсвующих битах числа битовая картат ниже
    # Поля
    + Info info - Структура содержит информацию о точке
    + double x = 0.0;
    + double y = 0.0;
    + double z = 0.0;
*/
struct Point
{
    Info info;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/*
    Структура содержащая количество точек на каждой из координатных линий + общая размерность + количество конечных элементов
    # Поля структуры:
    + int32_t Dim = 0;      - Размерность сетки(количество узлов) и СЛАУ в то же время
    + int32_t FEMCount = 0; - Количество конечных элементов
    + int32_t GlobalNx = 0; - Сумарное количество узлов по оси Х
    + int32_t GlobalNy = 0; - Сумарное количество узлов по оси У
    + int32_t GlobalNz = 0; - Сумарное количество узлов по оси Z
*/
struct Grid3D_Size
{
    int32_t Dim = 0;      // Размерность сетки(количество узлов) и СЛАУ в то же время
    int32_t FEMCount = 0; // Количество конечных элементов
    int32_t GlobalNx = 0; // Сумарное количество узлов по оси Х
    int32_t GlobalNy = 0; // Сумарное количество узлов по оси У
    int32_t GlobalNz = 0; // Сумарное количество узлов по оси Z

    friend ostream &operator<<(ostream &os, Grid3D_Size &Grid3D_Size_param);
};

/*
    Конечный элемент в виде прямой призмы
    В структуре содержится вся информация о конечном элементе, а именно
    + static const int32_t FinitElementSize - размер конечного элемента (равен 8)
    + static const int32_t BoundCount - количество границ всего их 6
    + bool isFictive = false - по умолчанию фиктивный
    + int32_t AreaInfo = -1; -  Номер формулы задающий парметры ДУ на заданной области
    + BoundS Bound[BoundCount] - Вся информация о границе конечного элемента
    + Point e[FinitElementSize];  - Точки на конечном элементе нужны для задания Базисных функций
    + int32_t GlobalIdx[FinitElementSize] - Глобальная нумерация на конечном элементе. Индекс в массиве а в этой ячейке хранится глобальный номер в матрице/сетке

    Описание структуры struct BoundS:
    +  bool IsBound = false;
    +  int32_t LocalIdx[4]{-1, -1, -1, -1};  - Локальная нумерация для удобного сбора Локальной матрицы
    +  int32_t GlobalIdx[4]{-1, -1, -1, -1}; - Глобальная нумерация для удобного соотнесения с глобальной матрицей
    +  int32_t BoundType = -1;               - Тип КУ
    +  int32_t BoundInfo = -1;               - Формула задающая КУ
    +  Всего на конечном элементе может быть 6 границ (8-и угольник)
        Со следующей локальной нумерацией
        {1,2,3,4}, {5,6,7,8}, {1,3,5,7}, {2,4,6,8}, {1,2,5,6}, {3,4,7,8} - в нумерации с 1
        на каждой границе могут быть заданы свои условия. При этом краевое условие определяется 2 числоми
        1 число - номер формул задающих КУ
        2 число - Тип КУ
        Заведем 6 структур для каждой границы отдельно
*/
struct Finit_Element_StreightQuadPrismatic
{
    int32_t FinitElementNum = -1; // Номер конечного элемента в сетке
    struct BoundS
    {
        bool IsBound = false;
        int32_t LocalIdx[4]{-1, -1, -1, -1};  // Локальная нумерация для удобного сбора Локальной матрицы
        int32_t GlobalIdx[4]{-1, -1, -1, -1}; // Глобальная нумерация для удобного соотнесения с глобальной матрицей
        int32_t BoundType = -1;               // Тип КУ
        int32_t BoundInfo = -1;               // Формула задающая КУ
    };
    // Локальные узлы на элементе в них лежат глобальные координаты, но для построения алгоритма обработки лучше возвращать такой елемент
    // Из такой структуры можно легко получить всю информацию

    static const int32_t FinitElementSize = 8; // Размер элемента
    static const int32_t BoundCount = 6;
    /* Информация о конечном элементе */
    bool isFictive = false; // по умолчанию фиктивный
    int32_t AreaInfo = -1;  // Номер формулы задающий парметры ДУ на заданной области
    /* Граница */
    /* Всего на конечном элементе может быть 6 границ (8-и угольник)
        Со следующей локальной нумерацией
        {1,2,3,4}, {5,6,7,8}, {1,3,5,7}, {2,4,6,8}, {1,2,5,6}, {3,4,7,8} - в нумерации с 1
        на каждой границе могут быть заданы свои условия. При этом краевое условие определяется 2 числоми
        1 число - номер формул задающих КУ
        2 число - Тип КУ
        Заведем 6 структур для каждой границы отдельно

        Структура в следующем формате:
        bool IsBound = false; - Граница ли это. По умолчанию нет
        int LocalIdx[4]; - Локальная нумерация
        int GlobalIdx[4]; - Соответсвующая глобальная нумерация для внесения результатов в Глобальную матрицу
        int BoundType = -1; - Тип КУ
        int BoundInfo = -1; - Набор фомул
    */
    BoundS Bound[BoundCount]; // Вся информация о границе конечного элемента

    Point e[FinitElementSize]; // Точки на конечном элементе нужны для задания Базисных функций

    // Глобальная нумерация на конечном элементе
    // Индекс в массиве а в этой ячейке хранится глобальный номер в матрице/сетке
    int32_t GlobalIdx[FinitElementSize];

    /*
        @param void
        @return void\\
        @result Вывод информации о конечном элементе 
    */
    void PrintElement() const
    {
        std::cout << "Finit Element: \n";
        for (int i = 0; i < 8; i++)
            std::cout << "Global idx = " << GlobalIdx[i] << " Cord = (" << e[i].x << ";" << e[i].z << ";" << e[i].y << ")\n";
        std::cout << "\n";

        std::cout << "Info: \n";
        std::cout << "Is Fictitios: " << isFictive << "\n";
        std::cout << "Type Area: " << AreaInfo << "\n";

        std::cout << "\n";
        std::cout << "Info About Bound: \n";
        for (int i = 0; i < BoundCount; i++)
        {
            std::cout << "Is Bound: " << Bound[i].IsBound << "\n";
            std::cout << "Bound Formula: " << Bound[i].BoundInfo << "  ";
            std::cout << "Bound Type: " << Bound[i].BoundType << "\n";
            std::cout << "Local Idx: {";
            for (int j = 0; j < 4; j++)
            {
                std::cout << Bound[i].LocalIdx[j] << ";";
            }
            std::cout << "} \n";

            std::cout << "Global Idx: {";
            for (int j = 0; j < 4; j++)
            {
                std::cout << Bound[i].GlobalIdx[j] << ";";
            }
            std::cout << "} \n\n";
        }
    }
};

class Grid3D_StreightQuadPrismatic : public GridI<BaseGrid3DStreightQuadPrismatic, Finit_Element_StreightQuadPrismatic>
{
private:
    /*   private Variables     */
    GridStatus Status;

    /*         все величины получаются после генерации сетки       */
    int32_t Dim = 0;      // Размерность сетки(количество узлов) и СЛАУ в то же время
    int32_t FEMCount = 0; // Количество конечных элементов
    int32_t GlobalNx = 0; // Сумарное количество узлов по оси Х
    int32_t GlobalNy = 0; // Сумарное количество узлов по оси У
    int32_t GlobalNz = 0; // Сумарное количество узлов по оси Z
    /**************************************************************/

    /* Базовая сетка */
    BaseGrid3DStreightQuadPrismatic baseGrid;

    // Массив точек получающийся при генерации конечных элементов
    vector<Point> Grid;

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

    /*
       @param void
       @return void
       @note Учет фиктивных узлов
   */
    void DivisionIntoSubAreas() noexcept;

    /*
        @param void
        @return void
        @note Функция учетка типа КУ и установка факта является ли элемент граничным
    */
    void DivisionIntoSubBounds() noexcept;

    /*
        @param int i - номер
        @param int axis - соответсвующая ость 0 - x, 1 - z 2 - y
        @return int Величина скачка в сетке
        @note Вернет число соответсвующее стартовой позиции по сути это скачок
    */
    int32_t Getlevel(int32_t i, int32_t axis) const noexcept;

    /**************************************************************/

protected:
public:
    /* Конструкторы объекта  */
    Grid3D_StreightQuadPrismatic() = default;

    /* Передача по ссылке разрешена но объект нельзя изменять */
    Grid3D_StreightQuadPrismatic(const Grid3D_StreightQuadPrismatic &) = default;

    /* Передача по сслке изменения разрешены */
    Grid3D_StreightQuadPrismatic(Grid3D_StreightQuadPrismatic &) = delete;

    /*
        @param const string &filename - Текстовый файл с разметкой
        @return Constructed Object Grid3D_StreightQuadPrismatic
        @note Создается инициализированный объект с полностью построенной сеткой
    */
    explicit Grid3D_StreightQuadPrismatic(const string &filename);

    /*
        @param const BaseGrid3DStreightQuadPrismatic& baseGrid_ - базовая сетка
        @return Constructed Object Grid3D_StreightQuadPrismatic
        @note Создается инициализированный объект с полностью построенной сеткой
    */
    explicit Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_);
    /**************************************************************/

    /* Методы интерфеса GridI */

    /*
        @param const string &filename - Текстовый файл с разметкой
        @return GridStatus\\
        @result Загрузка базовой сетки из файла
        @warning Валидация входных данных не предусмотрена
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus Load(const string &filename) noexcept;

    /*
        @param void
        @return GridStatus\\
        @result Генерация сетки
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus GenerateGrid() noexcept;

    /*
        @param const int coef - Коэффициент дробления
        @return GridStatus\\
        @result Дробление сетки в заданное количество раз
        @warning Производит исключительно установку новых параметров дробление. Для их применения нужно вызвать метод ReGenerateGrid()
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]] virtual GridStatus DivideGrid(const int32_t coef) noexcept;

    /*
        @param void
        @return GridStatus\\
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
    inline virtual BaseGrid3DStreightQuadPrismatic GetBaseGrid() const noexcept { return baseGrid; }

    /*
        @param int idx - Индекс центральной точки в глобальной нумерации
        @return Finit_Element_StreightQuadPrismatic - Структура содержащая всю необходимую информацию о конечном элементе
        @result Получить конечный элемент с полным его описанием (Границы, Область)
    */
    virtual Finit_Element_StreightQuadPrismatic GetElement(const int32_t idx) const noexcept;

    /*
        @param const BaseGrid3DStreightQuadPrismatic& baseGrid_ - Базовая сетка области
        @return GridStatus\\
        @result Установка параметров базовой сетки
        @warning Не изменяет уже построенной сетки. Для применения изменений нужно вызвать метод ReGenerateGrid()
        @note Результат работы метода нельзя игнорировать
    */
    [[nodiscard]]  GridStatus SetBaseGrid(const BaseGrid3DStreightQuadPrismatic &baseGrid_) noexcept;

    /*
        @param int idx - индекс точки в МКЭ сетке
        @return Point& - точка в МКЭ Области
        @result -
    */
    inline Point &operator[](const int32_t idx) noexcept { return Grid[static_cast<uint64_t>(idx)]; }

    /**************************************************************/

    /*  методы специфичные для данного класса */

    /*
        @param    const double x - Точка в области
        @param    const double y - Точка в области
        @param    const double z - Точка в области
        @return Finit_Element_StreightQuadPrismatic - Структура содержащая всю необходимую информацию о конечном элементе\\
        @result Получить конечный элемент с полным его описанием (Границы, Область)\\
        @result Если нужно получить конечный элемент по заданной координате, нужно например при расчете функции после получения решения
        @note Результат работы метода нельзя игнорировать
    */
    Finit_Element_StreightQuadPrismatic GetElement(const double x, const double y, const double z) const noexcept;

    /*
        @param void
        @return Grid3D_Size
    */
    inline Grid3D_Size GetGridSize() const { return Grid3D_Size{Dim, FEMCount, GlobalNx, GlobalNy, GlobalNz}; };

    /*
        @param int32_t level - Номер секущей плоскости
        @return void\\
        @result Вывод плоскости XZ по точкам
     */
    void PrintGridSlice(int32_t level) const;
    /**************************************************************/

    /* Копирование/звхват объекта запрещены */

    Grid3D_StreightQuadPrismatic(Grid3D_StreightQuadPrismatic &&) = delete;
    Grid3D_StreightQuadPrismatic &operator=(const Grid3D_StreightQuadPrismatic &) = delete;
    Grid3D_StreightQuadPrismatic &operator=(Grid3D_StreightQuadPrismatic &) = delete;
    Grid3D_StreightQuadPrismatic &operator=(Grid3D_StreightQuadPrismatic &&) = delete;
    /**************************************************************/

    ~Grid3D_StreightQuadPrismatic() = default;
};

#endif