#ifndef POINT_INFO_H_
#define POINT_INFO_H_

#include <cstdint>
/* Определеим Структуры и функции для управления информацией о точке в расчетной области */
/* В общем случае данный класс применим для объектов в 2D и 3D областях */

/* Развертка информации о типе КУ  */
struct BoundInfo
{
    int8_t size = 0;                // Сколько типов задано
    int8_t Cond[4]{0, 0, 0, 0};     // Набор формул определяющих КУ максимум 2
    int8_t TypeCond[4]{0, 0, 0, 0}; // Тип Краевого условия
};

/* Развертка информации об области */
struct AreaInfo
{
    int8_t size = 0;
    int8_t Cond[8]{0, 0, 0, 0, 0, 0, 0, 0};
};

struct Comand
{
    typedef const uint32_t ComandType;
    /* Clear comand */
    static ComandType BaseInfoClear = 0x00;
    static ComandType AreaInfoClear = 0x00000000;
    static ComandType BoundInfoClear = 0x0000;
    static ComandType TypeBoundCondClear = 0x00;

    /* Base Info comand */
    static ComandType SetZeroFiFictitious = 0xFE;
    static ComandType GetFiFictitious = 0x01;

    /* Area Comand */
    static ComandType SetZeroAreaCount = 0xE1;
    static ComandType GetAreaCount = 0x1E;
    static ComandType ShiftAreaCountBits = 1;
    static ComandType RightBoundAreaCount = 8;
    static ComandType RightBoundValAreaInfo = 15;
    static ComandType BaseAreaShift = 4;
    static ComandType GetAreaInfoBaseComand = 0x0000000F;

    /* Bound Comand */
    static ComandType SetZeroBoundCount = 0x1F;
    static ComandType ShiftBoundCountBits = 5;
    static ComandType GetBoundCount = 0xE0;
    static ComandType RightBoundBoundCount = 4;
    static ComandType RightBoundValBoundInfo = 15;
    static ComandType RightBoundValBoundType = 3;
    static ComandType BaseBoundShift = 4;
    static ComandType BaseTypeBoundShift = 2;
    static ComandType GetBoundInfoBaseComand = 0x000F;
    static ComandType GetTypeBoundBaseComand = 0x03;
};

struct Info
{
    uint64_t BaseInfo : 8;      /* Базовая информация */
    uint64_t AreaInfo : 32;     /* Информация об области определяющей параметры ДУ */
    uint64_t BoundInfo : 16;    /*Краевые условия просто набор формул */
    uint64_t TypeBoundCond : 8; /* Тип Краевого условия */
};

class InfoManeger
{
private:

    /* Установка Количества различныйх типов областей которым принадлежит точка
        @param
        Info& info
        uint8_t val - диапазон [0;8]
        ret void
    */
    static void SetAreaCountBits(Info &info, uint32_t val);

    /* Установка Количества различных типов границ к которым примыкает данная точка
        @param
        Info& info
        uint8_t val - диапазон [0;4]
        ret void
    */
    static void SetBoundCountBits(Info &info, uint32_t val);

    /* Получить количество различных областей к которым примыкает точка
        @param
        const Info &info
        ret uint32_t - информациоя о том к скольки областям примыкает  диапазон [0;8]
        0 - ни к одной
        1 - к одной
        ...
    */
    static uint32_t GetAreaCount(const Info &info);

    /* Получить количество различных областей к которым примыкает граница
        @param
        const Info &info
        ret uint32_t - информациоя о том к скольки Границам примыкает(различным)  диапазон [0;4]
        0 - ни к какой
        1 - к одной
        2 -
        3 -
        4 -
    */
    static uint32_t GetBoundCount(const Info &info);

public:

    /* Очистка структуры */
    /*
        @param
        Info& info
        @return void
        Очистка структуры
    */
    static void ClearInfo(Info &info);

    /* Установка фиктивный/не фиктивный узел
        @param
        Info& info
        uint8_t val true - не фиктивный false - фиктивный
        ret void
    */
    static void SetFictitious(Info &info, uint8_t val);

    /*
        @param
        Info& info
        uint32_t val - набор формул задающий значения на границе диапазон [0,15]
        uint32_t boundType - тип Ку диапазон  [0,3]
    */
    static void SetBoundInfo(Info &info, uint32_t val, uint32_t boundType);

    /*
        @param
        Info& info
        uint32_t val - набор формул задающий значения в области определяющей параметры ДУ  диапазон [0,15]
    */
    static void SetAreaInfo(Info &info, uint32_t val);

    /*
        @param
        const Info& info
        ret bool - false - фиктивный true - не фиктивный
    */
    static bool IsFiFictitious(const Info &info);

    /*
        @param
        const Info& info
        ret bool - false - не граница true - граница
    */
    static bool IsBound(const Info &info);

    /*
        @param
        const Info& info
        ret BoundInfo -структура содержащая информацию о типах Ку и формул которые их задают
    */
    static BoundInfo GetBoundInfo(const Info &info);

    /*
        @param
        const Info& info
        ret AreaInfo -структура содержащая информацию о типах областей которым принадлежит точка и как следствие это определяет парметры ДУ
    */
    static AreaInfo GetAreaInfo(const Info &info);
    /********************************************/

    /* Debug functions */
    static void PrintInfo(const Info &info);

    static void PrintBoundInfo(const BoundInfo &Bound);

    static void PrintAreaInfo(const AreaInfo &Area);

    ~InfoManeger() = default;
};

#endif
