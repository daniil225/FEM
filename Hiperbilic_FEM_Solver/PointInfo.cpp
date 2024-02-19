#include "PointInfo.h"
#include <iostream>

/* Privat Section */

void InfoManeger::SetAreaCountBits(Info &info, uint32_t val)
{
    if (val <= Comand::RightBoundAreaCount)
    {
        if (GetAreaCount(info) != 0)
            info.BaseInfo &= Comand::SetZeroAreaCount;
        info.BaseInfo |= static_cast<uint8_t>(val << Comand::ShiftAreaCountBits);
    }
    else
        throw "Error Bit Operation\n";
}

void InfoManeger::SetBoundCountBits(Info &info, uint32_t val)
{
    if (val <= Comand::RightBoundBoundCount)
    {
        if (GetBoundCount(info) != 0)
            info.BaseInfo &= Comand::SetZeroBoundCount;
        info.BaseInfo |= static_cast<uint8_t>(val << Comand::ShiftBoundCountBits);
    }
    else
        throw "Error Bit Operation\n";
}

uint32_t InfoManeger::GetAreaCount(const Info &info)
{
    return (info.BaseInfo & Comand::GetAreaCount) >> Comand::ShiftAreaCountBits;
}

uint32_t InfoManeger::GetBoundCount(const Info &info)
{
    return (info.BaseInfo & Comand::GetBoundCount) >> Comand::ShiftBoundCountBits;
}

/* Public section */

void InfoManeger::ClearInfo(Info &info)
{
    info.BaseInfo &= Comand::BaseInfoClear;
    info.AreaInfo &= Comand::AreaInfoClear;
    info.BoundInfo &= Comand::BoundInfoClear;
    info.TypeBoundCond &= Comand::TypeBoundCondClear;
}


void InfoManeger::SetFictitious(Info &info, uint8_t val)
{
    if (val <= 1)
    {
        if (IsFiFictitious(info))
            info.BaseInfo &= Comand::SetZeroFiFictitious;
        info.BaseInfo |= val;
    }
    else
        throw "Error Bit Operation\n";
}

void InfoManeger::SetBoundInfo(Info &info, uint32_t val, uint32_t boundType)
{
    if (val <= Comand::RightBoundValBoundInfo && boundType <= Comand::RightBoundValBoundType)
    {
        uint32_t Bnum = GetBoundCount(info);
        /* если не было значений */
        if (Bnum == 0)
        {
            info.BoundInfo |= static_cast<uint16_t>(val);
            info.TypeBoundCond |= static_cast<uint8_t>(boundType);
            SetBoundCountBits(info, 1);
        }
        else
        {
            BoundInfo Bound = GetBoundInfo(info);
            for (uint32_t i = 0; i < Bnum; i++)
            {
                // Проверяем дубликаты
                if ((Bound.Cond[i] ^ static_cast<int8_t>(val)) == 0) // Есть дубликат
                    return;
            }

            /* Новый элемент */
            info.BoundInfo |= static_cast<uint16_t>(val << Comand::BaseBoundShift * Bnum);
            info.TypeBoundCond |= static_cast<uint8_t>(boundType << Comand::BaseTypeBoundShift * Bnum);
            SetBoundCountBits(info, Bnum + 1);
        }
    }
    else
        throw "Error Bit Operation\n";
}

void InfoManeger::SetAreaInfo(Info &info, uint32_t val)
{
    if (val <= Comand::RightBoundValAreaInfo)
    {

        uint32_t Anum = GetAreaCount(info);
        if (Anum == 0)
        {
            info.AreaInfo |= val;
            SetAreaCountBits(info, 1);
        }
        else
        {
            AreaInfo Area = GetAreaInfo(info);
            for (uint32_t i = 0; i < Anum; i++)
            {
                // Проверяем дубликаты
                if ((Area.Cond[i] ^ static_cast<int8_t>(val)) == 0) // Есть дубликат
                    return;
            }

            /* Новый элемент */
            info.AreaInfo |= val << Comand::BaseAreaShift * Anum;
            SetAreaCountBits(info, Anum + 1);
        }
    }
    else
        throw "Error Bit Operation\n";
}

bool InfoManeger::IsFiFictitious(const Info &info)
{
    return info.BaseInfo & Comand::GetFiFictitious;
}

bool InfoManeger::IsBound(const Info &info)
{
    return GetBoundCount(info) == 0 ? false : true;
}

BoundInfo InfoManeger::GetBoundInfo(const Info &info)
{
    BoundInfo Bound;
    int32_t Bnum = static_cast<int32_t>(GetBoundCount(info));
    
    Bound.size = static_cast<int8_t>(Bnum);

    for (int32_t i = 0; i <= Bnum - 1; i++)
    {
        Bound.Cond[i] = static_cast<int8_t>((info.BoundInfo & (Comand::GetBoundInfoBaseComand 
                                            << static_cast<uint32_t>(i) * Comand::BaseBoundShift))
                                            >> static_cast<uint32_t>(i) * Comand::BaseBoundShift);

        Bound.TypeCond[i] = static_cast<int8_t>((info.TypeBoundCond & (Comand::GetTypeBoundBaseComand 
                                                << static_cast<uint32_t>(i) * Comand::BaseTypeBoundShift)) 
                                                >> static_cast<uint32_t>(i) * Comand::BaseTypeBoundShift);
    }

    return Bound;
}

AreaInfo InfoManeger::GetAreaInfo(const Info &info)
{
    AreaInfo Area;
    int32_t Anum = static_cast<int32_t>(GetAreaCount(info));
    Area.size = static_cast<int8_t>(Anum);

    for (int32_t i = 0; i <= Anum - 1; i++)
        Area.Cond[i] = static_cast<int8_t>((info.AreaInfo & (Comand::GetAreaInfoBaseComand 
                                            << Comand::BaseAreaShift * static_cast<uint32_t>(i)))
                                            >> Comand::BaseAreaShift * static_cast<uint32_t>(i));

    return Area;
}

void InfoManeger::PrintInfo(const Info &info)
{
    std::cout << "Info Struct\n";
    std::cout << "IsFiFictitious: " << IsFiFictitious(info) << "\n";
    std::cout << "AreaCount: " << GetAreaCount(info) << "\n";
    std::cout << "BoundCount: " << GetBoundCount(info) << "\n";
    PrintAreaInfo(GetAreaInfo(info));
    PrintBoundInfo(GetBoundInfo(info));
}

void InfoManeger::PrintBoundInfo(const BoundInfo &Bound)
{
    std::cout << "\nBound Info\n";
    for (int32_t i = 0; i < Bound.size; i++)
        std::cout << "K = " << i + 1 << " Num formula: " << static_cast<uint32_t>(Bound.Cond[i]) 
                  << " Type Bound: " << static_cast<uint32_t>(Bound.TypeCond[i]) << "\n";
}

void InfoManeger::PrintAreaInfo(const AreaInfo &Area)
{
    std::cout << "\nArea Info\n";
    for (int32_t i = 0; i < Area.size; i++)
        std::cout << "K = " << i + 1 << " Num formula: " << static_cast<uint32_t>(Area.Cond[i]) << "\n";
}