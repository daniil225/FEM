#include <iostream>
#include "PointInfo.h"


int main()
{
    Info info1;
    InfoManeger::ClearInfo(info1);

    InfoManeger::SetFictitious(info1, 1);
    InfoManeger::SetAreaInfo(info1, 3);
    InfoManeger::SetAreaInfo(info1, 1);
    InfoManeger::SetAreaInfo(info1, 2);
    InfoManeger::SetAreaInfo(info1, 1);

    InfoManeger::SetBoundInfo(info1, 1, 1);
    InfoManeger::SetBoundInfo(info1, 2, 3);
    InfoManeger::SetBoundInfo(info1, 1, 1);

    InfoManeger::PrintInfo(info1);

   
    return 0;
}