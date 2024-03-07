#include <iostream>
#include "Grid3D_StreightQuadPrismatic.h"

int main()
{

    Grid3D_StreightQuadPrismatic Grid;

    GridStatus Status = Grid.Load("Area1.txt");
    if(Status.GetState() == State::OK)
    {
        Status = Grid.GenerateGrid();
        if(Status.GetState() == State::OK)
        {
            //Grid3D_Size GridSize = Grid.GetGridSize();
            Grid.PrintGridSlice(0);
            //InfoManeger::PrintInfo(Grid[3].info);
            Status = Grid.DivideGrid(2);
            Status = Grid.ReGenerateGrid();
            if(Status.GetState() == State::OK)
                Grid.PrintGridSlice(0);
        }
        else
        {
            cout << Status.GetMsg();
        }
    }
    else
    {
        cout << Status.GetMsg();
    }
   
    return 0;
}