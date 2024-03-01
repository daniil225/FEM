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
            Grid.PrintGridSlice(1);
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