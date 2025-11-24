//
// Created by zw on 24-8-8.
//

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "globals.hh"

namespace B1
{

    constexpr G4int nx = 2;
    constexpr G4int ny = 2;
    constexpr G4int nz = 2;
    constexpr G4int nxy = nx*ny;
    constexpr G4int nxz = nx*nz;
    constexpr G4int Cells = nx * ny * nz;
}

#endif //CONSTANTS_H