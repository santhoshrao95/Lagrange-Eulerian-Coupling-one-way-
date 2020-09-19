#pragma once

#include"NSSolverRK3_CN.h"

#include "Particle.h"
#include<vector>
/** @brief This is a class for running 3d-Cavity flow problem.
 */
class CavityFlow{
    public:
    /** Default constructor */
    CavityFlow();
    /**Initializes all the coordinates. 32x32x32 grid with Re = 400, Time = 10sec. 
     */
    CavityFlow(const long double xL,const long double yL,const long double zL,
               const long double timeLength, const long double initialTimestep, const long double Re, const long double dia,
               size_t nx, size_t ny, size_t nz, size_t np);
    /**parameterized constructor which takes inputs from user 
     */
    /** Run program calls the nssolver
     */
    void run();
    private:
    const long double Lx, Ly, Lz, timeLength,initialTimestep, Re, dia;
    size_t nx, ny, nz, np;
    size_t nxp2, nyp2, nzp2;
    const long double dx = Lx / (nxp2-1), dy = Ly / (nyp2-1), dz = Lz / (nzp2-1);
    std::vector<std::vector<std::vector<long double>>> u;
    std::vector<std::vector<std::vector<long double>>> v;
    std::vector<std::vector<std::vector<long double>>> w;
    std::vector<std::vector<std::vector<long double>>> p;
    std::vector<long double> xe;
    std::vector<long double> ye; 
    std::vector<long double>ze;
    std::vector<long double>xc;
    std::vector<long double>yc;
    std::vector<long double>zc;
    std::vector<long double>xp;
    std::vector<long double>yp;
    std::vector<long double>zp;

};