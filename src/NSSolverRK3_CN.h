#pragma once
//#ifndef NSSolverRK3_CN_H
//#define NSSolverRK3_CN_H
#include "BC.h"
#include "gs.h"
#include "ThomasAlgo.h"
#include "output.h"
#include "Particle.h"
#include<fstream>
#include <vector>
//#endif
/** It is the main solver which will call other solvers like mg, Thomas etc.
  @param u - u velocity
  @param v - v velocity
  @param w - w velocity
  @param p - pressure
  @param nx - number of cells in x-direction
  @param ny - number of cells in x-direction
  @param nz - number of cells in x-direction
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param Time - Total physical time to carry out the simulation
  \f$dx = 1/(nxp2-1)\f$.
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
 */

void NSSolverRK3_CN(std::vector<std::vector<std::vector<long double>>> &un,
                    std::vector<std::vector<std::vector<long double>>> &vn,
                    std::vector<std::vector<std::vector<long double>>> &wn,
                    std::vector<std::vector<std::vector<long double>>> &pn,
                    size_t nx, size_t ny, size_t nz, size_t np,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Time,const long double initialTimestep, long double Re,
                    std::vector<long double> &xe, std::vector<long double> &ye, std::vector<long double> &ze,
                    std::vector<long double> &xc, std::vector<long double> &yc, std::vector<long double> &zc,
                    std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp);

/** Calculates global Delta T to proceed in time.
  @param u - u velocity
  @param v - v velocity
  @param w - w velocity
  @param p - pressure
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
    \f$dx = 1/(nxp2-1)\f$.
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
 */
long double calculateDelt(std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    const long double dx, const long double dy, const long double dz);


/** Performs rksubstep
  @param u - u velocity
  @param v - v velocity
  @param w - w velocity
  @param p - pressure
  @param nx - number of cells in x-direction
  @param ny - number of cells in x-direction
  @param nz - number of cells in x-direction
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
  @param Re - Reynold's number
 
 */
void rksubstep(long double a1, long double a2, long double a3, long double dt_rk, long double delT,
                    std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    std::vector<std::vector<std::vector<long double>>> &p,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Re);

/** Calculates explicit part of the x-mommentum equation
  @param un - u velocity in current step
  @param vn - v velocity in current step
  @param wn - w velocity in current step
  @param pn - pressure in current step
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
  @param Re - Reynold's number
 
 */
void uRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &uTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);
/** Calculates explicit part of the y-mommentum equation
  @param un - u velocity in current step
  @param vn - v velocity in current step
  @param wn - w velocity in current step
  @param pn - pressure in current step
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
  @param Re - Reynold's number
 
 */

void vRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &vTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);
/** Calculates explicit part of the z-mommentum equation
  @param un - u velocity in current step
  @param vn - v velocity in current step
  @param wn - w velocity in current step
  @param pn - pressure in current step
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param dx - considering uniform grid, distance between the two consecutive pointns in x-direction
  @param dy - considering uniform grid, distance between the two consecutive pointns in y-direction
  @param dz - considering uniform grid, distance between the two consecutive pointns in z-direction 
  @param Re - Reynold's number
 
 */

void wRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &wTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re);

/** Calculates explicit part of the x-mommentum equation
  @param v1 - u velocity in current step
  @param v2 - v velocity in current step
  @param v3 - w velocity in current step
  @param nxp2 - number of cells in x-direction including two ghost cells on either ends
  @param nyp2 - number of cells in y-direction including two ghost cells on either ends
  @param nzp2 - number of cells in z-direction including two ghost cells on either ends
  @param f1 - a1 
  @param f2 - a2
  @param f3 - delT_rk \f$ delT\_rk = c\_3m delT \f$ 
 */
void var_advance(std::vector<std::vector<std::vector<long double>>> &v1,
                    std::vector<std::vector<std::vector<long double>>> &v2,
                    std::vector<std::vector<std::vector<long double>>> &v3,
                    const long double f1, const long double f2, const long double f3,
                     const long double nxp2, const long double nyp2, const long double nzp2);
void psource(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt);
void psource2(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt);
void pgrad(std::vector<std::vector<std::vector<long double>>> &pn,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        const long double dx, const long double dy, const long double dz,
                        size_t nxp2, size_t nyp2, size_t nzp2, const long double dtrk);

void print(const long double lx, const long double ly, const long double lz,
            const long double dx, const long double dy, const long double dz,
            size_t nx, size_t ny, size_t nz,
            size_t nxp2, size_t nyp2, size_t nzp2, 
            std::vector<std::vector<std::vector<long double>>> &u,
            std::vector<std::vector<std::vector<long double>>> &v,
            std::vector<std::vector<std::vector<long double>>> &w,
            std::vector<std::vector<std::vector<long double>>> &p,
            size_t count, const long Time, const long delT);

void printlevel(size_t count, const long double totaltime, const long double delT);

void interpolate(std::vector<long double> &velP,
                  std::vector<std::vector<std::vector<long double>>> &veln,
                  std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp,
                  std::vector<long double> &xg, std::vector<long double> &yg, std::vector<long double> &zg,
                  size_t np);

void updateK(long double a1,
                std::vector<long double> &K,
                std::vector<long double> &velP,
                const long double delT, size_t np);

void advance_pa(long double a2,
                std::vector<long double> &K,
                std::vector<long double> &xp,
                size_t np, long double low, long double high, int direction);


void rkparticlesubstep(long double a1, long double a2, long double a3,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<long double> &ku, std::vector<long double> &kv, std::vector<long double> &kw, 
                        const long double delT, size_t np,
                        std::vector<long double> &xe, std::vector<long double> &ye, std::vector<long double> &ze,
                        std::vector<long double> &xc, std::vector<long double> &yc, std::vector<long double> &zc,
                        std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp,
                        size_t nxp2, size_t nyp2, size_t nzp2);