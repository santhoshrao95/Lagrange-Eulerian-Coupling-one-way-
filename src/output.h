#pragma once

#include <iostream>
#include<fstream>
#include<vector>

void printlevel(size_t count, const long double totaltime, const long double delT);

void print(const long double lx, const long double ly, const long double lz,
            const long double dx, const long double dy, const long double dz,
            size_t nx, size_t ny, size_t nz,
            size_t nxp2, size_t nyp2, size_t nzp2, 
            std::vector<std::vector<std::vector<long double>>> &u,
            std::vector<std::vector<std::vector<long double>>> &v,
            std::vector<std::vector<std::vector<long double>>> &w,
            std::vector<std::vector<std::vector<long double>>> &p,
            size_t count, const long double Time, const long double delT);

void printParticle(std::vector<long double> &xp,std::vector<long double> &yp,std::vector<long double> &zp,
                    size_t np, size_t count);