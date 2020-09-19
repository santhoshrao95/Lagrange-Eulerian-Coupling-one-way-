#pragma once

#include<vector>
#include<iostream>

void thomas_rhs(std::vector<std::vector<std::vector<long double>>> &VTemp,
                    std::vector<std::vector<std::vector<long double>>> &m_Vn,
                    std::vector<std::vector<std::vector<long double>>> &FTemp,
                    size_t nxp2, size_t nyp2, size_t nzp2);

void thomasSolver(std::vector<std::vector<std::vector<long double>>> &uTemp1,
                    std::vector<std::vector<std::vector<long double>>> &vTemp1,
                    std::vector<std::vector<std::vector<long double>>> &wTemp1,
                    std::vector<std::vector<std::vector<long double>>> &rTemp1,
                    std::vector<std::vector<std::vector<long double>>> &rTemp2,
                    std::vector<std::vector<std::vector<long double>>> &rTemp3,
                    size_t nx, size_t ny, size_t nz,
                    const long double dx, const long double dy, const long double dz,
                    const long double Re, const long double dtrk, size_t n);