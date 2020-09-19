//#pragma once

#include <vector>
#include <algorithm>
#include "NSSolverRK3_CN.h"

long double calculateDelt(std::vector<std::vector<std::vector<long double>>> &u,
                    std::vector<std::vector<std::vector<long double>>> &v,
                    std::vector<std::vector<std::vector<long double>>> &w,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    const long double dx, const long double dy, const long double dz)
{
    const long double cmax = 1.2; 
    long double umax = INT64_MIN, vmax = INT64_MIN, wmax = INT64_MIN;
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-2; i++)
      {
              if (umax<u[i][j][k])
                umax = u[i][j][k];
              /*if (vmax<v[i][j][k])
                vmax = v[i][j][k];
              if (wmax<w[i][j][k])
                wmax = w[i][j][k];*/
            }            
        }       
    }
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-2; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
              /*if (umax<u[i][j][k])
                umax = u[i][j][k];*/
              if (vmax<v[i][j][k])
                vmax = v[i][j][k];
              /*if (wmax<w[i][j][k])
                wmax = w[i][j][k];*/
            }            
        }       
    }
  for (size_t k = 1; k < nzp2-2; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
              /*if (umax<u[i][j][k])
                umax = u[i][j][k];
              if (vmax<v[i][j][k])
                vmax = v[i][j][k];*/
              if (wmax<w[i][j][k])
                wmax = w[i][j][k];
            }            
        }       
    }
    long double deltx = cmax * dx / umax;
    long double delty = cmax * dy / vmax;
    long double deltz = cmax * dz / wmax;

    return std::min(std::min(deltx, delty), deltz);
}

void uRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &uTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  //uTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-2; i++)
      {
        uaver2 = 0.5 * ( un[i+1][j][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i-1][j][k] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*uaver1 - uaver2*uaver2) * rdx;


        uaver2 = 0.5 * ( un[i][j][k+1] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i][j][k-1] );
        waver2 = 0.5 * ( wn[i+1][j][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i+1][j][k-1] + wn[i][j][k-1] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*waver1 - uaver2*waver2) * rdz;

        uaver2 = 0.5 * ( un[i][j+1][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i][j][k] + un[i][j-1][k] );
        waver2 = 0.5 * ( vn[i+1][j][k] + vn[i][j][k] );
        waver1 = 0.5 * ( vn[i+1][j-1][k] + vn[i][j-1][k] );
        uTemp[i][j][k] = uTemp[i][j][k] + ( uaver1*vaver1 - uaver2*vaver2) * rdy;

        uTemp[i][j][k] = uTemp[i][j][k] + rRe*( ( un[i+1][j][k] - un[i][j][k] )*rdx
                                              - ( un[i][j][k] - un[i-1][j][k] )*rdx )*rdx               
                                        + rRe*( ( un[i][j+1][k] - un[i][j][k] )*rdy  
                                              - ( un[i][j][k] - un[i][j-1][k] )*rdy )*rdy ;              
                                        /*+ rRe*( ( un[i][j][k+1] - un[i][j][k] )*rdy   
                                              - ( un[i][j][k] - un[i][j][k-1] )*rdy )*rdy ; */
      }
      
    }
    
  }

}

void vRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &vTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  // vTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-1; k++)
  {
    for (size_t j = 1; j < nyp2-2; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
        vaver2 = 0.5 * ( vn[i][j+1][k] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i][j-1][k] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( vaver1*vaver1 - vaver2*vaver2) * rdy;

        vaver2 = 0.5 * ( vn[i+1][j][k] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i-1][j][k] );
        uaver2 = 0.5 * ( un[i][j+1][k] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i-1][j+1][k] + un[i-1][j][k] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( uaver1*vaver1 - uaver2*vaver2) * rdx;

        vaver2 = 0.5 * ( vn[i][j][k+1] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j][k] + vn[i][j][k-1] );
        waver2 = 0.5 * ( wn[i][j+1][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j+1][k-1] + wn[i][j][k-1] );
        vTemp[i][j][k] = vTemp[i][j][k] + ( vaver1*waver1 - vaver2*waver2) * rdz;

        vTemp[i][j][k] = vTemp[i][j][k] + rRe*( ( vn[i+1][j][k] - vn[i][j][k] )*rdx
                                              - ( vn[i][j][k] - vn[i-1][j][k] )*rdx )*rdx 
                                        + rRe*( ( vn[i][j+1][k] - vn[i][j][k] )*rdy 
                                              - ( vn[i][j][k] - vn[i][j-1][k] )*rdy )*rdy ;
                                       /* + rRe*( ( vn[i][j][k+1] - vn[i][j][k] )*rdy 
                                              - ( vn[i][j][k] - vn[i][j][k-1] )*rdy )*rdy ; */
      }
      
    }
    
  }

}
void wRhs(std::vector<std::vector<std::vector<long double>>> &un,
            std::vector<std::vector<std::vector<long double>>> &vn,
            std::vector<std::vector<std::vector<long double>>> &wn,
            std::vector<std::vector<std::vector<long double>>> &wTemp,
            size_t nxp2, size_t nyp2, size_t nzp2,
            const long double dx, const long double dy, const long double dz,
            const long double Re)
{
  long double uaver1 = 0.0, uaver2 = 0.0, vaver1 = 0.0, vaver2 = 0.0, waver1 = 0.0, waver2 = 0.0;
  long double rdx = 1.0/dx , rdy = 1.0/dy, rdz = 1.0/dz, rRe = 1.0/Re;
  //wTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  for (size_t k = 1; k < nzp2-2; k++)
  {
    for (size_t j = 1; j < nyp2-1; j++)
    {
      for (size_t i = 1; i < nxp2-1; i++)
      {
        waver2 = 0.5 * ( wn[i][j][k+1] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i][j][k-1] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( waver1*waver1 - waver2*waver2) * rdz;

        waver2 = 0.5 * ( wn[i+1][j][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i-1][j][k] );
        uaver2 = 0.5 * ( un[i][j][k+1] + un[i][j][k] );
        uaver1 = 0.5 * ( un[i-1][j][k+1] + un[i-1][j][k] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( uaver1*waver1 - uaver2*waver2) * rdx;
        
        waver2 = 0.5 * ( wn[i][j+1][k] + wn[i][j][k] );
        waver1 = 0.5 * ( wn[i][j][k] + wn[i][j-1][k] );
        vaver2 = 0.5 * ( vn[i][j][k+1] + vn[i][j][k] );
        vaver1 = 0.5 * ( vn[i][j-1][k+1] + vn[i][j-1][k] );
        wTemp[i][j][k] = wTemp[i][j][k] + ( vaver1*waver1 - vaver2*waver2) * rdy;

        wTemp[i][j][k] = wTemp[i][j][k] + rRe*( ( wn[i+1][j][k] - wn[i][j][k] )*rdx
                                              - ( wn[i][j][k] - wn[i-1][j][k] )*rdx )*rdx               
                                        + rRe*( ( wn[i][j+1][k] - wn[i][j][k] )*rdy
                                              - ( wn[i][j][k] - wn[i][j-1][k] )*rdy )*rdy ;
                                      /*  + rRe*( ( wn[i][j][k+1] - wn[i][j][k] )*rdy 
                                              - ( wn[i][j][k] - wn[i][j][k-1] )*rdy )*rdy ; */
      }
      
    }
    
  }

}


void var_advance(std::vector<std::vector<std::vector<long double>>> &v1, 
                    std::vector<std::vector<std::vector<long double>>> &v2,
                    std::vector<std::vector<std::vector<long double>>> &v3,
                    const long double f1, const long double f2, const long double f3,
                    const long double nxp2, const long double nyp2, const long double nzp2)
                    //v1-->qU-->(EU*delT) v2-->UTemp1(c2*qU) v3-->FTemp(EU)
{
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2; i++)
      {
        v1[i][j][k] = f1*v1[i][j][k] + v3[i][j][k]*f3;
        v2[i][j][k] = f2*v1[i][j][k];
      } 
    } 
  }  
}

void psource(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt)
{
  const long double rdx = 1.0/dx;
  const long double rdy = 1.0/dy;
  const long double rdz = 1.0/dz;
  const long double rrkdt = 1.0/rkdt;

  for (size_t i = 0; i < nxp2; i++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t k = 0; k < nzp2; k++)
      {
        v1Temp[i][j][k] = v1Temp[i][j][k] + un[i][j][k];
        v2Temp[i][j][k] = v2Temp[i][j][k] + vn[i][j][k];
        v3Temp[i][j][k] = v3Temp[i][j][k] + wn[i][j][k];
      }
      
    }
    
  }
  
  for (size_t k = 1; k < nzp2; k++)
  {
    for (size_t j = 1; j < nyp2; j++)
    {
      for (size_t i = 1; i < nzp2; i++)
      {
        src[i][j][k] = src[i][j][k] +( (v1Temp[i][j][k] - v1Temp[i-1][j][k]) * rdx \
                                                  + (v2Temp[i][j][k] - v2Temp[i][j-1][k]) * rdy \
                                                  + (v3Temp[i][j][k] - v3Temp[i][j][k-1]) * rdz ) * rrkdt;
      }
      
    }
    
  }
}

void psource2(std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<std::vector<std::vector<long double>>> &v1Temp,
                        std::vector<std::vector<std::vector<long double>>> &v2Temp,
                        std::vector<std::vector<std::vector<long double>>> &v3Temp,
                        std::vector<std::vector<std::vector<long double>>> &src,
                        size_t nx, size_t ny, size_t nz,
                        size_t nxp2, size_t nyp2, size_t nzp2,
                        const long double dx, const long double dy, const long double dz, const long double rkdt)
{
  const long double rdx = 1.0/dx;
  const long double rdy = 1.0/dy;
  const long double rdz = 1.0/dz;
  const long double rrkdt = 1.0/rkdt;

  
  for (size_t k = 1; k < nzp2; k++)
  {
    for (size_t j = 1; j < nyp2; j++)
    {
      for (size_t i = 1; i < nzp2; i++)
      {
        src[i][j][k] = src[i][j][k] +( (v1Temp[i][j][k] - v1Temp[i-1][j][k]) * rdx \
                                                  + (v2Temp[i][j][k] - v2Temp[i][j-1][k]) * rdy \
                                                  + (v3Temp[i][j][k] - v3Temp[i][j][k-1]) * rdz ) * rrkdt;
      }
      
    }
    
  }
}

void pgrad(std::vector<std::vector<std::vector<long double>>> &pn,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        const long double dx, const long double dy, const long double dz,
                        size_t nxp2, size_t nyp2, size_t nzp2, const long double dtrk)
{
  const long double rdx = 1.0/dx, rdy = 1.0/dy, rdz = 1.0/dz;
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2-1; i++)
      {
        un[i][j][k] = un[i][j][k] - ( pn[i+1][j][k]-pn[i][j][k] )*rdx * dtrk;
      } 
    } 
  }  
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2-1; j++)
    {
      for (size_t i = 0; i < nxp2; i++)
      {
        vn[i][j][k] = vn[i][j][k] - ( pn[i][j+1][k]-pn[i][j][k] )*rdy * dtrk;
      } 
    } 
  }
  for (size_t k = 0; k < nzp2-1; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2-1; i++)
      {
        wn[i][j][k] = wn[i][j][k] - ( pn[i][j][k+1]-pn[i][j][k] )*rdz * dtrk;
      } 
    } 
  }

}


void rksubstep(long double a1, long double a2, long double a3, long double dt_rk, long double delT,
                    std::vector<std::vector<std::vector<long double>>> &m_un,
                    std::vector<std::vector<std::vector<long double>>> &m_vn,
                    std::vector<std::vector<std::vector<long double>>> &m_wn,
                    std::vector<std::vector<std::vector<long double>>> &m_pn,
                    std::vector<std::vector<std::vector<long double>>> &m_qu,
                    std::vector<std::vector<std::vector<long double>>> &m_qv,
                    std::vector<std::vector<std::vector<long double>>> &m_qw,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Re)
{
  std::vector<std::vector<std::vector<long double>>> FTemp;
  std::vector<std::vector<std::vector<long double>>> rhs;
  std::vector<std::vector<std::vector<long double>>> uTemp1;
  std::vector<std::vector<std::vector<long double>>> vTemp1;
  std::vector<std::vector<std::vector<long double>>> wTemp1;
  std::vector<std::vector<std::vector<long double>>> rTemp1;
  std::vector<std::vector<std::vector<long double>>> rTemp2;
  std::vector<std::vector<std::vector<long double>>> rTemp3;
  
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  rhs = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  uTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  vTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  wTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  rTemp1 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  rTemp2 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  rTemp3 = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  uRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qu, uTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  
  vRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qv, vTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  wRhs(m_un,m_vn,m_wn,FTemp,nxp2,nyp2,nzp2,dx,dy,dz,Re);
  var_advance(m_qw, wTemp1, FTemp, a1, a2, delT, nxp2, nyp2, nzp2);
  FTemp = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  applyFlowVelocityBC(uTemp1, vTemp1, wTemp1,  m_pn, nx, ny, nz, nxp2, nyp2, nzp2);

  thomas_rhs(uTemp1, m_un, rTemp1, nxp2,nyp2,nzp2);
  thomas_rhs(vTemp1, m_vn, rTemp2, nxp2,nyp2,nzp2);
  thomas_rhs(wTemp1, m_wn, rTemp3, nxp2,nyp2,nzp2);

  applyFlowVelocityBC(rTemp1, rTemp2, rTemp3, m_pn, nx,ny,nz,nxp2,nyp2,nzp2);

  thomasSolver(uTemp1,vTemp1,wTemp1,rTemp1, rTemp2, rTemp3, nx, ny, nz, dx, dy, dz, Re, dt_rk, 1);
  thomasSolver(uTemp1,vTemp1,wTemp1,rTemp1, rTemp2, rTemp3, nx, ny, nz, dx, dy, dz, Re, dt_rk, 2);
  thomasSolver(uTemp1,vTemp1,wTemp1,rTemp1, rTemp2, rTemp3, nx, ny, nz, dx, dy, dz, Re, dt_rk, 3);

  applyFlowVelocityBC(uTemp1, vTemp1, wTemp1,  m_pn, nx, ny, nz, nxp2, nyp2, nzp2);

  //psource(m_un, m_vn, m_wn, uTemp1, vTemp1, wTemp1, rhs, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, dt_rk);
  psource2(m_un, m_vn, m_wn, uTemp1, vTemp1, wTemp1, rhs, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, dt_rk);

  gs(m_pn,rhs,nx,ny,nz,nxp2,nyp2,nzp2,dx,dy,dz);

  pgrad(m_pn,uTemp1, vTemp1, wTemp1,dx,dy,dz,nxp2,nyp2,nzp2,dt_rk);

  applyFlowVelocityBC(uTemp1, vTemp1, wTemp1,m_pn,nx,ny,nz,nxp2,nyp2,nzp2);

  m_un = uTemp1;
  m_vn = vTemp1;
  m_wn = wTemp1;

  std::vector<std::vector<std::vector<long double>>>().swap(FTemp);
  std::vector<std::vector<std::vector<long double>>>().swap(rhs);
  std::vector<std::vector<std::vector<long double>>>().swap(uTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(vTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(wTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(rTemp1);
  std::vector<std::vector<std::vector<long double>>>().swap(rTemp2);
  std::vector<std::vector<std::vector<long double>>>().swap(rTemp3);

}

int binarySearch(std::vector<long double> &x, long double target)
{
  int lo = 0, n = x.size(), hi = n-1, mid = 0;
  while (lo<hi)
  {
    mid = lo + (hi-lo)/2;

    if(x[mid]<target)
    {
      lo = mid+1;
    }
    else
    {
      hi = mid;
    }
    
  }
  return hi;
  

}


void interpolate(std::vector<long double> &velP,
                  std::vector<std::vector<std::vector<long double>>> &veln,
                  std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp,
                  std::vector<long double> &xg, std::vector<long double> &yg, std::vector<long double> &zg,
                  size_t np)
{
  int  i = 0, j=0, k=0;
  long double a1=0.0,a2=0.0,b1=0.0,b2=0.0,c1=0.0,c2=0.0;
  for (size_t m = 0; m < np; m++)
  {
    
    i=binarySearch(xg,xp[m]);
    j=binarySearch(yg,yp[m]);
    k=binarySearch(zg,zp[m]);
    a1 = abs((xp[m]-xg[i-1]) / (xg[i]-xg[i-1])), a2 = abs(1.0-a1);
    b1 = abs((yp[m]-yg[j-1]) / (yg[j]-yg[j-1])), b2 = abs(1.0-b1);
    c1 = abs((zp[m]-zg[k-1]) / (zg[k]-zg[k-1])), c2 = abs(1.0-c1);

    velP[m] = (( (veln[i][j-1][k-1]*a1 + veln[i-1][j-1][k-1]*a2) * c2 + (veln[i][j-1][k]*a1 + veln[i-1][j-1][k]) * c1 ) * b2) + 
              (( (veln[i][j][k-1]*a1 + veln[i-1][j][k-1]*a2) * c2 + (veln[i][j][k]*a1 + veln[i-1][j][k]) * c1 ) * b1);
     

  }
  
  
}

void updateK(long double a1,
                std::vector<long double> &K,
                std::vector<long double> &velP,
                const long double delT, size_t np)
{
  for (size_t m = 0; m < np; m++)
  {
    K[m] = a1*K[m] + delT*velP[m];
  }
  
}

void advance_pa(long double a2,
                std::vector<long double> &K,
                std::vector<long double> &xp,
                size_t np, long double low, long double high, int direction)
{

  for (size_t m = 0; m < np; m++)
  {
    xp[m] = xp[m] + (a2 * K[m]) ; 
    std::cout <<"Advancing particle" << std::endl;

    switch (direction)
    {
    case 1:
      if(xp[m]<low) 
        xp[m] = low;
      else if(xp[m]>high)
        xp[m] = high;
      break;
    case 2:
      if(xp[m]<low) 
        xp[m] = low;
      else if(xp[m]>high)
        xp[m] = high;
      break;
    case 3:
      if(xp[m]<low) 
        xp[m] = low;
      else if(xp[m]>high)
        xp[m] = high;
      break;
    }

  }
  
}
void rkparticlesubstep(long double a1, long double a2, long double a3,
                        std::vector<std::vector<std::vector<long double>>> &un,
                        std::vector<std::vector<std::vector<long double>>> &vn,
                        std::vector<std::vector<std::vector<long double>>> &wn,
                        std::vector<long double> &ku, std::vector<long double> &kv, std::vector<long double> &kw, 
                        const long double delT, size_t np,
                        std::vector<long double> &xe, std::vector<long double> &ye, std::vector<long double> &ze,
                        std::vector<long double> &xc, std::vector<long double> &yc, std::vector<long double> &zc,
                        std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp,
                        size_t nxp2, size_t nyp2, size_t nzp2)
{
  std::vector<long double> iu;
  std::vector<long double> iv;
  std::vector<long double> iw;

  iu = std::vector<long double>(np, 0.0);
  iv = std::vector<long double>(np, 0.0);
  iw = std::vector<long double>(np, 0.0);

  long double lowX = xe[0], highX = xe[nxp2-2];
  long double lowY = ye[0], highY = ye[nyp2-2];
  long double lowZ = ze[0], highZ = ze[nzp2-2];

  interpolate(iu,un,xp,yp,zp,xe,yc,zc,np);
  updateK(a1,ku,iu,delT,np);

  interpolate(iv,vn,xp,yp,zp,xc,ye,zc,np);
  updateK(a1,kv,iv,delT,np);
  
  interpolate(iw,wn,xp,yp,zp,xc,yc,ze,np);
  updateK(a1,kw,iw,delT,np);
  
  advance_pa(a2,ku,xp,np,lowX,highX,1);
  advance_pa(a2,kv,yp,np,lowY,highY,2);
  advance_pa(a2,kw,zp,np,lowZ,highZ,3);


}


void NSSolverRK3_CN(std::vector<std::vector<std::vector<long double>>> &un,
                    std::vector<std::vector<std::vector<long double>>> &vn,
                    std::vector<std::vector<std::vector<long double>>> &wn,
                    std::vector<std::vector<std::vector<long double>>> &pn,
                    size_t nx, size_t ny, size_t nz, size_t np,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                    long double dx, long double dy, long double dz, long double Time,const long double initialTimestep, long double Re,
                    std::vector<long double> &xe, std::vector<long double> &ye, std::vector<long double> &ze,
                    std::vector<long double> &xc, std::vector<long double> &yc, std::vector<long double> &zc,
                    std::vector<long double> &xp, std::vector<long double> &yp, std::vector<long double> &zp)
{
  int count = 0;
  long double totaltime = 0.0;
  long double delT = 0.0;
  while(totaltime<Time)
  {
  applyFlowVelocityBC(un,vn,wn,pn,nx,ny,nz,nxp2,nyp2,nzp2);
  
  if(count == 0)
    delT = initialTimestep;
  else
    delT = calculateDelt(un,vn,wn,nxp2,nyp2,nzp2,dx,dy,dz);
  totaltime = totaltime + delT;
  printlevel(count,totaltime,delT);
  printParticle(xp,yp,zp,np,count);

  std::vector<std::vector<std::vector<long double>>> qu;
  std::vector<std::vector<std::vector<long double>>> qv;
  std::vector<std::vector<std::vector<long double>>> qw;

  qu = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qv = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qw = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  
  std::vector<long double> ku;
  std::vector<long double> kv;
  std::vector<long double> kw;

  ku = std::vector<long double>(np, 0.0);
  kv = std::vector<long double>(np, 0.0);
  kw = std::vector<long double>(np, 0.0);

  long double a1 = 0.0, a2 = 0.0, a3 = 0.0, dt_rk = 0.0;

  a1 = 0.0;
  a2 = 1.0;
  a3 = 1.0 / 3.0;
  dt_rk = 1.0 / 3.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  
  //if (count>80)
  //{
    a1 = 0.0;
    a2 = 1.0 / 3.0;
    a3 = 0.0;
    rkparticlesubstep(a1,a2,a3,un,vn,wn,ku,kv,kw,delT,np,xe,ye,ze,xc,yc,zc,xp,yp,zp,nxp2,nyp2,nzp2);
  //}
  
  

  a1 = -5.0/9.0;
  a2 = 15.0/16.0;
  a3 = 5.0/12.0;
  dt_rk = 5.0/12.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  
  //if (count>80)
  //{
    a1 = -5.0/9.0;
    a2 = 15.0/16.0;
    a3 = 1.0/3.0;
    rkparticlesubstep(a1,a2,a3,un,vn,wn,ku,kv,kw,delT,np,xe,ye,ze,xc,yc,zc,xp,yp,zp,nxp2,nyp2,nzp2);
  //}
  

  a1 = -153.0/128.0;
  a2 = 8.0/15.0;
  a3 = 1.0/4.0;
  dt_rk = 1.0 / 4.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  
  //if (count>5)
  //{
    a1 = -153.0/128.0;
    a2 = 8.0/15.0;
    a3 = 3.0/4.0;
    rkparticlesubstep(a1,a2,a3,un,vn,wn,ku,kv,kw,delT,np,xe,ye,ze,xc,yc,zc,xp,yp,zp,nxp2,nyp2,nzp2);
  //}
  

  std::vector<std::vector<std::vector<long double>>>().swap(qu);
  std::vector<std::vector<std::vector<long double>>>().swap(qv);
  std::vector<std::vector<std::vector<long double>>>().swap(qw);

  //if((count%5)==0)
  print(1.0,1.0,1.0,dx,dy,dz,nx,ny,nz,nxp2,nyp2,nzp2,un,vn,wn,pn,count,totaltime,delT);




  count++;
  }

  
    
}
/*void NSSolverRK3_CN(std::vector<std::vector<std::vector<long double>>> &un,
                    std::vector<std::vector<std::vector<long double>>> &vn,
                    std::vector<std::vector<std::vector<long double>>> &wn,
                    std::vector<std::vector<std::vector<long double>>> &pn,
                    size_t nx, size_t ny, size_t nz,
                    size_t nxp2, size_t nyp2, size_t nzp2,
                   long double dx, long double dy, long double dz, long double Time,const long double initialTimestep, long double Re)
{
  int count = 0;
  long double totaltime = 0.0;
  long double delT = 0.0;

  while(totaltime<Time)
  {
  applyFlowVelocityBC(un,vn,wn,pn,nx,ny,nz,nxp2,nyp2,nzp2);
  
  if(count == 0)
    delT = initialTimestep;
  else
    delT = calculateDelt(un,vn,wn,nxp2,nyp2,nzp2,dx,dy,dz);
  totaltime = totaltime + delT;
  printlevel(count,totaltime,delT);

  std::vector<std::vector<std::vector<long double>>> qu;
  std::vector<std::vector<std::vector<long double>>> qv;
  std::vector<std::vector<std::vector<long double>>> qw;

  qu = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qv = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
  qw = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

  long double a1 = 0.0, a2 = 0.0, a3 = 0.0, dt_rk = 0.0;

  a1 = 0.0;
  a2 = 1.0;
  a3 = 1.0 / 3.0;
  dt_rk = 1.0 / 3.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  a1 = -5.0/9.0;
  a2 = 15.0/16.0;
  a3 = 5.0/12.0;
  dt_rk = 5.0/12.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  a1 = -153.0/128.0;
  a2 = 8.0/15.0;
  a3 = 1.0/4.0;
  dt_rk = 1.0 / 4.0 * delT;
  rksubstep(a1, a2, a3, dt_rk, delT, un, vn, wn, pn, qu, qv, qw, nx, ny, nz, nxp2, nyp2, nzp2, dx, dy, dz, Re);

  std::vector<std::vector<std::vector<long double>>>().swap(qu);
  std::vector<std::vector<std::vector<long double>>>().swap(qv);
  std::vector<std::vector<std::vector<long double>>>().swap(qw);

  if((count%5)==0)
  print(1.0,1.0,1.0,dx,dy,dz,nx,ny,nz,nxp2,nyp2,nzp2,un,vn,wn,pn,count,totaltime,delT);




  count++;
  }

  
    
}
*/