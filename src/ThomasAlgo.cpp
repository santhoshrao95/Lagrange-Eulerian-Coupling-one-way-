#include"ThomasAlgo.h"

void thomas_rhs(std::vector<std::vector<std::vector<long double>>> &VTemp,
                    std::vector<std::vector<std::vector<long double>>> &m_Vn,
                    std::vector<std::vector<std::vector<long double>>> &FTemp,
                    size_t nxp2, size_t nyp2, size_t nzp2)
{
  for (size_t k = 0; k < nzp2; k++)
  {
    for (size_t j = 0; j < nyp2; j++)
    {
      for (size_t i = 0; i < nxp2; i++)
      {
        FTemp[i][j][k] = m_Vn[i][j][k] + VTemp[i][j][k];
      }      
    }    
  }
}
void thomasSolver(std::vector<std::vector<std::vector<long double>>> &uTemp1,
                    std::vector<std::vector<std::vector<long double>>> &vTemp1,
                    std::vector<std::vector<std::vector<long double>>> &wTemp1,
                    std::vector<std::vector<std::vector<long double>>> &rTemp1,
                    std::vector<std::vector<std::vector<long double>>> &rTemp2,
                    std::vector<std::vector<std::vector<long double>>> &rTemp3,
                    size_t nx, size_t ny, size_t nz,
                    const long double dx, const long double dy, const long double dz,
                    const long double Re, const long double dtrk, size_t n)
{
    const long double rRe = 1.0/Re, rdx = 1.0/dx, rdy = 1.0/dy, rdz = 1.0/dz;
    size_t nxp2 = nx+2, nyp2 = ny+2, nzp2 = nz+2;
    std::vector<long double> a;
    std::vector<long double> b;
    std::vector<long double> c;
    a = std::vector<long double>(nzp2,0.0);
    b = std::vector<long double>(nzp2,0.0);
    c = std::vector<long double>(nzp2,0.0);

    switch (n)  
    {
    case 1:
    {
    for (size_t i = 1; i < nxp2-2; i++)
    {
        for (size_t j = 1; j < nyp2-1; j++)
        {
            for (size_t k = 0; k < nzp2; k++)
            {
                a[k] = - dtrk * rRe * rdz * rdz;
                b[k] = 1.0 + (2.0 * dtrk * rRe * rdz * rdz);
                c[k] = - dtrk * rRe * rdz * rdz;
            }

            //for dirichilet bc
            a[0] = 0.0;
            b[0] = 1.0;
            c[0] = 0.0;
            a[nzp2-1] = 0.0;
            b[nzp2-1] = 1.0;
            c[nzp2-1] = 0.0;

            for (size_t k = 1; k < nzp2; k++)
            {
                long double factor = a[k]/b[k-1];
                b[k] = b[k] - c[k-1] * factor;
                rTemp1[i][j][k] = rTemp1[i][j][k] - rTemp1[i][j][k-1] * factor;  
            }

            uTemp1[i][j][nzp2-1] = rTemp1[i][j][nzp2-1] / b[nzp2-1];

            for (size_t k = (nzp2-2); k > 0; k--)
            {
                uTemp1[i][j][k] = (rTemp1[i][j][k] - (c[k] * uTemp1[i][j][k+1])) / b[k];
            }
            
        }
        
    }
    }
    break;

    case 2:
    {
        for (size_t i = 1; i < nxp2-1; i++)
        {
            for (size_t j = 1; j < nyp2-2; j++)
            {
                for (size_t k = 0; k < nzp2; k++)
                {
                    a[k] = - dtrk * rRe * rdz * rdz;
                    b[k] = 1.0 + (2.0 * dtrk * rRe * rdz * rdz);
                    c[k] = - dtrk * rRe * rdz * rdz;
                }

                //for dirichilet bc
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = 0.0;
                a[nzp2-1] = 0.0;
                b[nzp2-1] = 1.0;
                c[nzp2-1] = 0.0;

                for (size_t k = 1; k < nzp2; k++)
                {   
                    long double factor = a[k]/b[k-1];
                    b[k] = b[k] - c[k-1] * factor;
                    rTemp2[i][j][k] = rTemp2[i][j][k] - rTemp2[i][j][k-1] * factor;  
                }

                vTemp1[i][j][nzp2-1] = rTemp2[i][j][nzp2-1] / b[nzp2-1];

                for (size_t k = (nzp2-2); k > 0; k--)
                {
                    vTemp1[i][j][k] = (rTemp2[i][j][k] - (c[k] * vTemp1[i][j][k+1])) / b[k];
                }
            
            }
        
        }
    }
    break;
    case 3:
    {
        for (size_t i = 1; i < nxp2-1; i++)
        {
            for (size_t j = 1; j < nyp2-2; j++)
            {
                for (size_t k = 0; k < nzp2-1; k++)
                {
                    a[k] = - dtrk * rRe * rdz * rdz;
                    b[k] = 1.0 + (2.0 * dtrk * rRe * rdz * rdz);
                    c[k] = - dtrk * rRe * rdz * rdz;
                }

                //for dirichilet bc
                a[0] = 0.0;
                b[0] = 1.0;
                c[0] = 0.0;
                a[nzp2-2] = 0.0;
                b[nzp2-2] = 1.0;
                c[nzp2-2] = 0.0;

                for (size_t k = 1; k < nzp2-1; k++)
                {   
                    long double factor = a[k]/b[k-1];
                    b[k] = b[k] - c[k-1] * factor;
                    rTemp3[i][j][k] = rTemp3[i][j][k] - rTemp3[i][j][k-1] * factor;  
                }

                wTemp1[i][j][nzp2-2] = rTemp3[i][j][nzp2-2] / b[nzp2-2];

                for (size_t k = (nzp2-3); k > 0; k--)
                {
                    wTemp1[i][j][k] = (rTemp3[i][j][k] - (c[k] * wTemp1[i][j][k+1])) / b[k];
                }
            
            }
        
        }
    }
    break;
    default:
    {
        std::cout << "I'm out!" << std::endl;
    }
    break;
    }
    
}