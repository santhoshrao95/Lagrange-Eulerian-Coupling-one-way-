#include "output.h"

void printlevel(size_t count, const long double totaltime, const long double delT)
{
    std::ofstream output("whereiam", std::ios::app);
    output<<"I have entered count: " << count << " " << totaltime << " " << delT << std::endl;
    output.close();
}
void print(const long double lx, const long double ly, const long double lz,
            const long double dx, const long double dy, const long double dz,
            size_t nx, size_t ny, size_t nz,
            size_t nxp2, size_t nyp2, size_t nzp2, 
            std::vector<std::vector<std::vector<long double>>> &u,
            std::vector<std::vector<std::vector<long double>>> &v,
            std::vector<std::vector<std::vector<long double>>> &w,
            std::vector<std::vector<std::vector<long double>>> &p,
            size_t count, const long double Time, const long double delT)
{
  const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
  std::vector<double> xe(nxp2, 0.0),ye(nxp2, 0.0),ze(nxp2, 0.0),xc(nxp2, 0.0),yc(nxp2, 0.0),zc(nxp2, 0.0);
    for (size_t i = 0; i < nxp2; i++)
    {
        xe[i] = x0 + lx/double(nxp2-1) * double(i+1) - 0.5 * lx / double(nxp2-1);
        ye[i] = x0 + ly/double(nyp2-1) * double(i+1) - 0.5 * ly / double(nyp2-1);
        ze[i] = z0 + lz/double(nzp2-1) * double(i+1) - 0.5 * lz / double(nzp2-1);
    }
    
    for (size_t i = 1; i < nxp2; i++)
    {
        xc[i] = (xe[i]+xe[i-1]) / 2.0;
        yc[i] = (ye[i]+ye[i-1]) / 2.0;
        zc[i] = (ze[i]+ze[i-1]) / 2.0;
    }
    
    xc[0] = xe[0] - 0.5 * (xe[1]-xe[0]);
    yc[0] = ye[0] - 0.5 * (ye[1]-ye[0]);
    zc[0] = ze[0] - 0.5 * (ze[1]-ze[0]);


    std::string init("output");
    std::string add(std::to_string(count));
    size_t num =0, temp=count;
    while(temp!=0)
    {
      num++;
      temp = temp / 10;
    }
    if (num == 1)
    {
      init.append(".csv.000");
      init.append(add);
    }
    else if(num==2)
    {
      init.append(".csv.00");
      init.append(add);
    }
    else if (num==3)
    {
      init.append(".csv.0");
      init.append(add);
    }
    else
    {
      init.append(".csv.");
      init.append(add);
    }
    
        


/*std::string init("output");
    std::string add(std::to_string(count));
    //init.append(add);
    init.append(".csv.");
    init.append(add);
    //init.append(std::to_string(Time));*/


    std::ofstream output(init, std::ios::trunc);
    output << "x " <<"y " << "z " << "u " <<"v " << "w " << "p " <<std::endl; 
    for (size_t k = 0; k < nzp2; k++)
    {
        for (size_t j = 0; j < nyp2; j++)
        {
            for (size_t i = 0; i < nxp2; i++)
            {
                output << xc[i] << " " << yc[j] << " " << zc[k] << " "  << u[i][j][k] << " " << v[i][j][k] << " " << w[i][j][k] << " " << p[i][j][k] << " " << std::endl;  
            }
            
        }
        
    }
    output.close();

}

void printParticle(std::vector<long double> &xp,std::vector<long double> &yp,std::vector<long double> &zp,
                    size_t np, size_t count)
{
    /*std::string init("PatrticleOutput");
    std::string add(std::to_string(count));
    init.append(".csv.");
    init.append(add);*/


    std::string init("particleoutput");
    std::string add(std::to_string(count));
    size_t num =0, temp=count;
    while(temp!=0)
    {
      num++;
      temp = temp / 10;
    }
    if (num == 1)
    {
      init.append(".csv.000");
      init.append(add);
    }
    else if(num==2)
    {
      init.append(".csv.00");
      init.append(add);
    }
    else if (num==3)
    {
      init.append(".csv.0");
      init.append(add);
    }
    else
    {
      init.append(".csv.");
      init.append(add);
    }
    
    

    std::ofstream output(init, std::ios::trunc);
    for(size_t i=0;i<np;i++)
    {
        output << xp[i] << " " << yp[i] <<" " << zp[i] << std::endl;
    }
}
