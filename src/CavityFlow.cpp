#include"CavityFlow.h"
#include<vector>

//default constructor
CavityFlow::CavityFlow():Lx(1.0),Ly(1.0),Lz(1.0),timeLength(10.0), initialTimestep(0.001),
                                           Re(200.0),nx(32),ny(32),nz(32),dia(0.000001),
                                           nxp2(nx+2),nyp2(ny+2),nzp2(nz+2)
{
    u = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    v = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    w = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    p = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));

}

//parameterized constructor
CavityFlow::CavityFlow(const long double xL,const long double yL,const long double zL,
                                         const long double timeLength, const long double initialTimestep, const long double Re, const long double dia,
                                         size_t nx, size_t ny, size_t nz, size_t np):
                                         Lx(xL),Ly(yL),Lz(zL),timeLength(timeLength),
                                         initialTimestep(initialTimestep),Re(Re),dia(dia),
                                         nx(nx),ny(ny),nz(nz),np(np),nxp2(nx+2),nyp2(ny+2),nzp2(nz+2)
{
    u = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    v = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    w = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    p = std::vector<std::vector<std::vector<long double>>>(nxp2, std::vector<std::vector<long double>>(nyp2, std::vector<long double>(nzp2, 0.0)));
    xe = std::vector<long double>(nxp2,0.0);
    ye = std::vector<long double>(nyp2,0.0);
    ze = std::vector<long double>(nzp2,0.0);
    xc = std::vector<long double>(nxp2,0.0);
    yc = std::vector<long double>(nyp2,0.0);
    zc = std::vector<long double>(nzp2,0.0);
    xp = std::vector<long double>(np,0.0);
    yp = std::vector<long double>(np,0.0);
    zp = std::vector<long double>(np,0.0);

    const double x0 = 0.0, y0 = 0.0, z0 = 0.0;
    
    for (size_t i = 0; i < nxp2; i++)
    {
        xe[i] = x0 + Lx/double(nxp2-1) * double(i+1) - 0.5 * Lx / double(nxp2-1);
         ye[i] = x0 + Ly/double(nyp2-1) * double(i+1) - 0.5 * Ly / double(nyp2-1);
        ze[i] = z0 + Lz/double(nzp2-1) * double(i+1) - 0.5 * Lz / double(nzp2-1);
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

    

    /*for (size_t i = 0; i < np; i++)
    {
        xp[i] = xc[10];
        yp[i] = yc[10];
        zp[i] = zc[i+10];
    }*/

    int len = sqrt(np);

    for (size_t i = 0; i < len; i++)
    {
        for (size_t j = 0; j < len; j++)
        {
            xp[i*len + j] = xc[j+1];
            yp[i*len + j] = yc[len/2];
            zp[i*len + j] = zc[i+1];
        }
        
    }
    

    
    
    
}
void CavityFlow::run()
{
    //Particle pa(66,66,66,10,0.000001,1.0,1.0,1.0);
    NSSolverRK3_CN(u,v,w,p,nx,ny,nz,np,nxp2,nyp2,nzp2,dx,dy,dz,timeLength,initialTimestep,Re,xe,ye,ze,xc,yc,zc,xp,yp,zp);
}