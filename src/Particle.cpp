#include"Particle.h"

Particle::Particle(size_t nxp2, size_t nyp2, size_t nzp2, size_t np, const long double d,
        const long double lx, const long double ly, const long double lz):
        nxp2(nxp2),nyp2(nyp2),nzp2(nzp2),np(np),d(d),lx(lx),ly(ly),lz(lz){
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

    std::vector<double> xp(np, 0.0),yp(np, 0.0),zp(np, 0.0);

    for (size_t i = 0; i < np; i++)
    {
        xp[i] = xc[i+20];
        yp[i] = yc[i+20];
        zp[i] = zc[i+20];
    }
    

    
    
}


/*s mul(const double a, const s& v){
    return s(a*v.x,a*v.y,a*v.z);
}
s operator * (const double a, const s& v){
    return mul(a,v);
}*/