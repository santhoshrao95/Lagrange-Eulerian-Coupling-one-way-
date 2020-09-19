#pragma once
#include<vector>
/*struct s{
        long double x;
        long double y;
        long double z;

        s() {}
        s(long double x, long double y, long double z):x(x),y(y),z(z) {}
    };*/
class Particle{
    
    public:
    
    
    //std::vector<s> p,k;
    Particle(size_t nxp2, size_t nyp2, size_t nzp2, size_t np, const long double d,
        const long double lx, const long double ly, const long double lz);
    Particle();
    size_t nx,ny,nz,nxp2,nyp2,nzp2,np;
    const long double lx,ly,lz,d; 
    std::vector<long double> xe,ye,ze,xc,yc,zc,xp,yp,zp;

    
};

//s mul(const double a, const s& v);
//s operator * (const double a, const s& v);