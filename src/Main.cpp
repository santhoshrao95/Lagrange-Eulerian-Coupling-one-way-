#include <math.h>
#include <sstream>
#include <vector>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#include "CavityFlow.h"
int main(int argc, char* argv[])
{
    CavityFlow flow(1.0,1.0,1.0,120,0.01818,400.0,0.000001,32,32,32,1024);

    flow.run();

    return 0;
}