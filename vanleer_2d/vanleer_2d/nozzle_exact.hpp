//
//  nozzle_exact.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef nozzle_exact_hpp
#define nozzle_exact_hpp

#include <stdio.h>
#include <fstream>

class nozzle_exact {
    
    
private:
   //private variables that has been used in the main fn has _Var to distinquish it from the public variable in the main
    double _A_x , _P_0, _T_0, _Th_area, T, P, rho, u; //initialize area, reference (pressure, temp), throttle area, temperature(x), pressure(x), rho(x), u(x)
    double M; //mach number
    double epsi; //to make calculations easier
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double phi; //to make calculations easier
    double A_par; //A/A^*
    double F_M; //F(M) M is the mach number
    double dF_M; //d(F(M)/dM
    int iter = 0; //variable to track num of iterations
    
    std::fstream myfile; //creation of file to write the results
   
public:
    //constructor to open the file for write
    nozzle_exact(){
        
        myfile.open("/Users/cringedaddy/CFD class/project_roe/hw_3_roe/hw_3_roe/cfd_project_roe/cfd_project_roe/results/uexact.csv", std::ios::trunc | std::ios::out);
        myfile<<"#########x########   "<<"#########u_exact########"<<std::endl;
    };
    // destructor to close the file
    ~nozzle_exact()
    {
        myfile.close();
        
    };
    //initialization functions that takes Area(x),pressure, temperature, throttle area
    void initialization(double, double, double, double);
    //to run iterations
    void run_simulation(double);
    //to print the results in file
    void print_(double);
   
    
    
};


#include <cmath>
using std::acos;
using std::cos;
using std::sin;

// Set the real precision for the program
using real_t = double;

// Define constants module
namespace constants {

    const real_t one = 1.0;
    const real_t two = 2.0;
    const real_t three = 3.0;
    const real_t four = 4.0;
    const real_t five = 5.0;
    const real_t six = 6.0;
}

// Define set_precision module
namespace set_precision {

    const int sngl = std::numeric_limits<real_t>::digits10;
    const int dbl = 2 * sngl;
    const int dp = dbl;
}

// Define set_inputs module
namespace set_inputs {

    using namespace constants;
    using namespace set_precision;

    real_t pi;

    void initialize_constants() {
        pi = acos(-one);
    }
}

// Define MMS_constants module
namespace MMS_constants {

    using namespace constants;
    using namespace set_precision;

    const real_t rho0 = 1.0;
    const real_t rhox = 0.15;
    const real_t rhoy = -0.1;
    const real_t uvel0 = 800.0;
    const real_t uvelx = 50.0;
    const real_t uvely = -30.0;
    const real_t vvel0 = 800.0;
    const real_t vvelx = -75.0;
    const real_t vvely = 40.0;
    const real_t wvel0 = 0.0;
    const real_t wvelx = 0.0;
    const real_t wvely = 0.0;
    const real_t press0 = 100000.0;
    const real_t pressx = 20000.0;
    const real_t pressy = 50000.0;
}



#endif /* nozzle_exact_hpp */


