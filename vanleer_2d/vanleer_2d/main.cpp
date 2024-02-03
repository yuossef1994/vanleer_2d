//
//  main.cpp
//  homework_3
//
//  Created by Youssef Z on 2/26/23.
//

#include <iostream>
#include "Dv_nozzle.hpp"

#include <math.h>

int main(int argc, const char * argv[]) {
    
    
 
    
   
   
    
    
    
    
    
    int max_iter=600000;
    double L2norm;
    
    Dv_nozzle nozzle;
    
    nozzle.mesh_nodes();
    
    nozzle.initialization(300000, 600);
    
    nozzle.set_boundary_cond();
    nozzle.time_step();
    nozzle.euler_explicit();
    nozzle.rL2initial();
    nozzle.L2norm(1);
    
    for(int i=1;i<=max_iter;i++)
    {
        
       
        nozzle.set_boundary_cond();
        nozzle.time_step();
        nozzle.euler_explicit();
        
        L2norm= nozzle.L2norm(i);
        
        if(i%10==0)std::cout<< " iter number  "<<i <<"  L2 norm is  " <<L2norm<<std::endl;
        
        
        
        if(L2norm < 1e-8){
            
            std::cout<< "solution converged"<<std::endl;
            break;
        }
       
    }
    
    
    nozzle.print_res();
    nozzle.DE();
   // std::cout<<"h is = "<<abs((x_f[1]-x_f[0]))<<std::endl;
  
    return 0;
}
