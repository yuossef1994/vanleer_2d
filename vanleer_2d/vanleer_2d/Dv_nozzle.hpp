//
//  Dv_nozzle.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef Dv_nozzle_hpp
#define Dv_nozzle_hpp

#include <stdio.h>
#include <fstream>
#include <fstream>
#include <iostream>
#include<string>
#include <cmath>
#include <limits>

class Dv_nozzle {
    
    
private:

    
  
    
    
    double*** epsi_plus_y;
    double*** epsi_plus_x;
    double*** epsi_minus_y;
    double*** epsi_minus_x;
    double*** _V;
    double***  _U;
    double*** _F;
    double*** _F_eta;
    static int const _imax_f=65, _imax_c=64, _imin=0;
    static int const  i_n_max=65, j_n_max=65, k_n_max=2;
    int im_f, im_c;
    
    double A, B, C, D, E, F;  // just variables for calculations
    double k_2=0.2, k_4=0.02;
    double x_f[_imax_f];
    double y_f[_imax_f];
    double epsi_plus_x_ghost[4][_imax_f+2];
    double epsi_minus_x_ghost[4][_imax_f+2];
    double r_1[4],r_2[4],r_3[4],r_4[4];
    double lambda[4];
    double eigen_max[2];
    double R;
    double rho_roe, u_roe,v_roe, ht_roe, a_roe;
    double T_vel_L,T_vel_R,T_vel_D,T_vel_U,T_vel_roe;
    double nx,ny;
    
    double delta_u,delta_v, delta_P, delta_rho;
    //variables for grid nodes and geometry terms
    
    double x_n[i_n_max][j_n_max];
    double y_n[i_n_max][j_n_max];
  
    
   
   
    
   
    
    
    
    
    
    
    double _F_P_kasi[4] , _F_C_kasi[4];
    double _F_P_eta[4] , _F_C_eta[4];
    
    
    
    
    double delta_w[4];
    
    double avg_area;
    double avg_d_area;
    double x_c[_imax_c][_imax_c];
    double y_c[_imax_c][_imax_c];
    double delta_t[_imax_c][_imax_c];
    double delta_t_gl=999;
    double CFL=0.8;
    double area_f[_imax_f][_imax_f];
    double _area, _d_area;
    double mu[_imax_c];
    double epsilon_2;
    double _R[4]={0,0,0,0};
    double _R_iter[4]={0,0,0,0};
    double eigen_v_avg;
    double _delta_x, _delta_y;
    double _P_0,  _T_0, mach_in, mach_out, T_in, u_in, P_in, rho_in, et_in, T_out, u_out, P_out, rho_out, et_out, v_in,v_out;
    double _M_0;
    double _mach[_imax_c][_imax_c];
    double T[_imax_c][_imax_c];
    double et[_imax_c][_imax_c];
    double ht[_imax_c][_imax_c];
    
   
    double _U_ghost_outflow[4][_imax_c];
    
    double _V_ghost_wall_0[6][_imax_c];
    
    double _V_ghost_wall_jmax[6][_imax_c];
    
    double _V_ghost_outflow[6][_imax_c];
    double _V_ghost_outflow_2[6][_imax_c];
    double _V_ghost_outflow_3[6][_imax_c];
    double _U_ghost_inflow[4][_imax_c];
    double _V_ghost_inflow[6][_imax_c];
    double _V_ghost_inflow_2[6][_imax_c];
    double _V_ghost_inflow_3[6][_imax_c];
    double  _V_ghost_wall_0_2[6][_imax_c];
    double  _V_ghost_wall_0_3[6][_imax_c];
    
    double  _V_ghost_wall_jmax_2[6][_imax_c];
    double  _V_ghost_wall_jmax_3[6][_imax_c];
    
    //upwind schemes variables
    
    double kappa=-1;
    double epsilon_upwind=1;
    
    
    double c_plus;
    double c_minus;
    double alpha_plus;
    double alpha_minus;
    double beta_R;
    double beta_L;
    double beta_U;
    double beta_D;
    double a[_imax_c][_imax_c];
    double a_L;
    double a_R;
    double a_D;
    double a_U;
    
    double D_plus;
    double D_minus;
    double P_plus;
    double P_minus;
    double V_L[4];
    double V_R[4];
    double V_D[4];
    double V_U[4];
    double T_L;
    double T_R;
    double T_D;
    double T_U;
    double r_plus;
    double r_minus;
    
    double ht_L,ht_R,ht_D,ht_U;
    
    
    double mach_knight;
    double M_plus;
    double M_minus;
    double M_L;
    double M_R;
    double M_D;
    double M_U;
    
    
    
    
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double epsi_1;
    double rL2norm[4];
    double _rL2initial[4]={0,0,0,0};
    std::fstream L_vs_Iter; //creation of file to write the results
    
    std::fstream rho_vs_x; //creation of file to write the results
    std::fstream u_vs_x; //creation of file to write the results
    std::fstream v_vs_x; //creation of file to write the results
    std::fstream p_vs_x; //creation of file to write the results
    
public:
    
    
    // Allocate memory for the rows
    
    //constructor to open the file for write
    Dv_nozzle(){
        
        
        
        rho_vs_x.open("/Users/cringedaddy/CFD class/vanleer/vanleer_2d/vanleer_2d/results/rho.csv", std::ios::trunc | std::ios::out);
        
        u_vs_x.open("/Users/cringedaddy/CFD class/vanleer/vanleer_2d/vanleer_2d/results/u.csv", std::ios::trunc | std::ios::out);
        v_vs_x.open("/Users/cringedaddy/CFD class/vanleer/vanleer_2d/vanleer_2d/results/v.csv", std::ios::trunc | std::ios::out);
        p_vs_x.open("/Users/cringedaddy/CFD class/vanleer/vanleer_2d/vanleer_2d/results/p.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter.open("/Users/cringedaddy/CFD class/vanleer/vanleer_2d/vanleer_2d/results/L2.csv", std::ios::trunc | std::ios::out);
        L_vs_Iter<<"#iter#"<<"##L2_1#"<<"##L2_2#"<<"##L2_3#"<<std::endl;
        
    };
    // destructor to close the file
    ~Dv_nozzle()
    {
        
        
     
        for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < _imax_f+2; j++) {
                        delete[] epsi_plus_y[i][j];
                        delete[] epsi_plus_x[i][j];
                        delete[] epsi_minus_y[i][j];
                        delete[] epsi_minus_x[i][j];
                        delete[] _V[i][j];
                        delete[]  _U[i][j];
                        delete[] _F[i][j];
                        delete[] _F_eta[i][j];
                    }
                       delete[] epsi_plus_y[i];
                       delete[] epsi_plus_x[i];
                       delete[] epsi_minus_y[i];
                       delete[] epsi_minus_x[i];
                       delete[] _V[i];
                       delete[]  _U[i];
                       delete[] _F[i];
                       delete[] _F_eta[i];
                }
                delete[] epsi_plus_y;
             
                delete[] epsi_plus_x;
                delete[] epsi_minus_y;
                delete[] epsi_minus_x;
                delete[] _V;
                delete[]  _U;
                delete[] _F;
                delete[] _F_eta;
      
        
        
        
        
       
        
        
        
        rho_vs_x.close();
        u_vs_x.close();
        p_vs_x.close();
        L_vs_Iter.close();
        
    };
    
public:
    //void set_geometry(int,double,double);
    void initialization(double,double);
    void set_boundary_cond();
    void euler_explicit();
    double area(double);
    double d_area(double);
    void time_step();
    void print_res();
    void rL2initial();
    double L2norm(int);
    void vanleer();
    void vanleer_kasi(int,int);
    
    void vanleer_eta(int ,int);

    void mesh_nodes();
    // to calculate face area, cell volume (for the eigen values) and outpointing normals
    double face_area(double,double,double,double);
    double n_x(double, double,double);
    double n_y(double,double,double);
    double cell_vol(int,int);
    double cell_vol_eta(double,double);
    double cell_vol_kasi(double,double);
    double  Delta_ghost_x_y;
    double  y_c_ghost;
    double x_c_ghost;
    void DE();

private:
    const int sngl = std::numeric_limits<float>::digits10;
    const int dbl = std::numeric_limits<double>::digits10;
    const int dp = dbl;

    const double pi = std::acos(-1.0);
    const double Pi = std::acos(-1.0);
    const double rho0 = 1.0;
    const double rhox = 0.15;
    const double rhoy = -0.1;
    const double uvel0 = 800.0;
    const double uvelx = 50.0;
    const double uvely = -30.0;
    const double vvel0 = 800.0;
    const double vvelx = -75.0;
    const double vvely = 40;
    const double wvel0 = 0.0;
    const double wvelx = 0.0;
    const double wvely = 0.0;
    const double press0 = 100000.0;
    const double pressx = 20000.0;
    const double pressy = 50000.0;
    const double one = 1.0;
    const double two = 2.0;
    const double three = 3.0;
    const double four = 4.0;
    const double five = 5.0;
    const double six = 6.0;

public:
    double rho(double length, double x, double y) {
        return rho0 + rhoy * std::cos((pi*y) / (2.*length)) + rhox * std::sin((pi*x) / length);
    }

    double uvel(double length, double x, double y) {
        return uvel0 + uvely * std::cos((3.*pi*y) / (5.*length)) + uvelx * std::sin((3.*pi*x) / (2.*length));
    }

    double vvel(double length, double x, double y) {
        return vvel0 + vvelx * std::cos((pi*x) / (2.*length)) + vvely * std::sin((2.*pi*y) / (3.*length));
    }

    double press(double length, double x, double y) {
        return press0 + pressx * std::cos((2.*pi*x) / length) + pressy * std::sin((pi*y) / length);
    }

    double rmassconv(double length, double x, double y) {
            double result_mass= (3*pi*uvelx*cos((3*pi*x)/(2.*length))*(rho0 + rhoy*cos((pi*y)/(2.*length)) + rhox*sin((pi*x)/length)))/(2.*length) +
           (2*pi*vvely*cos((2*pi*y)/(3.*length))*(rho0 + rhoy*cos((pi*y)/(2.*length)) + rhox*sin((pi*x)/length)))/(3.*length) +
           (pi*rhox*cos((pi*x)/length)*(uvel0 + uvely*cos((3*pi*y)/(5.*length)) + uvelx*sin((3*pi*x)/(2.*length))))/length -
        (pi*rhoy*sin((pi*y)/(2.*length))*(vvel0 + vvelx*cos((pi*x)/(2.*length)) + vvely*sin((2*pi*y)/(3.*length))))/(2.*length);
        
        return result_mass;
        /*
         
         (three*pi*uvelx*cos((three*pi*x)/(two*length)) *
         (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /
         (two*length) + (two*pi*vvely*cos((two*pi*y)/(three*length)) *
         (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length))) /
         (three*length) + (pi*rhox*cos((pi*x)/length) *
         (uvel0 + uvely*cos((three*pi*y)/(five*length)) + uvelx*sin((three*pi*x)/
         (two*length))))/length - (pi*rhoy*sin((pi*y)/(two*length)) *
         (vvel0 + vvelx*cos((pi*x)/(two*length)) + vvely*sin((two*pi*y) /
         (three*length))))/(two*length);
         */
        
        
        
        
        
        }
    
    
    


double xmtmconv(double L, double x, double y) {
    

    double result =  (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*
                       (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                            rhox*sin((Pi*x)/L))*
                           (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                             uvelx*sin((3*Pi*x)/(2.*L))))/L +
                        (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*
                           (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                             rhox*sin((Pi*x)/L))*
                          (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                             uvelx*sin((3*Pi*x)/(2.*L))))/(3.*L) +
    (pi*rhox*cos((pi*x)/L)*pow((uvel0 + uvely*cos((three*pi*y) / (five*L)) +uvelx*sin((three*pi*x)/(two*L))),2))/L -
                        (2*Pi*pressx*sin((2*Pi*x)/L))/L -
                        (Pi*rhoy*(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                             uvelx*sin((3*Pi*x)/(2.*L)))*
                           sin((Pi*y)/(2.*L))*
                           (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                            vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) -
                        (3*Pi*uvely*(rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                             rhox*sin((Pi*x)/L))*sin((3*Pi*y)/(5.*L))*
                          (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                           vvely*sin((2*Pi*y)/(3.*L))))/(5.*L);/*(three * pi * uvelx * std::cos((three * pi * x) / (two * length)) *
            (rho0 + rhoy * std::cos((pi * y) / (two * length)) + rhox * std::sin((pi * x) / length)) *
            (uvel0 + uvely * std::cos((three * pi * y) / (five * length)) +
            uvelx * std::sin((three * pi * x) / (two * length)))) / length +

            (two * pi * vvely * std::cos((two * pi * y) / (three * length)) *
            (rho0 + rhoy * std::cos((pi * y) / (two * length)) + rhox * std::sin((pi * x) / length)) *
            (uvel0 + uvely * std::cos((three * pi * y) / (five * length)) +
            uvelx * std::sin((three * pi * x) / (two * length)))) / (three * length) +

            (pi * rhox * std::cos((pi * x) / length) *
            std::pow((uvel0 + uvely * std::cos((three * pi * y) / (five * length)) +
            uvelx * std::sin((three * pi * x) / (two * length))), 2)) / length -

            (two * pi * pressx * std::sin((two * pi * x) / length)) / length -

            (pi * rhoy * (uvel0 + uvely * std::cos((three * pi * y) / (five * length)) +
            uvelx * std::sin((three * pi * x) / (two * length))) *
            std::sin((pi * y) / (two * length)) *
            (vvel0 + vvelx * std::cos((pi * x) / (two * length)) +
            vvely * std::sin((two * pi * y) / (three * length)))) / (two * length) -

            (three * pi * uvely * (rho0 + rhoy * std::cos((pi * y) / (two * length)) +
            rhox * std::sin((pi * x) / length)) *
            std::sin((three * pi * y) / (five * length)) *
            (vvel0 + vvelx * std::cos((pi * x) / (two * length)) +
            vvely * std::sin((two * pi * y) / (three * length)))) / (five * length);*/

    return result;
}

    double ymtmconv(double L, double x, double y)
    {
        
        
        double ymtmconv =     (Pi*pressy*cos((Pi*y)/L))/L -
        (Pi*vvelx*sin((Pi*x)/(2.*L))*
         (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
          rhox*sin((Pi*x)/L))*
         (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
          uvelx*sin((3*Pi*x)/(2.*L))))/(2.*L) +
        (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*
         (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
          rhox*sin((Pi*x)/L))*
         (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
          vvely*sin((2*Pi*y)/(3.*L))))/(2.*L) +
        (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*
         (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
          rhox*sin((Pi*x)/L))*
         (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
          vvely*sin((2*Pi*y)/(3.*L))))/(3.*L) +
        (Pi*rhox*cos((Pi*x)/L)*
         (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
          uvelx*sin((3*Pi*x)/(2.*L)))*
         (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
          vvely*sin((2*Pi*y)/(3.*L))))/L -
        (Pi*rhoy*sin((Pi*y)/(2.*L))*
         pow((vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
        vvely*sin((2*Pi*y)/(3.*L))),2))/(2.*L);
        
        
        
        
        
        /*(pi*pressy*cos((pi*y)/length))/length -
        (pi*vvelx*sin((pi*x)/(two*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +
        rhox*sin((pi*x)/length))*(uvel0 + uvely*cos((three*pi*y)/(five*length)) +
        uvelx*sin((three*pi*x)/(two*length))))/(two*length) +
        (three*pi*uvelx*cos((three*pi*x)/(two*length)) *
        (rho0 + rhoy*cos((pi*y)/(two*length)) + rhox*sin((pi*x)/length)) *
        (vvel0 + vvelx*cos((pi*x)/(two*length)) +
        vvely*sin((two*pi*y)/(three*length))))/(two*length) +
        (four*pi*vvely*cos((two*pi*y) /
        (three*length))*(rho0 + rhoy*cos((pi*y)/(two*length)) +
        rhox*sin((pi*x)/length))*(vvel0 + vvelx*cos((pi*x)/(two*length)) +
        vvely*sin((two*pi*y)/(three*length))))/(three*length) +
        (pi*rhox*cos((pi*x)/length) *
        (uvel0 + uvely*cos((three*pi*y)/(five*length)) +
        uvelx*sin((three*pi*x)/(two*length))) *
        (vvel0 + vvelx*cos((pi*x)/(two*length)) +
        vvely*sin((two*pi*y)/(three*length))))/length -
        (pi*rhoy*sin((pi*y)/(two*length)) *
        pow((vvel0 + vvelx*cos((pi*x)/(two*length)) +
        vvely*sin((two*pi*y)/(three*length))), 2))/(two*length);*/

    
        return  ymtmconv;
 }
    
    double energyconv(double L, double x, double y)
    {
        double energyconv=(uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                              uvelx*sin((3*Pi*x)/(2.*L)))*
                             ((-2*Pi*pressx*sin((2*Pi*x)/L))/L +
                               (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                  rhox*sin((Pi*x)/L))*
                                ((-2*Pi*pressx*sin((2*Pi*x)/L))/
                                   ((-1 + gamma)*L*
                                     (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                        rhox*sin((Pi*x)/L))) +
                                   ((3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*
                                        (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                          uvelx*sin((3*Pi*x)/(2.*L))))/L -
                                     (Pi*vvelx*sin((Pi*x)/(2.*L))*
                                        (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                           vvely*sin((2*Pi*y)/(3.*L))))/L)/2. -
                                  (Pi*rhox*cos((Pi*x)/L)*
                                     (press0 + pressx*cos((2*Pi*x)/L) +
                                       pressy*sin((Pi*y)/L)))/
                                    ((-1 + gamma)*L*
                                     pow((rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                        rhox*sin((Pi*x)/L)),2))) +
                                (Pi*rhox*cos((Pi*x)/L)*
                                  ((pow(wvel0,2) +
                                       pow((uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                           uvelx*sin((3*Pi*x)/(2.*L))),2) +
                                        pow((vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                          vvely*sin((2*Pi*y)/(3.*L))),2))/2. +
                                     (press0 + pressx*cos((2*Pi*x)/L) +
                                       pressy*sin((Pi*y)/L))/
                                      ((-1 + gamma)*
                                        (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                          rhox*sin((Pi*x)/L)))))/L) +
                             (3*Pi*uvelx*cos((3*Pi*x)/(2.*L))*
                                (press0 + pressx*cos((2*Pi*x)/L) +
                                  pressy*sin((Pi*y)/L) +
                                  (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                     rhox*sin((Pi*x)/L))*
                                  ((pow(wvel0,2) +
                                       pow((uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                           uvelx*sin((3*Pi*x)/(2.*L))),2) +
                                        pow((vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                          vvely*sin((2*Pi*y)/(3.*L))),2))/2. +
                                     (press0 + pressx*cos((2*Pi*x)/L) +
                                        pressy*sin((Pi*y)/L))/
                                    ((-1 + gamma)*
                                      (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                        rhox*sin((Pi*x)/L))))))/(2.*L) +
                           (2*Pi*vvely*cos((2*Pi*y)/(3.*L))*
                              (press0 + pressx*cos((2*Pi*x)/L) +
                                pressy*sin((Pi*y)/L) +
                                (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                   rhox*sin((Pi*x)/L))*
                                 ((pow(wvel0,2) +
                                      pow((uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                         uvelx*sin((3*Pi*x)/(2.*L))),2) +
                                      pow((vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                         vvely*sin((2*Pi*y)/(3.*L))),2))/2. +
                                   (press0 + pressx*cos((2*Pi*x)/L) +
                                      pressy*sin((Pi*y)/L))/
                                    ((-1 + gamma)*
                                      (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                        rhox*sin((Pi*x)/L))))))/(3.*L) +
                           (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                              vvely*sin((2*Pi*y)/(3.*L)))*
                            ((Pi*pressy*cos((Pi*y)/L))/L -
                              (Pi*rhoy*sin((Pi*y)/(2.*L))*
                                 ((pow(wvel0,2) +
                                      pow((uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                         uvelx*sin((3*Pi*x)/(2.*L))),2) +
                                      pow((vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                         vvely*sin((2*Pi*y)/(3.*L))),2))/2. +
                                   (press0 + pressx*cos((2*Pi*x)/L) +
                                      pressy*sin((Pi*y)/L))/
                                    ((-1 + gamma)*
                                      (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                        rhox*sin((Pi*x)/L)))))/(2.*L) +
                              (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                 rhox*sin((Pi*x)/L))*
                               ((Pi*pressy*cos((Pi*y)/L))/
                                  ((-1 + gamma)*L*
                                    (rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                      rhox*sin((Pi*x)/L))) +
                                 ((-6*Pi*uvely*
                                       (uvel0 + uvely*cos((3*Pi*y)/(5.*L)) +
                                         uvelx*sin((3*Pi*x)/(2.*L)))*
                                       sin((3*Pi*y)/(5.*L)))/(5.*L) +
                                    (4*Pi*vvely*cos((2*Pi*y)/(3.*L))*
                                       (vvel0 + vvelx*cos((Pi*x)/(2.*L)) +
                                         vvely*sin((2*Pi*y)/(3.*L))))/(3.*L))/
                                  2. + (Pi*rhoy*sin((Pi*y)/(2.*L))*
                                    (press0 + pressx*cos((2*Pi*x)/L) +
                                      pressy*sin((Pi*y)/L)))/
                                  (2.*(-1 + gamma)*L*
                                    pow((rho0 + rhoy*cos((Pi*y)/(2.*L)) +
                                     rhox*sin((Pi*x)/L)),2))));
       
        
        return  energyconv;
    }
    
    
    
};









#endif /* Dv_nozzle_hpp */
