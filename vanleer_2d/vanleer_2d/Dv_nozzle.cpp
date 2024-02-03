//
//  Dv_nozzle.cpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#include "Dv_nozzle.hpp"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;




    void Dv_nozzle::initialization (double P_0, double T_0)
{
        
        
        
        _P_0= P_0;
        _T_0= T_0;
        int k = 4;
        int x = _imax_f+2;
        int y = _imax_f+2;
        
        
       
        _F= new double**[4];
        epsi_plus_y= new double**[4];
        epsi_plus_x= new double**[4];
        epsi_minus_y= new double**[4];
        epsi_minus_x= new double**[4];
        _V= new double**[4];
        _U= new double**[4];
        _F_eta= new double**[4];
       
        
        
        
        for (int i = 0; i <k; i++) {
            
            epsi_plus_y[i] = new double*[x];
            epsi_plus_x[i] = new double*[x];
            epsi_minus_y[i] = new double*[x];
            epsi_minus_x[i] = new double*[x];
            _V[i] = new double*[x];
            _U[i] = new double*[x];
            _F[i] = new double*[x];
            _F_eta[i] = new double*[x];
            
            
                   for (int j = 0; j < y; j++) {
                    
                      
                       epsi_plus_y[i][j] = new double[y];
                       epsi_plus_x[i][j] = new double[y];
                       epsi_minus_y[i][j] = new double[y];
                       epsi_minus_x[i][j] = new double[y];
                       _V[i][j] = new double[y];
                       _U[i][j] = new double[y];
                       _F[i][j] = new double[y];
                       _F_eta[i][j] = new double[y];
                   }
               }
        
        
        
     
        
        
        
        
        
        
        
        
        
        
        
       
       
        
        
        // 0 is rho, 1 is velocity,  2 is pressure
        //values calculated at cell as cell average
        for (int j= _imin; j<=im_c; j++)
        {
            for (int i= _imin; i<= im_c; i++)
            {
                
                
                
               
                
                
                _V[1][i][j] = uvel0;
                _V[2][i][j] = vvel0;
                
                _V[3][i][j]= press0;
                
                _V[0][i][j]= rho0;
                
                
             
                
                T[i][j]=_V[3][i][j]/(R_air*_V[0][i][j]);
                
                et[i][j] = (R_air/(gamma-1))*T[i][j] + 0.5* (_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j]);
                
                _U[0][i][j]= _V[0][i][j];
                
                _U[1][i][j]= _V[0][i][j]*_V[1][i][j];
                
                _U[2][i][j]= _V[0][i][j]*_V[2][i][j];
                
                _U[3][i][j]=_V[0][i][j]*et[i][j];
                
                
            
                
                
            }
         
            
        }
      
        
        //calculate epsi- kasi direction kasi corresponds to x kinda
        
        for (int j=0;j<=im_c;j++)
        {
            for (int i=2;i<=im_f-2;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i-1][j]),1e-6),_V[0][i][j]-_V[0][i-1][j]);
                r_plus= (_V[0][i+1][j]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i-1][j]-_V[0][i-2][j])/(denm);
                epsi_plus_x[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i-1][j]),1e-6),_V[1][i][j]-_V[1][i-1][j]);
                r_plus= (_V[1][i+1][j]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i-1][j]-_V[1][i-2][j])/(denm);
                epsi_plus_x[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i-1][j]),1e-6),_V[2][i][j]-_V[2][i-1][j]);
                r_plus= (_V[2][i+1][j]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i-1][j]-_V[2][i-2][j])/(denm);
                epsi_plus_x[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i-1][j]),1e-6),_V[3][i][j]-_V[3][i-1][j]);
                r_plus= (_V[3][i+1][j]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i-1][j]-_V[3][i-2][j])/(denm);
                epsi_plus_x[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
        
        
        //calculate epsi-eta direction  corresponds to y kinda
        
        for (int j=2;j<=im_f-2;j++)
        {
            for (int i=0;i<=im_c;i++)
            {
                
              
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i][j-1]),1e-6),_V[0][i][j]-_V[0][i][j-1]);
                r_plus= (_V[0][i][j+1]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i][j-1]-_V[0][i][j-2])/(denm);
                epsi_plus_y[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
               
                
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i][j-1]),1e-6),_V[1][i][j]-_V[1][i][j-1]);
                r_plus= (_V[1][i][j+1]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i][j-1]-_V[1][i][j-2])/(denm);
                epsi_plus_y[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i][j-1]),1e-6),_V[2][i][j]-_V[2][i][j-1]);
                r_plus= (_V[2][i][j+1]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i][j-1]-_V[2][i][j-2])/(denm);
                epsi_plus_y[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i][j-1]),1e-6),_V[3][i][j]-_V[3][i][j-1]);
                r_plus= (_V[3][i][j+1]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i][j-1]-_V[3][i][j-2])/(denm);
                epsi_plus_y[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                
            }
            
        }
        
        vanleer();
        
        
        
        
        
        
    }



void Dv_nozzle::set_boundary_cond()
{
   
    for (int j=0;j<=im_c;j++)
        
    {
        //inflow boundary conditions
       // mach_in = 0.5*(3.0*_mach[0][j] - _mach[1][j] );
       // epsi_1= 1 + ((gamma-1)*mach_in*mach_in)/2;
       // T_in = _T_0/epsi_1;
        
     //   Ghost cells at inlet


      //  ----- first ghost
        Delta_ghost_x_y = 0.5*(y_n[0][j+1] -y_n[0][j]);

        y_c_ghost= y_n[0][j]+ Delta_ghost_x_y;

        x_c_ghost=x_n[0][j]- Delta_ghost_x_y;

        _V_ghost_inflow[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_inflow[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_inflow[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_inflow[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T
        
       // ------second ghost
        Delta_ghost_x_y = 0.5*(y_n[0][j+1]-y_n[0][j]);

        y_c_ghost= y_n[0][j]+ Delta_ghost_x_y;

        x_c_ghost= x_c_ghost - 2*Delta_ghost_x_y;
       
        _V_ghost_inflow_2[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_inflow_2[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_inflow_2[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_inflow_2[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_2[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

        
       
        
        
        
      //  -----third ghost
        Delta_ghost_x_y = 0.5*(y_n[0][j+1]- y_n[0][j]);

        y_c_ghost= y_n[0][j]+ Delta_ghost_x_y;

        x_c_ghost= x_c_ghost - 2*Delta_ghost_x_y;

        _V_ghost_inflow_3[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_inflow_3[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_inflow_3[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_inflow_3[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_3[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

   

        
        
        
        
        /*
        u_in =  uvel0;
        v_in= vvel0;
        P_in= press0;
        
        rho_in= rho0;
        T_in=P_in/(R_air*rho_in);
        
        */
        
        
        // et_in = (R_air/(gamma-1))*T_in + 0.5*u_in*u_in;
        
        //extrapolate U and V for 1st ghost cell
        /*
         _U_ghost_inflow[0]=2*_U[0][_imin]-_U[0][_imin-1];
         _U_ghost_inflow[1]=2*_U[1][_imin]-_U[1][_imin-1];
         _U_ghost_inflow[2]=2*_U[2][_imin]-_U[2][_imin-1];
         */
        
        /*
        //first ghoast cell
        
        
        
        
        _V_ghost_inflow[0][j]=2*rho_in-_V[0][_imin][j];
        _V_ghost_inflow[1][j]=2*u_in-_V[1][_imin][j];
        _V_ghost_inflow[2][j]=2*v_in-_V[2][_imin][j];
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_inflow[3][j]=2*P_in-_V[3][_imin][j];
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow[5][j]=2*T_in-T[_imin][j];    //T
        
        //second ghost cell
        _V_ghost_inflow_2[0][j]=2*_V_ghost_inflow[0][j]-_V[0][_imin][j];
        _V_ghost_inflow_2[1][j]=2*_V_ghost_inflow[1][j]-_V[1][_imin][j];
        
        _V_ghost_inflow_2[2][j]=2*_V_ghost_inflow[2][j]-_V[2][_imin][j];
        _V_ghost_inflow_2[3][j]=2*_V_ghost_inflow[3][j]-_V[3][_imin][j];
        //  _V_ghost_inflow_2[2][j]=0;
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_2[5][j]=2*_V_ghost_inflow[4][j]-T[_imin][j];
        
        //third ghost cell
        _V_ghost_inflow_3[0][j]=2*_V_ghost_inflow_2[0][j]-_V_ghost_inflow[0][j];
        _V_ghost_inflow_3[1][j]=2*_V_ghost_inflow_2[1][j]-_V_ghost_inflow[1][j];
        
        _V_ghost_inflow_3[2][j]=2*_V_ghost_inflow_2[2][j]-_V_ghost_inflow[2][j];
        _V_ghost_inflow_3[3][j]=2*_V_ghost_inflow_2[3][j]-_V_ghost_inflow[3][j];
        // _V_ghost_inflow_3[2][j]=0;
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_inflow_3[5][j]=2*_V_ghost_inflow_2[4][j]-_V_ghost_inflow[4][j];
        
        
        */
        
        
        
        
        
        
        //calculate epsi at 1
        
        
        
        //for density
        double denm= copysign(max(abs(_V[0][1][j]-_V[0][0][j]),1e-6),_V[0][1][j]-_V[0][0][j]);
        r_plus= (_V[0][2][j]-_V[0][1][j])/(denm);
        r_minus= (_V[0][0][j]-_V_ghost_inflow[0][j])/(denm);
        epsi_plus_x[0][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[0][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u- velocity
        denm= copysign(max(abs(_V[1][1][j]-_V[1][0][j]),1e-6),_V[1][1][j]-_V[1][0][j]);
        r_plus= (_V[1][2][j]-_V[1][1][j])/(denm);
        r_minus= (_V[1][0][j]-_V_ghost_inflow[1][j])/(denm);
        epsi_plus_x[1][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[1][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v velocity
        denm= copysign(max(abs(_V[2][1][j]-_V[2][0][j]),1e-6),_V[2][1][j]-_V[2][0][j]);
        r_plus= (_V[2][2][j]-_V[2][1][j])/(denm);
        r_minus= (_V[2][0][j]-_V_ghost_inflow[2][j])/(denm);
        epsi_plus_x[2][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[2][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][1][j]-_V[3][0][j]),1e-6),_V[3][1][j]-_V[3][0][j]);
        r_plus= (_V[3][2][j]-_V[3][1][j])/(denm);
        r_minus= (_V[3][0][j]-_V_ghost_inflow[3][j])/(denm);
        epsi_plus_x[3][1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[3][1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        //calculate epsi at 0
        
        //for density
        denm= copysign(max(abs(_V[0][0][j]-_V_ghost_inflow[0][j]),1e-6),_V[0][0][j]-_V_ghost_inflow[0][j]);
        r_plus= (_V[0][1][j]-_V[0][0][j])/(denm);
        r_minus= (_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j])/(denm);
        epsi_plus_x[0][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[0][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u velocity
        denm= copysign(max(abs(_V[1][0][j]-_V_ghost_inflow[1][j]),1e-6),_V[1][0][j]-_V_ghost_inflow[1][j]);
        r_plus= (_V[1][1][j]-_V[1][0][j])/(denm);
        r_minus= (_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j])/(denm);
        epsi_plus_x[1][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[1][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v- velocity
        denm= copysign(max(abs(_V[2][0][j]-_V_ghost_inflow[2][j]),1e-6),_V[2][0][j]-_V_ghost_inflow[2][j]);
        r_plus= (_V[2][1][j]-_V[2][0][j])/(denm);
        r_minus= (_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j])/(denm);
        epsi_plus_x[2][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[2][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][0][j]-_V_ghost_inflow[3][j]),1e-6),_V[3][0][j]-_V_ghost_inflow[3][j]);
        r_plus= (_V[3][1][j]-_V[3][0][j])/(denm);
        r_minus= (_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j])/(denm);
        epsi_plus_x[3][0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x[3][0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        //calculate epsi at imin-1
        
        //for density
        denm= copysign(max(abs(_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]),1e-6),_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]);
        r_plus= (_V[0][0][j]-_V_ghost_inflow[0][j])/(denm);
        r_minus= (_V_ghost_inflow_2[0][j]-_V_ghost_inflow_3[0][j])/(denm);
        epsi_plus_x_ghost[0][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[0][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u-velocity
        denm= copysign(max(abs(_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]),1e-6),_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]);
        r_plus= (_V[1][0][j]-_V_ghost_inflow[1][j])/(denm);
        r_minus= (_V_ghost_inflow_2[1][j]-_V_ghost_inflow_3[1][j])/(denm);
        epsi_plus_x_ghost[1][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[1][j]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v -velocity
        
        denm= copysign(max(abs(_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]),1e-6),_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]);
        r_plus= (_V[2][0][j]-_V_ghost_inflow[2][j])/(denm);
        r_minus= (_V_ghost_inflow_2[2][j]-_V_ghost_inflow_3[2][j])/(denm);
        epsi_plus_x_ghost[2][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[2][j]= (r_minus+abs(r_minus))/(1+r_minus);
        
        //for pressure
        denm= copysign(max(abs(_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]),1e-6),_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]);
        r_plus= (_V[3][0][j]-_V_ghost_inflow[3][j])/(denm);
        r_minus= (_V_ghost_inflow_2[3][j]-_V_ghost_inflow_3[3][j])/(denm);
        epsi_plus_x_ghost[3][j]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_x_ghost[3][j]= (r_minus+abs(r_minus))/(1+r_minus);
        /*
        //flux at face 0
        
        
        
        double a2=face_area(x_n[0][j+1], x_n[0][j],y_n[0][j+1], y_n[0][j]);
        nx=n_x(y_n[0][j+1], y_n[0][j], a2);
        ny=n_y(x_n[0][j+1], x_n[0][j], a2);
        T_vel_roe=abs((nx*u_in))+abs((ny*v_in));
        
        double ht_in=(gamma/(gamma-1))*(P_in/rho_in) + (u_in*u_in+v_in*v_in)*0.5;
        
        
        _F[0][0][j]=rho_in*T_vel_roe ;
        _F[1][0][j]= rho_in*u_in*T_vel_roe + P_in*nx;
        _F[2][0][j]= rho_in*v_in*T_vel_roe + P_in*ny;
        _F[3][0][j]= rho_in*ht_in*T_vel_roe;
        
        */
        /*
        if((j>1)&&(j<im_f-1))
        {
            
            roe_boundary_F_eta(0, j);
            
        }
        */
        
        
        
        
        //flux at face 1
        
        
        
        //density
        V_L[0]= _V[0][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j]) + (1+kappa)*epsi_minus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j])        );
        
        V_R[0]= _V[0][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][2][j]*(_V[0][2][j]-_V[0][1][j]) + (1+kappa)*epsi_plus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j])        );
        //velocity
        V_L[1]= _V[1][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j]) + (1+kappa)*epsi_minus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j])        );
        
        V_R[1]= _V[1][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][2][j]*(_V[1][2][j]-_V[1][1][j]) + (1+kappa)*epsi_plus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j])        );
        //pressure
        V_L[2]= _V[2][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j]) + (1+kappa)*epsi_minus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j])        );
        
        V_R[2]= _V[2][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][2][j]*(_V[2][2][j]-_V[2][1][j]) + (1+kappa)*epsi_plus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j])        );
        //pressure
        V_L[3]= _V[3][0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j]) + (1+kappa)*epsi_minus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j])        );
        
        V_R[3]= _V[3][1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][2][j]*(_V[3][2][j]-_V[3][1][j]) + (1+kappa)*epsi_plus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j])        );
        
        ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + (V_L[1]*V_L[1]+V_L[2]*V_L[2])*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + (V_R[1]*V_R[1]+V_R[2]*V_R[2])*0.5;
        
        
        
        
        
        vanleer_kasi(1,j);
        
       
        
        
        
  
        
        // flux at 0
        
        
        //density
        V_L[0]= _V_ghost_inflow[0][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[0][j]*(_V_ghost_inflow[0][j]-_V_ghost_inflow_2[0][j]) + (1+kappa)*epsi_minus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j]));
        
        V_R[0]= _V[0][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][1][j]*(_V[0][1][j]-_V[0][0][j]) + (1+kappa)*epsi_plus_x[0][0][j]*(_V[0][0][j]-_V_ghost_inflow[0][j])        );
        
        
        
        //u-velocity
        V_L[1]= _V_ghost_inflow[1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[1][j]*(_V_ghost_inflow[1][j]-_V_ghost_inflow_2[1][j]) + (1+kappa)*epsi_minus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j]));
        
        V_R[1]= _V[1][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][1][j]*(_V[1][1][j]-_V[1][0][j]) + (1+kappa)*epsi_plus_x[1][0][j]*(_V[1][0][j]-_V_ghost_inflow[1][j])        );
        //v-velocity
        V_L[2]= _V_ghost_inflow[2][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[2][j]*(_V_ghost_inflow[2][j]-_V_ghost_inflow_2[2][j]) + (1+kappa)*epsi_minus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j]));
        
        V_R[2]= _V[2][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][1][j]*(_V[2][1][j]-_V[2][0][j]) + (1+kappa)*epsi_plus_x[2][0][j]*(_V[2][0][j]-_V_ghost_inflow[2][j])        );
        //pressure
        V_L[3]= _V_ghost_inflow[3][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x_ghost[3][j]*(_V_ghost_inflow[3][j]-_V_ghost_inflow_2[3][j]) + (1+kappa)*epsi_minus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j]));
        
        V_R[3]= _V[3][0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][1][j]*(_V[3][1][j]-_V[3][0][j]) + (1+kappa)*epsi_plus_x[3][0][j]*(_V[3][0][j]-_V_ghost_inflow[3][j])        );
        
        //total enthalpy
        ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + (V_L[1]*V_L[1]+V_L[2]*V_L[2])*0.5;
        ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + (V_R[1]*V_R[1]+V_R[2]*V_R[2])*0.5;
        
        
        
        
        vanleer_kasi(0,j);

       
  /*
        double a2=face_area(x_n[0][j+1], x_n[0][j],y_n[0][j+1], y_n[0][j]);
        nx=n_x(y_n[0][j+1], y_n[0][j], a2);
        ny=n_y(x_n[0][j+1], x_n[0][j], a2);
        T_vel_roe=abs((nx*u_in))+abs((ny*v_in));
        
        double ht_in=(gamma/(gamma-1))*(P_in/rho_in) + (u_in*u_in+v_in*v_in)*0.5;
        
        
        _F[0][0][j]=rho_in*T_vel_roe ;
        _F[1][0][j]= rho_in*u_in*T_vel_roe + P_in*nx;
        _F[2][0][j]= rho_in*v_in*T_vel_roe + P_in*ny;
        _F[3][0][j]= rho_in*ht_in*T_vel_roe;
        
 */
       
        /*
        _F[0][1][j]= 2*_F[0][2][j]-_F[0][3][j];
        _F[1][1][j]= 2*_F[1][2][j]-_F[1][3][j];
        _F[2][1][j]= 2*_F[2][2][j]-_F[2][3][j];
        _F[3][1][j]= 2*_F[3][2][j]-_F[3][3][j];
        */
        
        
       
        
        //******************************************
        // outflow boundary conditions
        
        
           
            // P_out = 125000;
            P_out= 0.5*(3*_V[3][im_c][j]-_V[3][im_c-1][j]);
            // std::cout<<"back pressure is "<<P_out<<std::endl;
            
            
            rho_out= 0.5*(3*_U[0][im_c][j]-_U[0][im_c-1][j]);
            
            
            
           u_out = (0.5*(3*_U[1][im_c][j]-_U[1][im_c-1][j]))/rho_out;
          v_out = (0.5*(3*_U[2][im_c][j]-_U[2][im_c-1][j]))/rho_out;
       // u_out = (0.5*(3*_V[1][im_c][j]-_V[1][im_c-1][j]));
    //    v_out = (0.5*(3*_V[2][im_c][j]-_V[2][im_c-1][j]));
        T_out=P_out/(rho_out*R_air);
          
            
            
           // ///here
            
        //Ghost cells at outlet


     //   ----- first ghost
        Delta_ghost_x_y = 0.5*(y_n[i_n_max-1][j+1] - y_n[i_n_max-1][j]);

        y_c_ghost= y_n[i_n_max-1][j]+ Delta_ghost_x_y;

        x_c_ghost=x_n[i_n_max-1][j]+ Delta_ghost_x_y;

       
        _V_ghost_outflow[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_outflow[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_outflow[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_outflow[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_outflow[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T
        
        //cout<<" density at lastt cell "<<im_c<<","<<j<<" "<<rho(1,x_c[im_c][j],y_c[im_c][j])<<endl;

       // cout<<" last cell center is "<<x_c[im_c][j]<<","<<y_c[im_c][j]<<" "<<" at  cell "<<im_c<<","<<j<<endl;
       // cout<<" press at first ghost cell at jmax cell "<<j<<"  "<<press(1,x_c_ghost,y_c_ghost)<<endl;

      //  cout<<" last cell center is "<<x_c_ghost<<","<<y_c_ghost<<" "<<" at  cell along i  "<<j<<endl;

     //   ------second ghost
        Delta_ghost_x_y = 0.5*(y_n[i_n_max-1][j+1] - y_n[i_n_max-1][j]);

        y_c_ghost= y_n[i_n_max-1][j]+ Delta_ghost_x_y;

        x_c_ghost= x_c_ghost + 2*Delta_ghost_x_y;

        
        _V_ghost_outflow_2[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_outflow_2[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_outflow_2[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_outflow_2[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_outflow[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

   //     -----third ghost
        Delta_ghost_x_y = 0.5*(y_n[i_n_max-1][j+1] -y_n[i_n_max-1][j]);

        y_c_ghost= y_n[i_n_max-1][j]+ Delta_ghost_x_y;

        x_c_ghost= x_c_ghost + 2*Delta_ghost_x_y;
        _V_ghost_outflow_3[0][j]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_outflow_3[1][j]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_outflow_3[2][j]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_outflow_3[3][j]=press(1,x_c_ghost,y_c_ghost);
        
        //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
        _V_ghost_outflow_3[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

            
            /*
            //ghost cells
            //1st ghost
            //
            //
            _V_ghost_outflow[0][j]=2*_V[0][im_c][j]-_V[0][im_c-1][j];
            _V_ghost_outflow[1][j]=2*_V[1][im_c][j]-_V[1][im_c-1][j];
            
            _V_ghost_outflow[2][j]=2*_V[2][im_c][j]-_V[2][im_c-1][j];
            
            _V_ghost_outflow[3][j]=2*P_out-_V[3][im_c][j];
            
            //_V_ghost_outflow[3]=2*et[_imin]-et[_imin];  //et
            _V_ghost_outflow[5][j]=2*T[im_c][j]-T[im_c-1][j];   //T
            
            ///second ghost cell
            ///
            ///
            _V_ghost_outflow_2[0][j]=2*_V_ghost_outflow[0][j]- _V[0][im_c][j];
            _V_ghost_outflow_2[1][j]=2*_V_ghost_outflow[1][j]- _V[1][im_c][j];
            _V_ghost_outflow_2[2][j]=2*_V_ghost_outflow[2][j]- _V[2][im_c][j];
            _V_ghost_outflow_2[3][j]=2*_V_ghost_outflow[3][j]-_V[3][im_c][j];
            
            
            //  _V_ghost_outflow_2[3]=2*_V_ghost_outflow[3]-et[im_c];
            _V_ghost_outflow_2[5][j]=2*_V_ghost_outflow[5][j]-T[im_c][j];
            
            //third ghost cell
            ///
            ///
            _V_ghost_outflow_3[0][j]=2*_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j];
            _V_ghost_outflow_3[1][j]=2*_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j];
            _V_ghost_outflow_3[2][j]=2*_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j];
            _V_ghost_outflow_3[3][j]=2*_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j];
            //  _V_ghost_outflow_2[3]=2*_V_ghost_outflow[3]-et[im_c];
            _V_ghost_outflow_3[5][j]=2*_V_ghost_outflow_2[5][j]-_V_ghost_outflow[5][j];
            
            */
            
        
        
        
        
        
        
        
        
        
            
            //calculate epsi at imf-1 or imc
            
            
            //for density
            denm= copysign(max(abs(_V[0][im_c][j]-_V[0][im_c-1][j]),1e-6),_V[0][im_c][j]-_V[0][im_c-1][j]);
            r_plus= (_V_ghost_outflow[0][j]-_V[0][im_c][j])/(denm);
            r_minus= (_V[0][im_c-1][j]-_V[0][im_c-2][j])/(denm);
            epsi_plus_x[0][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for u-velocity
            denm= copysign(max(abs(_V[1][im_c][j]-_V[1][im_c-1][j]),1e-6),_V[1][im_c][j]-_V[1][im_c-1][j]);
            r_plus= (_V_ghost_outflow[1][j]-_V[1][im_c][j])/(denm);
            r_minus= (_V[1][im_c-1][j]-_V[1][im_c-2][j])/(denm);
            epsi_plus_x[1][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for v-velocity
            denm= copysign(max(abs(_V[2][im_c][j]-_V[2][im_c-1][j]),1e-6),_V[2][im_c][j]-_V[2][im_c-1][j]);
            r_plus= (_V_ghost_outflow[2][j]-_V[2][im_c][j])/(denm);
            r_minus= (_V[2][im_c-1][j]-_V[2][im_c-2][j])/(denm);
            epsi_plus_x[2][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm= copysign(max(abs(_V[3][im_c][j]-_V[3][im_c-1][j]),1e-6),_V[3][im_c][j]-_V[3][im_c-1][j]);
            r_plus= (_V_ghost_outflow[3][j]-_V[3][im_c][j])/(denm);
            r_minus= (_V[3][im_c-1][j]-_V[3][im_c-2][j])/(denm);
            epsi_plus_x[3][im_c][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[3][im_c][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            
            
            
            //calculate epsi for imf or imc+1
            
            
            //for density
            denm= copysign(max(abs(_V_ghost_outflow[0][j]-_V[0][im_c][j]),1e-6),_V_ghost_outflow[0][j]-_V[0][im_c][j]);
            r_plus= (_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j])/(denm);
            r_minus= (_V[0][im_c][j]-_V[0][im_c-1][j])/(denm);
            epsi_plus_x[0][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for u-velocity
            denm= copysign(max(abs(_V_ghost_outflow[1][j]-_V[1][im_c][j]),1e-6),_V_ghost_outflow[1][j]-_V[1][im_c][j]);
            r_plus= (_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j])/(denm);
            r_minus= (_V[1][im_c]-_V[1][im_c-1])/(denm);
            epsi_plus_x[1][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for v-velocity
            denm= copysign(max(abs(_V_ghost_outflow[2][j]-_V[2][im_c][j]),1e-6),_V_ghost_outflow[2][j]-_V[2][im_c][j]);
            r_plus= (_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j])/(denm);
            r_minus= (_V[2][im_c][j]-_V[2][im_c-1][j])/(denm);
            epsi_plus_x[2][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            //for pressure
            denm= copysign(max(abs(_V_ghost_outflow[3][j]-_V[3][im_c][j]),1e-6),_V_ghost_outflow[3][j]-_V[3][im_c][j]);
            r_plus= (_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j])/(denm);
            r_minus= (_V[3][im_c][j]-_V[3][im_c-1][j])/(denm);
            epsi_plus_x[3][im_c+1][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[3][im_c+1][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            //calculate epsi for imf+1 or imc+2
            
            
            //for density
            denm=copysign(max(abs(_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]),1e-6),_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]);
            
            r_plus= (_V_ghost_outflow_3[0][j]-_V_ghost_outflow_2[0][j])/(denm);
            r_minus= (_V_ghost_outflow[0][j]-_V[0][im_c][j])/(denm);
            
            epsi_plus_x[0][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[0][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for u-velocity
            denm
            =copysign(max(abs(_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]),1e-6),_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]);
            
            
            r_plus= (_V_ghost_outflow_3[1][j]-_V_ghost_outflow_2[1][j])/(denm);
            r_minus= (_V_ghost_outflow[1][j]-_V[1][im_c][j])/(denm);
            
            epsi_plus_x[1][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[1][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for v-velocity
            denm
            =copysign(max(abs(_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]),1e-6),_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]);
            
            
            r_plus= (_V_ghost_outflow_3[2][j]-_V_ghost_outflow_2[2][j])/(denm);
            r_minus= (_V_ghost_outflow[2][j]-_V[2][im_c][j])/(denm);
            
            epsi_plus_x[2][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[2][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm
            =copysign(max(abs(_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j]),1e-6),_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j]);
            
            
            r_plus= (_V_ghost_outflow_3[3][j]-_V_ghost_outflow_2[3][j])/(denm);
            r_minus= (_V_ghost_outflow[3][j]-_V[3][im_c][j])/(denm);
            
            epsi_plus_x[3][im_c+2][j]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_x[3][im_c+2][j]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            //calculate fluxes at imf-1 or imc
            
            
           
            
            //density
            V_L[0]= _V[0][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][im_c-1][j]*(_V[0][im_c-1][j]-_V[0][im_c-2][j]) + (1+kappa)*epsi_minus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j])        );
            
            V_R[0]= _V[0][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j]) + (1+kappa)*epsi_plus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j])        );
            //u-velocity
            V_L[1]= _V[1][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][im_c-1][j]*(_V[1][im_c-1][j]-_V[1][im_c-2][j]) + (1+kappa)*epsi_minus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j])        );
            
            V_R[1]=_V[1][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][im_c+1][j]*(_V_ghost_outflow[1][j]-_V[1][im_c][j]) + (1+kappa)*epsi_plus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j])        );
            //v-veloity
            V_L[2]= _V[2][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][im_c-1][j]*(_V[2][im_c-1][j]-_V[2][im_c-2][j]) + (1+kappa)*epsi_minus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j])        );
            
            V_R[2]=_V[2][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j]) + (1+kappa)*epsi_plus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j])        );
        //pressure
        V_L[3]= _V[3][im_c-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][im_c-1][j]*(_V[3][im_c-1][j]-_V[3][im_c-2][j]) + (1+kappa)*epsi_minus_x[3][im_c][j]*(_V[3][im_c][j]-_V[3][im_c-1][j])        );
        
        V_R[3]=_V[3][im_c][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][im_c+1][j]*(_V_ghost_outflow[3][j]-_V[3][im_c][j]) + (1+kappa)*epsi_plus_x[3][im_c][j]*(_V[3][im_c][j]-_V[3][im_c-1][j])        );
            
            ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + (V_L[1]*V_L[1]+V_L[2]*V_L[2])*0.5;
            ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + (V_R[1]*V_R[1]+V_R[2]*V_R[2])*0.5;
            
            vanleer_kasi(im_f-1,j);
            
      
            
            
        
            
            
            // flux at im_f
            
            
            //density
            V_L[0]= _V[0][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][im_c][j]*(_V[0][im_c][j]-_V[0][im_c-1][j]) + (1+kappa)*epsi_minus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j])        );
            
            V_R[0]= _V_ghost_outflow[0][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][im_c+2][j]*(_V_ghost_outflow_2[0][j]-_V_ghost_outflow[0][j]) + (1+kappa)*epsi_plus_x[0][im_c+1][j]*(_V_ghost_outflow[0][j]-_V[0][im_c][j])        );
            
            
            //u-velocity
            V_L[1]= _V[1][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][im_c][j]*(_V[1][im_c][j]-_V[1][im_c-1][j]) + (1+kappa)*epsi_minus_x[1][im_c+1][j]*(_V_ghost_outflow[1][j]-_V[1][im_c][j])        );
            
            V_R[1]= _V_ghost_outflow[1][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][im_c+2][j]*(_V_ghost_outflow_2[1][j]-_V_ghost_outflow[1][j]) + (1+kappa)*epsi_plus_x[1][im_c+1][j]*(_V_ghost_outflow[1][j]-_V[1][im_c][j])        );
            //v-velocity
            V_L[2]= _V[2][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][im_c][j]*(_V[2][im_c][j]-_V[2][im_c-1][j]) + (1+kappa)*epsi_minus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j])        );
            
            V_R[2]= _V_ghost_outflow[2][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][im_c+2][j]*(_V_ghost_outflow_2[2][j]-_V_ghost_outflow[2][j]) + (1+kappa)*epsi_plus_x[2][im_c+1][j]*(_V_ghost_outflow[2][j]-_V[2][im_c][j])        );
        //pressure
        V_L[3]= _V[3][im_c][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][im_c][j]*(_V[3][im_c][j]-_V[3][im_c-1][j]) + (1+kappa)*epsi_minus_x[2][im_c+1][j]*(_V_ghost_outflow[3][j]-_V[3][im_c][j])        );
        
        V_R[3]= _V_ghost_outflow[3][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][im_c+2][j]*(_V_ghost_outflow_2[3][j]-_V_ghost_outflow[3][j]) + (1+kappa)*epsi_plus_x[3][im_c+1][j]*(_V_ghost_outflow[3][j]-_V[3][im_c][j])        );
            
           ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + (V_L[1]*V_L[1]+V_L[2]*V_L[2])*0.5;
           ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + (V_R[1]*V_R[1]+V_R[2]*V_R[2])*0.5;
            
            
            vanleer_kasi(im_f,j);
            
       
        
        
        
     /*
        a2=face_area(x_n[im_f][j+1], x_n[im_f][j],y_n[im_f][j+1], y_n[im_f][j]);
        nx=n_x(y_n[im_f][j+1], y_n[im_f][j], a2);
        ny=n_y(x_n[im_f][j+1], x_n[im_f][j], a2);
        T_vel_roe=abs(nx*u_out)+abs(ny*v_out);
        
       double ht_out=(gamma/(gamma-1))*(P_out/rho_out) + (u_out*u_out+v_out*v_out)*0.5;
        
        
        _F[0][im_f][j]=rho_out*T_vel_roe ;
        _F[1][im_f][j]= rho_out*u_out*T_vel_roe + P_out*nx;
        _F[2][im_f][j]= rho_out*v_out*T_vel_roe + P_out*ny;
        _F[3][im_f][j]= rho_out*ht_out*T_vel_roe;
        
        
       
        
        _F[0][1][j]= 2*_F[0][2][j]-_F[0][3][j];
        _F[1][1][j]= 2*_F[1][2][j]-_F[1][3][j];
        _F[2][1][j]= 2*_F[2][2][j]-_F[2][3][j];
        _F[3][1][j]= 2*_F[3][2][j]-_F[3][3][j];
        
        
        
        */
        
    }

   
    for(int i=0;i<=im_c;i++)
        
        
    {
        //for for eta faces at cell j=0
        
        
        //zero flux implementation
        
       
       /*
        double P_wall =_V[3][i][0] ;//- 0.5*(_V[3][i][1]-_V[3][i][0]);
        
        
        double a4=face_area(x_n[i+1][0], x_n[i][0],y_n[i+1][0], y_n[i][0]);
        nx= n_x(y_n[i+1][0], y_n[i][0], a4);
        ny= n_y(x_n[i+1][0], x_n[i][0], a4);
        
        _F_eta[0][i][0]=0;
        
        _F_eta[1][i][0]=(nx*P_wall);
        
        _F_eta[2][i][0]=(ny*P_wall);
        
        _F_eta[3][i][0]=0;
        
        
        */
        
        /*
        //ghost cell implementation
        double a4=face_area(x_n[i+1][0], x_n[i][0],y_n[i+1][0], y_n[i][0]);
        nx= n_x(y_n[i+1][0], y_n[i][0], a4);
        ny= n_y(x_n[i+1][0], x_n[i][0], a4);
        
        A=_V[1][i][0]*ny+_V[2][i][0]*nx;
        
        B= -_V[1][i][0]*nx-_V[2][i][0]*ny;
        
        
        
        
        
        
        _V_ghost_wall_0[2][i]=(B*ny-A*nx     )/(  ny*ny-nx*nx );
        
        _V_ghost_wall_0[1][i]=(A-_V_ghost_wall_0[2][i]*nx    )/(  ny );
        
        _V_ghost_wall_0[5][i]=T[i][0];
        
        _V_ghost_wall_0[3][i]=2*_V[3][i][0]-_V[3][i][1];
       
        
        _V_ghost_wall_0[0][i]= _V_ghost_wall_0[3][i]/( _V_ghost_wall_0[5][i]*R_air);
        
        
        
        //density
         V_D[0]= _V_ghost_wall_0[0][i];
        
         V_U[0]= _V[0][i][0];
        //x-velocity
         V_D[1]= _V_ghost_wall_0[1][i];
        
         V_U[1]= _V[1][i][0];
        //y-velocity
         V_D[2]= _V_ghost_wall_0[2][i];
        
         V_U[2]= _V[2][i][0];
        //pressure
         V_D[3]= _V_ghost_wall_0[3][i];
        
         V_U[3]= _V[3][i][0];
        
        
        //total enthalpy
        ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
        ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
        roe_boundary_F_eta(i, 0);
        */
        
        //   Ghost cells at wall 0


         //  ----- first ghost
           Delta_ghost_x_y = 0.5*(x_n[i+1][0] -x_n[i][0]);

           y_c_ghost= y_n[i][0]- Delta_ghost_x_y;

           x_c_ghost=x_n[i][0]+ Delta_ghost_x_y;

           _V_ghost_wall_0[0][i]=rho(1,x_c_ghost,y_c_ghost);
           
           _V_ghost_wall_0[1][i]=uvel(1,x_c_ghost,y_c_ghost);
           _V_ghost_wall_0[2][i]=vvel(1,x_c_ghost,y_c_ghost);
           //   _V_ghost_inflow[2][j]=0;
           _V_ghost_wall_0[3][i]=press(1,x_c_ghost,y_c_ghost);
        
        //second ghost
        Delta_ghost_x_y = 0.5*(x_n[i+1][0] -x_n[i][0]);

        y_c_ghost= y_n[i][0]- 2*Delta_ghost_x_y;

        x_c_ghost=x_n[i][0]+ Delta_ghost_x_y;

        _V_ghost_wall_0_2[0][i]=rho(1,x_c_ghost,y_c_ghost);
        
        _V_ghost_wall_0_2[1][i]=uvel(1,x_c_ghost,y_c_ghost);
        _V_ghost_wall_0_2[2][i]=vvel(1,x_c_ghost,y_c_ghost);
        //   _V_ghost_inflow[2][j]=0;
        _V_ghost_wall_0_2[3][i]=press(1,x_c_ghost,y_c_ghost);
        
        //calculate epsi at j=1
        
        
        
        //for density
        double denm= copysign(max(abs(_V[0][i][1]-_V[0][i][0]),1e-6),_V[0][i][1]-_V[0][i][0]);
        r_plus= (_V[0][i][2]-_V[0][i][1])/(denm);
        r_minus= (_V[0][i][0]-_V_ghost_wall_0[0][i])/(denm);
        epsi_plus_y[0][i][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[0][i][1]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u- velocity
        denm= copysign(max(abs(_V[1][i][1]-_V[1][i][0]),1e-6),_V[1][i][1]-_V[1][i][0]);
        r_plus= (_V[1][i][2]-_V[1][i][1])/(denm);
        r_minus= (_V[1][i][0]-_V_ghost_wall_0[1][i])/(denm);
        epsi_plus_y[1][i][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[1][i][1]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v velocity
        denm= copysign(max(abs(_V[2][i][1]-_V[2][i][0]),1e-6),_V[2][i][1]-_V[2][i][0]);
        r_plus= (_V[2][i][2]-_V[2][i][1])/(denm);
        r_minus= (_V[2][i][0]-_V_ghost_wall_0[2][i])/(denm);
        epsi_plus_y[2][i][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[2][i][1]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][i][1]-_V[3][i][0]),1e-6),_V[3][i][1]-_V[3][i][0]);
        r_plus= (_V[3][i][2]-_V[3][i][1])/(denm);
        r_minus= (_V[3][i][0]-_V_ghost_wall_0[3][i])/(denm);
        epsi_plus_y[3][i][1]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[3][i][1]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        //calculate epsi at 0
        
        //for density
        denm= copysign(max(abs(_V[0][i][0]-_V_ghost_wall_0[0][i]),1e-6),_V[0][i][0]-_V_ghost_wall_0[0][i]);
        r_plus= (_V[0][i][1]-_V[0][i][0])/(denm);
        r_minus= (_V_ghost_wall_0[0][i]-_V_ghost_wall_0_2[0][i])/(denm);
        epsi_plus_y[0][i][0]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[0][i][0]= (r_minus+abs(r_minus))/(1+r_minus);
        //for u velocity
        denm= copysign(max(abs(_V[1][i][0]-_V_ghost_wall_0[1][i]),1e-6),_V[1][i][0]-_V_ghost_wall_0[1][i]);
        r_plus= (_V[1][i][1]-_V[1][i][0])/(denm);
        r_minus= (_V_ghost_wall_0[1][i]-_V_ghost_wall_0_2[1][i])/(denm);
        epsi_plus_y[1][i][0]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[1][i][0]= (r_minus+abs(r_minus))/(1+r_minus);
        //for v- velocity
        denm= copysign(max(abs(_V[2][i][0]-_V_ghost_wall_0[2][i]),1e-6),_V[2][i][0]-_V_ghost_wall_0[2][i]);
        r_plus= (_V[2][i][1]-_V[2][i][0])/(denm);
        r_minus= (_V_ghost_wall_0[2][i]-_V_ghost_wall_0_2[2][i])/(denm);
        epsi_plus_y[2][i][0]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[2][i][0]= (r_minus+abs(r_minus))/(1+r_minus);
        //for pressure
        denm= copysign(max(abs(_V[3][i][0]-_V_ghost_wall_0[3][i]),1e-6),_V[3][i][0]-_V_ghost_wall_0[3][i]);
        r_plus= (_V[3][i][1]-_V[3][i][0])/(denm);
        r_minus= (_V_ghost_wall_0[3][i]-_V_ghost_wall_0_2[3][i])/(denm);
        epsi_plus_y[3][i][0]= (r_plus+abs(r_plus))/(1+r_plus);
        epsi_minus_y[3][i][0]= (r_minus+abs(r_minus))/(1+r_minus);
        
        
        
        
       
        
        
        
        
        
        //for eta faces at cell j=0
         //density
        V_D[0]= _V_ghost_wall_0[0][i]+(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_0[0][i]-_V_ghost_wall_0_2[0][i]) + (1+kappa)*epsi_minus_y[0][i][0]*(_V[0][i][0]-_V_ghost_wall_0[0][i]));
       
        V_U[0]= _V[0][i][0]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][1]*(_V[0][i][1]-_V[0][i][0]) + (1+kappa)*epsi_plus_y[0][i][0]*(_V[0][i][0]-_V_ghost_wall_0[0][i])        );
         //x-velocity
        V_D[1]=  _V_ghost_wall_0[1][i]+(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_0[1][i]-_V_ghost_wall_0_2[1][i]) + (1+kappa)*epsi_minus_y[1][i][0]*(_V[1][i][0]-_V_ghost_wall_0[1][i]));;
       
        V_U[1]= _V[1][i][0]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][1]*(_V[1][i][1]-_V[1][i][0]) + (1+kappa)*epsi_plus_y[1][i][0]*(_V[1][i][0]-_V_ghost_wall_0[1][i])        );
         //y-velocity
        V_D[2]=  _V_ghost_wall_0[2][i]+(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_0[2][i]-_V_ghost_wall_0_2[2][i]) + (1+kappa)*epsi_minus_y[2][i][0]*(_V[2][i][0]-_V_ghost_wall_0[2][i]));;
       
        V_U[2]= _V[2][i][0]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][1]*(_V[2][i][1]-_V[2][i][0]) + (1+kappa)*epsi_plus_y[2][i][0]*(_V[2][i][0]-_V_ghost_wall_0[2][i])        );
         //pressure
        V_D[3]=  _V_ghost_wall_0[3][i]+(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_0[3][i]-_V_ghost_wall_0_2[3][i]) + (1+kappa)*epsi_minus_y[3][i][0]*(_V[3][i][0]-_V_ghost_wall_0[3][i]));;
       
        V_U[3]= _V[3][i][0]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][1]*(_V[3][i][1]-_V[3][i][0]) + (1+kappa)*epsi_plus_y[3][i][0]*(_V[3][i][0]-_V_ghost_wall_0[3][i])        );
         //total enthalpy
         ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
         ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
         vanleer_eta(i, 0);
      
        /*
           //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
           _V_ghost_inflow[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

          // ------second ghost
           Delta_ghost_x_y = 0.5*(y_n[0][j+1]-y_n[0][j]);

           y_c_ghost= y_n[0][j]+ Delta_ghost_x_y;

           x_c_ghost= x_c_ghost - Delta_ghost_x_y;
          
           _V_ghost_inflow_2[0][j]=rho(1,x_c_ghost,y_c_ghost);
           
           _V_ghost_inflow_2[1][j]=uvel(1,x_c_ghost,y_c_ghost);
           _V_ghost_inflow_2[2][j]=vvel(1,x_c_ghost,y_c_ghost);
           //   _V_ghost_inflow[2][j]=0;
           _V_ghost_inflow_2[3][j]=press(1,x_c_ghost,y_c_ghost);
           
           //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
           _V_ghost_inflow_2[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

           
          
           
           
           
         //  -----third ghost
           Delta_ghost_x_y = 0.5*(y_n[0][j+1]- y_n[0][j]);

           y_c_ghost= y_n[0][j]+ Delta_ghost_x_y;

           x_c_ghost= x_c_ghost - Delta_ghost_x_y;

           _V_ghost_inflow_3[0][j]=rho(1,x_c_ghost,y_c_ghost);
           
           _V_ghost_inflow_3[1][j]=uvel(1,x_c_ghost,y_c_ghost);
           _V_ghost_inflow_3[2][j]=vvel(1,x_c_ghost,y_c_ghost);
           //   _V_ghost_inflow[2][j]=0;
           _V_ghost_inflow_3[3][j]=press(1,x_c_ghost,y_c_ghost);
           
           //_V_ghost_inflow[4]=2*et[_imin]-et[_imin];  //et
           _V_ghost_inflow_3[5][j]= _V_ghost_inflow[3][j]/( _V_ghost_inflow[0][j]*R_air); //T

      
 */
        
        //density

     
      
        
        
       //for eta faces at cell j=1
        //density
        V_D[0]= _V[0][i][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][0]*(_V[0][i][0]-_V_ghost_wall_0[0][i]) + (1+kappa)*epsi_minus_y[0][i][1]*(_V[0][i][1]-_V[0][i][0])        );
      
        V_U[0]= _V[0][i][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][2]*(_V[0][i][2]-_V[0][i][1]) + (1+kappa)*epsi_plus_y[0][i][1]*(_V[0][i][1]-_V[0][i][0])        );
        //x-velocity
        V_D[1]= _V[1][i][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][0]*(_V[1][i][0]-_V_ghost_wall_0[1][i]) + (1+kappa)*epsi_minus_y[1][i][1]*(_V[1][i][1]-_V[1][i][0])        );
      
        V_U[1]= _V[1][i][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][2]*(_V[1][i][2]-_V[1][i][1]) + (1+kappa)*epsi_plus_y[1][i][1]*(_V[1][i][1]-_V[1][i][0])        );
        //y-velocity
        V_D[2]= _V[2][i][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][0]*(_V[2][i][0]-_V_ghost_wall_0[2][i]) + (1+kappa)*epsi_minus_y[2][i][1]*(_V[2][i][1]-_V[2][i][0])        );
      
        V_U[2]= _V[2][i][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][2]*(_V[2][i][2]-_V[2][i][1]) + (1+kappa)*epsi_plus_y[2][i][1]*(_V[2][i][1]-_V[2][i][0])        );
        //pressure
        V_D[3]= _V[3][i][0]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][0]*(_V[3][i][0]-_V_ghost_wall_0[3][i]) + (1+kappa)*epsi_minus_y[3][i][1]*(_V[3][i][1]-_V[3][i][0])        );
      
        V_U[3]= _V[3][i][1]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][2]*(_V[3][i][2]-_V[3][i][1]) + (1+kappa)*epsi_plus_y[3][i][1]*(_V[3][i][1]-_V[3][i][0])        );
        //total enthalpy
        ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
        ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
        vanleer_eta(i, 1);
        
        
        //for eta faces at cell im_c
        
        //   Ghost cells at wall jmax


         //  ----- first ghost
           Delta_ghost_x_y = 0.5*(x_n[i+1][i_n_max-1] -x_n[i][i_n_max-1]);

           y_c_ghost= y_n[i][i_n_max-1]+ Delta_ghost_x_y;

           x_c_ghost=x_n[i][i_n_max-1]+ Delta_ghost_x_y;

           _V_ghost_wall_jmax[0][i]=rho(1,x_c_ghost,y_c_ghost);
           
           _V_ghost_wall_jmax[1][i]=uvel(1,x_c_ghost,y_c_ghost);
           _V_ghost_wall_jmax[2][i]=vvel(1,x_c_ghost,y_c_ghost);
           //   _V_ghost_inflow[2][j]=0;
           _V_ghost_wall_jmax[3][i]=press(1,x_c_ghost,y_c_ghost);
        //  ----- 2nd ghost
          Delta_ghost_x_y = 0.5*(x_n[i+1][i_n_max-1] -x_n[i][i_n_max-1]);

          y_c_ghost= y_n[i][i_n_max-1]+ 2*Delta_ghost_x_y;

          x_c_ghost=x_n[i][i_n_max-1]+ Delta_ghost_x_y;

          _V_ghost_wall_jmax_2[0][i]=rho(1,x_c_ghost,y_c_ghost);
          
          _V_ghost_wall_jmax_2[1][i]=uvel(1,x_c_ghost,y_c_ghost);
          _V_ghost_wall_jmax_2[2][i]=vvel(1,x_c_ghost,y_c_ghost);
          //   _V_ghost_inflow[2][j]=0;
          _V_ghost_wall_jmax_2[3][i]=press(1,x_c_ghost,y_c_ghost);
      
        
        
        
        
        
            
            //calculate epsi at imf-1 or imc
            
            
            //for density
            denm= copysign(max(abs(_V[0][i][im_c]-_V[0][i][im_c-1]),1e-6),_V[0][i][im_c]-_V[0][i][im_c-1]);
            r_plus= (_V_ghost_wall_jmax[0][i]-_V[0][i][im_c])/(denm);
            r_minus= (_V[0][i][im_c-1]-_V[0][i][im_c-2])/(denm);
            epsi_plus_y[0][i][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[0][i][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
            //for u-velocity
            denm= copysign(max(abs(_V[1][i][im_c]-_V[1][i][im_c-1]),1e-6),_V[1][i][im_c]-_V[1][i][im_c-1]);
            r_plus= (_V_ghost_wall_jmax[1][i]-_V[1][i][im_c])/(denm);
            r_minus= (_V[1][i][im_c-1]-_V[1][i][im_c-2])/(denm);
            epsi_plus_y[1][i][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[1][i][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
            //for v-velocity
            denm= copysign(max(abs(_V[2][i][im_c]-_V[2][i][im_c-1]),1e-6),_V[2][i][im_c]-_V[2][i][im_c-1]);
            r_plus= (_V_ghost_wall_jmax[2][i]-_V[2][i][im_c])/(denm);
            r_minus= (_V[2][i][im_c-1]-_V[2][i][im_c-2])/(denm);
            epsi_plus_y[2][i][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[2][i][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
            //for pressure
            denm= copysign(max(abs(_V[3][i][im_c]-_V[3][i][im_c-1]),1e-6),_V[3][i][im_c]-_V[3][i][im_c-1]);
            r_plus= (_V_ghost_wall_jmax[3][i]-_V[3][i][im_c])/(denm);
            r_minus= (_V[3][i][im_c-1]-_V[3][i][im_c-2])/(denm);
            epsi_plus_y[3][i][im_c]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[3][i][im_c]= (r_minus+abs(r_minus))/(1+r_minus);
            
            
            
            
            
            //calculate epsi for imf or imc+1
            
            
            //for density
            denm= copysign(max(abs(_V_ghost_wall_jmax[0][i]-_V[0][i][im_c]),1e-6),_V_ghost_wall_jmax[0][i]-_V[0][i][im_c]);
            r_plus= (_V_ghost_wall_jmax_2[0][i]-_V_ghost_wall_jmax[0][i])/(denm);
            r_minus= (_V[0][i][im_c]-_V[0][i][im_c-1])/(denm);
            epsi_plus_y[0][i][im_c+1]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[0][i][im_c+1]= (r_minus+abs(r_minus))/(1+r_minus);
            //for velocity
            denm= copysign(max(abs(_V_ghost_wall_jmax[1][i]-_V[1][i][im_c]),1e-6),_V_ghost_wall_jmax[1][i]-_V[1][i][im_c]);
            r_plus= (_V_ghost_wall_jmax_2[1][i]-_V_ghost_wall_jmax[1][i])/(denm);
            r_minus= (_V[1][i][im_c]-_V[1][i][im_c-1])/(denm);
            epsi_plus_y[1][i][im_c+1]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[1][i][im_c+1]= (r_minus+abs(r_minus))/(1+r_minus);
            //for v-velocity
            denm= copysign(max(abs(_V_ghost_wall_jmax[2][i]-_V[2][i][im_c]),1e-6),_V_ghost_wall_jmax[2][i]-_V[2][i][im_c]);
            r_plus= (_V_ghost_wall_jmax_2[2][i]-_V_ghost_wall_jmax[2][i])/(denm);
            r_minus= (_V[2][i][im_c]-_V[2][i][im_c-1])/(denm);
            epsi_plus_y[2][i][im_c+1]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[2][i][im_c+1]= (r_minus+abs(r_minus))/(1+r_minus);
            
            //for pressure
            denm= copysign(max(abs(_V_ghost_wall_jmax[3][i]-_V[3][i][im_c]),1e-6),_V_ghost_wall_jmax[3][i]-_V[3][i][im_c]);
            r_plus= (_V_ghost_wall_jmax_2[3][i]-_V_ghost_wall_jmax[3][i])/(denm);
            r_minus= (_V[3][i][im_c]-_V[3][i][im_c-1])/(denm);
            epsi_plus_y[3][i][im_c+1]= (r_plus+abs(r_plus))/(1+r_plus);
            epsi_minus_y[3][i][im_c+1]= (r_minus+abs(r_minus))/(1+r_minus);
            
      
        
      
        
        //density
       
        
      
        
        
        //for eta faces at imf-1
        
     
         //density
        V_D[0]=_V[0][i][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][im_c-1]*(_V[0][i][im_c-1]-_V[0][i][im_c-2]) + (1+kappa)*epsi_minus_y[0][i][im_c]*(_V[0][i][im_c]-_V[0][i][im_c-1])        );
       
  
        V_U[0]=  _V[0][i][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][im_c+1]*(_V_ghost_wall_jmax[0][i]-_V[0][i][im_c]) + (1+kappa)*epsi_plus_y[0][i][im_c]*(_V[0][i][im_c]-_V[0][i][im_c-1])        );
         //x-velocity
        V_D[1]= _V[1][i][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][im_c-1]*(_V[1][i][im_c-1]-_V[1][i][im_c-2]) + (1+kappa)*epsi_minus_y[1][i][im_c]*(_V[1][i][im_c]-_V[1][i][im_c-1])        );
       
        V_U[1]=  _V[1][i][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][im_c+1]*(_V_ghost_wall_jmax[1][i]-_V[1][i][im_c]) + (1+kappa)*epsi_plus_y[1][i][im_c]*(_V[1][i][im_c]-_V[1][i][im_c-1])        );
         //y-velocity
        V_D[2]= _V[2][i][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][im_c-1]*(_V[2][i][im_c-1]-_V[2][i][im_c-2]) + (1+kappa)*epsi_minus_y[2][i][im_c]*(_V[2][i][im_c]-_V[2][i][im_c-1])        );
       
        V_U[2]=_V[2][i][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][im_c+1]*(_V_ghost_wall_jmax[2][i]-_V[2][i][im_c]) + (1+kappa)*epsi_plus_y[2][i][im_c]*(_V[2][i][im_c]-_V[2][i][im_c-1])        );
         //pressure
        V_D[3]= _V[3][i][im_c-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][im_c-1]*(_V[3][i][im_c-1]-_V[3][i][im_c-2]) + (1+kappa)*epsi_minus_y[3][i][im_c]*(_V[3][i][im_c]-_V[3][i][im_c-1])        );
       
        V_U[3]=_V[3][i][im_c]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][im_c+1]*(_V_ghost_wall_jmax[3][i]-_V[3][i][im_c]) + (1+kappa)*epsi_plus_y[3][i][im_c]*(_V[3][i][im_c]-_V[3][i][im_c-1])        );
         //total enthalpy
         ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
         ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
         vanleer_eta(i, im_f-1);
        
        
        
        //zero flux implementation
        /*
        P_wall =_V[3][i][im_c];// - 0.5*(_V[3][i][im_c-1]-_V[3][i][im_c]);
        
        
        a4=face_area(x_n[i+1][im_f], x_n[i][im_f],y_n[i+1][im_f], y_n[i][im_f]);
        nx= n_x(y_n[i+1][im_f], y_n[i][im_f], a4);
        ny= n_y(x_n[i+1][im_f], x_n[i][im_f], a4);
        
        _F_eta[0][i][im_f]=0;
        
        _F_eta[1][i][im_f]=(nx*P_wall);
        
        _F_eta[2][i][im_f]=(ny*P_wall);
        
        _F_eta[3][i][im_f]=0;
        */
        // ghost cells implementation
       /*
        double a3=face_area(x_n[i+1][im_f], x_n[i][im_f],y_n[i+1][im_f], y_n[i][im_f]);
        nx= n_x(y_n[i+1][im_f], y_n[i][im_f], a3);
        ny= n_y(x_n[i+1][im_f], x_n[i][im_f], a3);
        
        A=_V[1][i][im_f]*ny+_V[2][i][im_f]*nx;
        
        B= -_V[1][i][im_f]*nx-_V[2][i][im_f]*ny;
        
        
        
        
        
        
        _V_ghost_wall_jmax[2][i]=(B*ny-A*nx     )/(  ny*ny-nx*nx );
        
        _V_ghost_wall_jmax[1][i]=(A-_V_ghost_wall_jmax[2][i]*nx    )/(  ny );
        
        _V_ghost_wall_jmax[5][i]=T[i][im_c];
        
        
        _V_ghost_wall_jmax[3][i]=2*_V[3][i][im_c]-_V[3][i][im_c-1];
        
        _V_ghost_wall_jmax[0][i]= _V_ghost_wall_jmax[3][i]/( _V_ghost_wall_jmax[5][i]*R_air);
      
       
        
        //density
        V_D[0]= _V[0][i][im_c];
      
        V_U[0]= _V_ghost_wall_jmax[0][i];
        //x-velocity
           V_D[1]= _V[1][i][im_c];
      
        V_U[1]= _V_ghost_wall_jmax[1][i];
        //y-velocity
           V_D[2]= _V[2][i][im_c];
      
        V_U[2]= _V_ghost_wall_jmax[2][i];
        //pressure
           V_D[3]= _V[3][i][im_c];
      
        V_U[3]= _V_ghost_wall_jmax[3][i];
   
    //total enthalpy
    ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
    ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
        
    roe_boundary_F_eta(i, im_f);
   */
        // for cell at j = im_c -1
       
      
       
        
        
        
        //density
        
        //density
     
        V_D[0]= _V[0][i][im_c]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][im_c]*(_V[0][i][im_c]-_V[0][i][im_c-1]) + (1+kappa)*epsi_minus_y[0][i][im_c+1]*(_V_ghost_wall_jmax[0][i]-_V[0][i][im_c])        );
      
        V_U[0]= _V_ghost_wall_jmax[0][i]-(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_jmax_2[0][i]-_V_ghost_wall_jmax[0][i]) + (1+kappa)*epsi_plus_y[0][i][im_c+1]*(_V_ghost_wall_jmax[0][i]-_V[0][i][im_c])        );
        //x-velocity
        V_D[1]= _V[1][i][im_c]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][im_c]*(_V[1][i][im_c]-_V[1][i][im_c-1]) + (1+kappa)*epsi_minus_y[1][i][im_c+1]*(_V_ghost_wall_jmax[1][i]-_V[1][i][im_c])        );
      
        V_U[1]= _V_ghost_wall_jmax[1][i]-(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_jmax_2[1][i]-_V_ghost_wall_jmax[1][i]) + (1+kappa)*epsi_plus_y[1][i][im_c+1]*(_V_ghost_wall_jmax[1][i]-_V[1][i][im_c])        );
        //y-velocity
        V_D[2]= _V[2][i][im_c]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][im_c]*(_V[2][i][im_c]-_V[2][i][im_c-1]) + (1+kappa)*epsi_minus_y[2][i][im_c+1]*(_V_ghost_wall_jmax[2][i]-_V[2][i][im_c])        );
      
        V_U[2]= _V_ghost_wall_jmax[2][i]-(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_jmax_2[2][i]-_V_ghost_wall_jmax[2][i]) + (1+kappa)*epsi_plus_y[2][i][im_c+1]*(_V_ghost_wall_jmax[2][i]-_V[2][i][im_c])        );
        //pressure
        V_D[3]= _V[3][i][im_c]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][im_c]*(_V[3][i][im_c]-_V[3][i][im_c-1]) + (1+kappa)*epsi_minus_y[3][i][im_c+1]*(_V_ghost_wall_jmax[3][i]-_V[3][i][im_c])        );
      
        V_U[3]= _V_ghost_wall_jmax[3][i]-(epsilon_upwind/4)*((1-kappa)*(_V_ghost_wall_jmax_2[3][i]-_V_ghost_wall_jmax[3][i]) + (1+kappa)*epsi_plus_y[3][i][im_c+1]*(_V_ghost_wall_jmax[3][i]-_V[3][i][im_c])        );
        //total enthalpy
        ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
        ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
      
      vanleer_eta(i, im_f);
        
        
        
       
        
        
    }
    
    
    
}

void Dv_nozzle::euler_explicit()
{
    double a1, a2,a3,a4;
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            //a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            
            
            
            
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
           /*
            double n_x_i,n_x_iplus,n_y_i,n_y_iplus;
            double n_x_j,n_x_jplus,n_y_j,n_y_jplus;
            
            n_x_iplus=n_x(y_n[i+1][j+1],y_n[i+1][j],a1);
            n_y_iplus=n_y(x_n[i+1][j+1],x_n[i+1][j],a1);
            n_x_i=n_x(y_n[i][j+1],y_n[i][j],a1);
            n_y_i=n_y(x_n[i][j+1],x_n[i][j],a1);
            
            n_x_jplus=n_x(y_n[i+1][j+1],y_n[i][j+1],a1);
            n_y_jplus=n_y(x_n[i+1][j+1],x_n[i][j+1],a1);
            n_x_j=n_x(y_n[i+1][j],y_n[i][j],a1);
            n_y_j=n_y(x_n[i+1][j],x_n[i][j],a1);
            */
            
            
            
           
            
            //continuity eqn
            
            _U[0][i][j]=((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4)*(delta_t_gl/(cell_vol(i,j))) + _U[0][i][j] + rmassconv(1, x_c[i][j], y_c[i][j])*(delta_t_gl);
            
            
            //x-mtm
            
            _U[1][i][j]=((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4)*(delta_t_gl/(cell_vol(i,j))) + _U[1][i][j] + xmtmconv(1, x_c[i][j], y_c[i][j])*(delta_t_gl);
            
            //y-mtm
            _U[2][i][j]=((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4)*(delta_t_gl/(cell_vol(i,j))) + _U[2][i][j] + ymtmconv(1, x_c[i][j], y_c[i][j])*(delta_t_gl);
            
            
            //energy eqn
            
            _U[3][i][j]=((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4)*(delta_t_gl/(cell_vol(i,j))) + _U[3][i][j] + energyconv(1, x_c[i][j], y_c[i][j])*(delta_t_gl);
            
        //    cout<< "scell vol at cell  ("<<i<<","<<j<<") = "<<cell_vol(x_c[i][j], y_c[i][j])<<endl;
      //   cout<< "source term for ymtm eqn at cell  ("<<i<<","<<j<<") = "<<ymtmconv(1, x_c[i][j], y_c[i][j])<<endl;
            
        }
    }
    
    for(int i=0;i<=im_f;i++)
    {
        for(int j=0;j<=im_f-1;j++)
        {
            
            
       // cout<< "F_kasi at  face("<<i<<","<<j<<") = "<<_F[0][i][j]<<endl;
            
        }
    }
    
    for (int j=0;j<=im_c;j++)
    {
       
        for (int i=0;i<=im_c;i++)
        {
            
            _V[0][i][j]=_U[0][i][j];
            _V[1][i][j]= _U[1][i][j]/_V[0][i][j];
            _V[2][i][j]= _U[2][i][j]/_V[0][i][j];
            //  if (_V[1][i]<0)_V[1][i]=1;
            et[i][j] = _U[3][i][j]/_V[0][i][j];
              if (et[i][j]<200)et[i][j]=200;
            T[i][j] = (et[i][j] - 0.5* (_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j] ))*((gamma-1)/R_air);
            if (T[i][j]<50)T[i][j]=50;
            
            
            _V[3][i][j]=_V[0][i][j]*R_air*T[i][j];
              if (_V[3][i][j]<3000)_V[3][i][j]=3000;
            ht[i][j]=(gamma/(gamma-1))*(_V[3][i][j]/_V[0][i][j]) + (_V[2][i][j]*_V[2][i][j]+_V[1][i][j]*_V[1][i][j])*0.5;
            _mach[i][j]=sqrt(_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j])/sqrt(gamma*R_air*abs(T[i][j]));
            
            
        // cout<< "flux at  cell ("<<i<<","<<j<<") = "<<_F_eta[0][j][i]<<endl;
       // cout<< "flux at  cell ("<<i<<","<<1<<") = "<<_F_eta[0][i][1]<<endl;
       //cout<< "density at  cell ("<<i<<","<<j<<") = "<<_V[0][i][j]<<endl;
            
        }
        
    }
        //calculate epsi- kasi direction kasi corresponds to x kinda
        
        for (int j=2;j<im_c-1;j++)
        {
           
           // if(std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])))<=1e-3)break;
          
            for (int i=2;i<im_c;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i-1][j]),1e-6),_V[0][i][j]-_V[0][i-1][j]);
                r_plus= (_V[0][i+1][j]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i-1][j]-_V[0][i-2][j])/(denm);
                epsi_plus_x[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i-1][j]),1e-6),_V[1][i][j]-_V[1][i-1][j]);
                r_plus= (_V[1][i+1][j]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i-1][j]-_V[1][i-2][j])/(denm);
                epsi_plus_x[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i-1][j]),1e-6),_V[2][i][j]-_V[2][i-1][j]);
                r_plus= (_V[2][i+1][j]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i-1][j]-_V[2][i-2][j])/(denm);
                epsi_plus_x[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i-1][j]),1e-6),_V[3][i][j]-_V[3][i-1][j]);
                r_plus= (_V[3][i+1][j]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i-1][j]-_V[3][i-2][j])/(denm);
                epsi_plus_x[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_x[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
         /*
                if(std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])))<=1e-7){
                    
                    epsi_minus_x[0][i][j]=1;
                    epsi_minus_x[1][i][j]=1;
                    epsi_minus_x[2][i][j]=1;
                    epsi_minus_x[3][i][j]=1;
                    
                    epsi_plus_x[0][i][j]=1;
                    epsi_plus_x[1][i][j]=1;
                    epsi_plus_x[2][i][j]=1;
                    epsi_plus_x[3][i][j]=1;
                    
               
                }
          */
         
            }
            
        }
        
        
        //calculate epsi- eta direction kasi corresponds to y kinda
        
        for (int j=2;j<im_c;j++)
        {
           // if(std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])))<=1e-3)break;
            for (int i=2;i<im_c-1;i++)
            {
                //for density
                double denm= copysign(max(abs(_V[0][i][j]-_V[0][i][j-1]),1e-6),_V[0][i][j]-_V[0][i][j-1]);
                r_plus= (_V[0][i][j+1]-_V[0][i][j])/(denm);
                r_minus= (_V[0][i][j-1]-_V[0][i][j-2])/(denm);
                epsi_plus_y[0][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[0][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for u velocity
                denm= copysign(max(abs(_V[1][i][j]-_V[1][i][j-1]),1e-6),_V[1][i][j]-_V[1][i][j-1]);
                r_plus= (_V[1][i][j+1]-_V[1][i][j])/(denm);
                r_minus= (_V[1][i][j-1]-_V[1][i][j-2])/(denm);
                epsi_plus_y[1][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[1][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                
                //for v velocity
                denm= copysign(max(abs(_V[2][i][j]-_V[2][i][j-1]),1e-6),_V[2][i][j]-_V[2][i][j-1]);
                r_plus= (_V[2][i][j+1]-_V[2][i][j])/(denm);
                r_minus= (_V[2][i][j-1]-_V[2][i][j-2])/(denm);
                epsi_plus_y[2][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[2][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
                //for pressure
                denm= copysign(max(abs(_V[3][i][j]-_V[3][i][j-1]),1e-6),_V[3][i][j]-_V[3][i][j-1]);
                r_plus= (_V[3][i][j+1]-_V[3][i][j])/(denm);
                r_minus= (_V[3][i][j-1]-_V[3][i][j-2])/(denm);
                epsi_plus_y[3][i][j]= (r_plus+abs(r_plus))/(1+r_plus);
                epsi_minus_y[3][i][j]= (r_minus+abs(r_minus))/(1+r_minus);
           
                /*
                if(std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])))<=1e-7){
                    
                    epsi_minus_y[0][i][j]=1;
                    epsi_minus_y[1][i][j]=1;
                    epsi_minus_y[2][i][j]=1;
                    epsi_minus_y[3][i][j]=1;
                    
                    epsi_plus_y[0][i][j]=1;
                    epsi_plus_y[1][i][j]=1;
                    epsi_plus_y[2][i][j]=1;
                    epsi_plus_y[3][i][j]=1;
                    
                  
                }
            */
            
            }
            
        }
   
    
    
    
    
    vanleer();
    
    
    
    
    
}
double Dv_nozzle::area(double x)
{
    _area = 0.2 + 0.4*(1 + sin(3.14159265359*(x-0.5)));
    
    
    return _area;
    
}
double Dv_nozzle::d_area(double x)
{
    _d_area = 0.4*3.14159265359*cos(3.14159265359*(x-0.5));
    
    
    return _d_area;
    
}
void Dv_nozzle::time_step()
{
     double a1, a2, a3, a4;
                
     for (int i=0;i<=im_c;i++)
         for (int j =0;j<=im_c;j++)
         {
             
             
             
              a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
              a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
              a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
              a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
             
             
             
             
             
         //for two x faces
             
             //for kase
             nx=n_x(y_n[i+1][j+1], y_n[i+1][j],a1);
             
             nx= 0.5*(nx + n_x(y_n[i][j+1],y_n[i][j],a2));
             
             
             ny=n_y(x_n[i+1][j+1], x_n[i+1][j], a1);
             
             ny= 0.5*(ny + n_y(x_n[i][j+1], x_n[i][j], a2));
             a[i][j]=sqrt(abs((gamma-1)*(ht[i][j]-0.5*(_V[1][i][j]*_V[1][i][j]+_V[2][i][j]*_V[2][i][j]))));
             
             
            
             eigen_max[0]=abs(_V[1][i][j]*nx+_V[2][i][j]*ny) + a[i][j];
             
             // for eta
             nx=-n_x(y_n[i+1][j+1], y_n[i][j+1], a3);
             
             nx= 0.5*(nx  -n_x(y_n[i+1][j],y_n[i][j],a4));
             
             
             ny=-n_y(x_n[i+1][j+1], x_n[i][j+1], a3);
             
             ny= 0.5*(ny - n_y(x_n[i+1][j], x_n[i][j], a4));
             
             
             eigen_max[1]=abs(_V[1][i][j]*nx+_V[2][i][j]*ny) + a[i][j];
             
             
             
             delta_t[i][j] =( CFL*cell_vol(i,j)) /(eigen_max[0]*cell_vol_kasi(a1, a2)  + eigen_max[1]*cell_vol_eta(a3, a4)   );
             
             
             
             delta_t_gl=min(delta_t[i][j],delta_t_gl);
             
         //cout<<"delta t  is "<<delta_t[i][j]<<endl;
         //    cout<<"cell vol is "<<cell_vol(i,j)<<endl;
             
         }
     
     
 }
void Dv_nozzle::DE()
{
    double DE_1 =0;
    double DE_inf=-99999;
    for (int i= 0;i<=im_c;i++)
    {
        for (int j=0;j<=im_c;j++)
       
        
        {
            
            DE_1 = DE_1 + abs( rho(1, x_c[i][j], y_c[i][j])-_V[0][i][j]    );
            DE_inf = max(DE_inf,abs( rho(1, x_c[i][j], y_c[i][j])-_V[0][i][j]    ));
            
            
            
        }
            
            

        
     
    }
    
    DE_1 = DE_1/(_imax_c*_imax_c);
    
    cout<<"L1 norm is "<<DE_1<<endl;
    cout<<"L_inf norm is "<<DE_inf<<endl;
    
}
double Dv_nozzle::L2norm(int iter){
    
    
    rL2norm[0]=0;
    rL2norm[1]=0;
    rL2norm[2]=0;
    rL2norm[3]=0;
    
    double a1, a2,a3,a4;
    
    
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
        
            _R_iter[0]=(((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4)  +rmassconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j))      );
            
            rL2norm[0] = rL2norm[0] +_R_iter[0]*_R_iter[0];
            
            
            _R_iter[1]=((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4 + xmtmconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j))  )     ;
            
            rL2norm[1] = rL2norm[1] +  _R_iter[1]* _R_iter[1];
            
            _R_iter[2]=(((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4)+ymtmconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j)) )  ;
            
            rL2norm[2] = rL2norm[2] +  _R_iter[2]* _R_iter[2];
            
            _R_iter[3]=(((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4)+energyconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j)))  ;
            
            rL2norm[3] = rL2norm[3] +  _R_iter[3]* _R_iter[3];
            
           // std::cout<<"for cell "<<i<<","<<j<< " R eqn1 "<< _R_iter[0]<<std::endl;
            //std::cout<<"for cell "<<i<<","<<j<< " R eqn2 "<< _R_iter[1]<<std::endl;
           //  std::cout<<"for cell "<<i<<","<<j<< " R eqn3 "<< _R_iter[2]<<std::endl;
             
            
        }
    }
    double imax_c = _imax_c*1.0;
    
    rL2norm[0]=sqrt(rL2norm[0]/(imax_c*imax_c))/_rL2initial[0];
    rL2norm[1]=sqrt(rL2norm[1]/(imax_c*imax_c))/_rL2initial[1];
    rL2norm[2]=sqrt(rL2norm[2]/(imax_c*imax_c))/_rL2initial[2];
    rL2norm[3]=sqrt(rL2norm[3]/(imax_c*imax_c))/_rL2initial[3];
    
    /*
     std::cout<< " R eqn1 "<< _R_iter[0]<<std::endl;
     std::cout<< " R eqn2 "<< _R_iter[1]<<std::endl;
     std::cout<< " R eqn3 "<< _R_iter[2]<<std::endl;
     std::cout<< " L eqn1 "<< rL2norm[0]<<std::endl;
     std::cout<< " L eqn2 "<< rL2norm[1]<<std::endl;
     std::cout<< " L eqn3 "<< rL2norm[2]<<std::endl;
     std::cout<< " L inti eqn1 "<< _rL2initial[0]<<std::endl;
     std::cout<< " L inti eqn2 "<< _rL2initial[1]<<std::endl;
     std::cout<< " L inti eqn3 "<< _rL2initial[2]<<std::endl;
     std::cout<< " u is  "<< _V[1][9]<<std::endl;
     */
    if(iter%10==0)L_vs_Iter<<std::setprecision(5)<<iter<<","<<std::setprecision(5)<<rL2norm[0]<<","<<std::setprecision(5)<<rL2norm[1]<<","<<std::setprecision(5)<<rL2norm[2]<<std::endl;
    
 return std::max(rL2norm[3],max(rL2norm[0],std::max(rL2norm[1],rL2norm[2])));
    
    
}
void Dv_nozzle::rL2initial()


{
    double a1, a2,a3,a4;
    _rL2initial[0]=0;
    _rL2initial[1]=0;
    _rL2initial[2]=0;
    _rL2initial[3]=0;
    for (int j=0;j<=im_c;j++)
    {
        for (int i=0;i<=im_c;i++)
        {
            
            
            
            
            
            a1=face_area(x_n[i+1][j+1], x_n[i+1][j],y_n[i+1][j+1], y_n[i+1][j]);
            
            a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
            a3=face_area(x_n[i+1][j+1], x_n[i][j+1],y_n[i+1][j+1], y_n[i][j+1]);
            a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
        
            _R[0]= (((-_F[0][i+1][j])*a1+a2*_F[0][i][j]-a3*_F_eta[0][i][j+1]+(_F_eta[0][i][j])*a4)+rmassconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j))  )  ;
            
            
            _rL2initial[0] = _rL2initial[0] + _R[0]*_R[0];
            
            
            _R[1]= (((-_F[1][i+1][j])*a1+a2*_F[1][i][j]-a3*_F_eta[1][i][j+1]+(_F_eta[1][i][j])*a4)+xmtmconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j))  )  ;
            
            _rL2initial[1] = _rL2initial[1] + _R[1]*_R[1];
            
            _R[2]= (((-_F[2][i+1][j])*a1+a2*_F[2][i][j]-a3*_F_eta[2][i][j+1]+(_F_eta[2][i][j])*a4)+ymtmconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j)) )  ;
            
            _rL2initial[2] = _rL2initial[2] + _R[2]*_R[2];
            _R[3]= (((-_F[3][i+1][j])*a1+a2*_F[3][i][j]-a3*_F_eta[3][i][j+1]+(_F_eta[3][i][j])*a4)+energyconv(1, x_c[i][j], y_c[i][j])*(cell_vol(i,j)))  ;
            
            _rL2initial[3] = _rL2initial[3] + _R[3]*_R[3];
            
           
            std::cout<<" cell  ("<<i<<","<<j<<")"<<std::endl;
            std::cout<<" r1 initial is "<<_rL2initial[0]<<std::endl;
            std::cout<<" r2 initial is "<<_rL2initial[1]<<std::endl;
            std::cout<<" r3 initial is "<<_rL2initial[2]<<std::endl;
            std::cout<<" r4 initial is "<<_rL2initial[3]<<std::endl;
            
        }
    }
    double imax_c = _imax_c*1.0;
    _rL2initial[0] = sqrt(_rL2initial[0]/(imax_c*imax_c));
    _rL2initial[1] = sqrt(_rL2initial[1]/(imax_c*imax_c));
    _rL2initial[2] = sqrt(_rL2initial[2]/(imax_c*imax_c));
    _rL2initial[3] = sqrt(_rL2initial[3]/(imax_c*imax_c));
    
    std::cout<<" r initial is "<<_rL2initial[0]<<std::endl;
    
    
}
void Dv_nozzle::vanleer()
{
    
    //calculate the flux for kasi
for (int j=0;j<=im_c;j++){
       for (int i=2;i<=im_f-2;i++)
       {
           
           
           //density
           V_L[0]= _V[0][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[0][i-1][j]*(_V[0][i-1][j]-_V[0][i-2][j]) + (1+kappa)*epsi_minus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
           
           V_R[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[0][i+1][j]*(_V[0][i+1][j]-_V[0][i][j]) + (1+kappa)*epsi_plus_x[0][i][j]*(_V[0][i][j]-_V[0][i-1][j])        );
           //x-velocity
           V_L[1]= _V[1][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[1][i-1][j]*(_V[1][i-1][j]-_V[1][i-2][j]) + (1+kappa)*epsi_minus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
           
           V_R[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[1][i+1][j]*(_V[1][i+1][j]-_V[1][i][j]) + (1+kappa)*epsi_plus_x[1][i][j]*(_V[1][i][j]-_V[1][i-1][j])        );
           //y-velocity
           V_L[2]= _V[2][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[2][i-1][j]*(_V[2][i-1][j]-_V[2][i-2][j]) + (1+kappa)*epsi_minus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
           
           V_R[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[2][i+1][j]*(_V[2][i+1][j]-_V[2][i][j]) + (1+kappa)*epsi_plus_x[2][i][j]*(_V[2][i][j]-_V[2][i-1][j])        );
           //pressure
           V_L[3]= _V[3][i-1][j]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_x[3][i-1][j]*(_V[3][i-1][j]-_V[3][i-2][j]) + (1+kappa)*epsi_minus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
           
           V_R[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_x[3][i+1][j]*(_V[3][i+1][j]-_V[3][i][j]) + (1+kappa)*epsi_plus_x[3][i][j]*(_V[3][i][j]-_V[3][i-1][j])        );
           //total enthalpy
           ht_L=(gamma/(gamma-1))*(V_L[3]/V_L[0]) + ((V_L[1]*V_L[1])+(V_L[2]*V_L[2]))*0.5;
           ht_R=(gamma/(gamma-1))*(V_R[3]/V_R[0]) + ((V_R[1]*V_R[1])+(V_R[2]*V_R[2]))*0.5;
           
           vanleer_kasi(i, j);
       
       }
   }
    
    
    for (int j=2;j<=im_f-2;j++){
        for (int i=0;i<=im_c;i++)
        {
            
            // flux for eta faces
            
            
            //density
            V_D[0]= _V[0][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[0][i][j-1]*(_V[0][i][j-1]-_V[0][i][j-2]) + (1+kappa)*epsi_minus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
            
            V_U[0]= _V[0][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[0][i][j+1]*(_V[0][i][j+1]-_V[0][i][j]) + (1+kappa)*epsi_plus_y[0][i][j]*(_V[0][i][j]-_V[0][i][j-1])        );
            //x-velocity
            V_D[1]= _V[1][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[1][i][j-1]*(_V[1][i][j-1]-_V[1][i][j-2]) + (1+kappa)*epsi_minus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
            
            V_U[1]= _V[1][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[1][i][j+1]*(_V[1][i][j+1]-_V[1][i][j]) + (1+kappa)*epsi_plus_y[1][i][j]*(_V[1][i][j]-_V[1][i][j-1])        );
            //y-velocity
            V_D[2]= _V[2][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[2][i][j-1]*(_V[2][i][j-1]-_V[2][i][j-2]) + (1+kappa)*epsi_minus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
            
            V_U[2]= _V[2][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[2][i][j+1]*(_V[2][i][j+1]-_V[2][i][j]) + (1+kappa)*epsi_plus_y[2][i][j]*(_V[2][i][j]-_V[2][i][j-1])        );
            //pressure
            V_D[3]= _V[3][i][j-1]+(epsilon_upwind/4)*((1-kappa)*epsi_plus_y[3][i][j-1]*(_V[3][i][j-1]-_V[3][i][j-2]) + (1+kappa)*epsi_minus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
            
            V_U[3]= _V[3][i][j]-(epsilon_upwind/4)*((1-kappa)*epsi_minus_y[3][i][j+1]*(_V[3][i][j+1]-_V[3][i][j]) + (1+kappa)*epsi_plus_y[3][i][j]*(_V[3][i][j]-_V[3][i][j-1])        );
            //total enthalpy
            ht_D=(gamma/(gamma-1))*(V_D[3]/V_D[0]) + ((V_D[1]*V_D[1])+(V_D[2]*V_D[2]))*0.5;
            ht_U=(gamma/(gamma-1))*(V_U[3]/V_U[0]) + ((V_U[1]*V_U[1])+(V_U[2]*V_U[2]))*0.5;
            
            
            
           
            
            vanleer_eta(i, j);
            
        }
        
    }
    
}
void Dv_nozzle::vanleer_kasi (int i,int j)

{
    
    
    //flux in kasi direction
    
     
    double a2=face_area(x_n[i][j+1], x_n[i][j],y_n[i][j+1], y_n[i][j]);
    nx=n_x(y_n[i][j+1], y_n[i][j], a2);
    ny=n_y(x_n[i][j+1], x_n[i][j], a2);
    
    
    a_L=sqrt(gamma*abs(V_L[3]/V_L[0]));
    a_R=sqrt(gamma*abs(V_R[3]/V_R[0]));
    
    T_vel_L=V_L[1]*nx+V_L[2]*ny;
    T_vel_R=V_R[1]*nx+V_R[2]*ny;
    
    M_L= T_vel_L/a_L;
    
    M_R=T_vel_R/a_R;
    mach_knight=(V_L[1]+V_R[1])/(a_L+a_R);
    M_plus= 0.25*(M_L+1)*(M_L+1);
    M_minus= -0.25*(M_R-1)*(M_R-1);
    
    
    
    
    
    beta_L=-max(0,1-int(abs(M_L)));
    
    beta_R=-max(0,1-int(abs(M_R)));
    
    
    alpha_plus=0.5*(1+copysign(1,M_L));
    alpha_minus=0.5*(1-copysign(1,M_R));
    c_plus = alpha_plus*(1+beta_L)*M_L -beta_L*M_plus;
    c_minus = alpha_minus*(1+beta_R)*M_R -beta_R*M_minus;
    
    
    
    
    
    P_plus = M_plus*(2-M_L);
    P_minus = M_minus*(-2-M_R);
    
    D_plus=alpha_plus*(1+beta_L)-beta_L*P_plus;
    
    D_minus=alpha_minus*(1+beta_R)-beta_R*P_minus;
    
    
    //calculate F_P
    _F_P_kasi[0]=0;
    
    _F_P_kasi[1]=D_plus*nx*V_L[3]+D_minus*nx*V_R[3];
    _F_P_kasi[2]=D_plus*ny*V_L[3]+D_minus*ny*V_R[3];
    
    _F_P_kasi[3]=0;
    
    //calculate F_c
    
    _F_C_kasi[0]= V_L[0]*a_L*c_plus*1 + V_R[0]*a_R*c_minus*1;
    
    _F_C_kasi[1]= V_L[0]*a_L*c_plus*V_L[1] + V_R[0]*a_R*c_minus*V_R[1];
    _F_C_kasi[2]= V_L[0]*a_L*c_plus*V_L[2] + V_R[0]*a_R*c_minus*V_R[2];
    _F_C_kasi[3]= V_L[0]*a_L*c_plus*ht_L + V_R[0]*a_R*c_minus*ht_R;
    
    
    _F[0][i][j]=_F_C_kasi[0]+_F_P_kasi[0];
    _F[1][i][j]=_F_C_kasi[1]+_F_P_kasi[1];
    _F[2][i][j]=_F_C_kasi[2]+_F_P_kasi[2];
    _F[3][i][j]=_F_C_kasi[3]+_F_P_kasi[3];

   
    
}
                        


                
                


void Dv_nozzle::mesh_nodes()
{
    im_f = _imax_f-1;
    im_c = _imax_c-1;
 // double imax=65;
 // double jmax=65;
 
    std::ifstream infile("/Users/cringedaddy/CFD class/Project/Project_Files/Grids/curviliniear-grids/curv2d65.grd");

    int nzones, imax, jmax,kmax;
      infile >> nzones >> imax >> jmax>>kmax;


       

       // Read in x-coordinate
    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                infile >> x_n[i][j];
            }
        }
        
    }
       // Read in y-coordinate
    for (int k = 0; k < kmax; k++) {
        for (int j = 0; j < jmax; j++) {
            for (int i = 0; i < imax; i++) {
                infile >> y_n[i][j];
            }
        }
    }

       // Close the file
   infile.close();

   /*
    double x_min=0;
    double x_max=1;
    
  
    
    
    for (int j = 0; j < jmax; j++) {
        for (int i = 0; i < imax; i++) {
            x_n[i][j]= x_min + (float(i)/float(imax-1))*(x_max- x_min) ;
            y_n[i][j]= x_min + (float(j)/float(jmax-1))*(x_max- x_min) ;
       //     std::cout<<"point is "<<" ( "<<x_n[i][j]<<","<<y_n[i][j]<<" ) "<<std::endl;
          //  std::cout<<x_n[i][j]<<"             "<<y_n[i][j]<<std::endl;
         //   std::cout<<"y_n is "<<y_n[i][j]<<std::endl;
          //  std::cout<<"x_n is "<<x_n[i][j]<<std::endl;
        }
    }
    */
       // Do whatever you need to do with the grid data here

    for (int j = 0; j < jmax-1; j++) {
        for (int i = 0; i < imax-1; i++) {
            x_c[i][j]= 0.25*(x_n[i][j]+x_n[i+1][j]+x_n[i+1][j+1]+x_n[i][j+1]);
            y_c[i][j]= 0.25*(y_n[i][j]+y_n[i+1][j]+y_n[i+1][j+1]+y_n[i][j+1]);
          //  std::cout<<"y_c is "<<y_c[i][j]<<std::endl;
         //   std::cout<<"x_c is "<<x_c[i][j]<<std::endl;
           
        }
    }
 
    
    
    
       
   }
  




    
    
    



double Dv_nozzle:: face_area(double x_1,double x_2,double y_1,double y_2)


{
    
    double f_area;
    
    f_area= sqrt((x_1-x_2)*(x_1-x_2)+(y_1-y_2)*(y_1-y_2));
    
    return f_area;
    
    
}
double Dv_nozzle:: n_x(double y_1, double y_2,double f_area)

{
    double n_x;
    
    n_x = (y_1-y_2)/f_area;
    
    return n_x;
    
}
double Dv_nozzle:: n_y(double x_1, double x_2,double f_area)

{
    double n_y;
    
    n_y = -((x_1-x_2))/f_area;
    
    return n_y;
    
}
double Dv_nozzle:: cell_vol(int i, int j)


{
    double vol;
    
    vol=(x_n[i+1][j+1]-x_n[i][j])*(y_n[i+1][j]-y_n[i][j+1])-(y_n[i+1][j+1]-y_n[i][j])*(x_n[i+1][j]-x_n[i][j+1]);
    
    vol=0.5*abs(vol);
    
    
    
    
  //cout<<"vol is "<<vol<<endl;
    
    return vol;
    
}

double Dv_nozzle:: cell_vol_kasi(double area_1,double area_2)


{
    double vol;
    
    vol = 0.5*(area_1+area_2);
    
    return vol;
    
}
double Dv_nozzle:: cell_vol_eta(double area_1,double area_2)


{
    double vol;
    
    vol = 0.5*(area_1+area_2);
    
    return vol;
    
}






void Dv_nozzle::vanleer_eta (int i,int j)
{
    
    
    
    // flux for eta faces
         
         
    
         
         
         
    double a4=face_area(x_n[i+1][j], x_n[i][j],y_n[i+1][j], y_n[i][j]);
    
    
    nx=-n_x(y_n[i+1][j], y_n[i][j], a4);
    ny=-n_y(x_n[i+1][j], x_n[i][j], a4);
    
    
   
    
    
    a_D=sqrt(gamma*abs(V_D[3]/V_D[0]));
    a_U=sqrt(gamma*abs(V_U[3]/V_U[0]));
    
    T_vel_D=V_D[1]*nx+V_D[2]*ny;
    T_vel_U=V_U[1]*nx+V_U[2]*ny;
    
    M_D= T_vel_D/a_D;
    
    M_U=T_vel_U/a_U;
    mach_knight=(V_D[1]+V_U[1])/(a_D+a_U);
    M_plus= 0.25*(M_D+1)*(M_D+1);
    M_minus= -0.25*(M_U-1)*(M_U-1);
    
    
    
    
    
    beta_D=-max(0,1-int(abs(M_D)));
    
    beta_U=-max(0,1-int(abs(M_U)));
    
    
    alpha_plus=0.5*(1+copysign(1,M_D));
    alpha_minus=0.5*(1-copysign(1,M_U));
    c_plus = alpha_plus*(1+beta_D)*M_D -beta_D*M_plus;
    c_minus = alpha_minus*(1+beta_U)*M_U -beta_U*M_minus;
    
    
    
    
    
    P_plus = M_plus*(2-M_D);
    P_minus = M_minus*(-2-M_U);
    
    D_plus=alpha_plus*(1+beta_D)-beta_D*P_plus;
    
    D_minus=alpha_minus*(1+beta_U)-beta_U*P_minus;
    
    
    //calculate F_P
    _F_P_eta[0]=0;
    
    _F_P_eta[1]=D_plus*nx*V_D[3]+D_minus*nx*V_U[3];
    _F_P_eta[2]=D_plus*ny*V_D[3]+D_minus*ny*V_U[3];
    
    _F_P_eta[3]=0;
    
    //calculate F_c
    
    _F_C_eta[0]= V_D[0]*a_D*c_plus*1 + V_U[0]*a_U*c_minus*1;
    
    _F_C_eta[1]= V_D[0]*a_D*c_plus*V_D[1] + V_U[0]*a_U*c_minus*V_U[1];
    _F_C_eta[2]= V_D[0]*a_D*c_plus*V_D[2] + V_U[0]*a_U*c_minus*V_U[2];
    _F_C_eta[3]= V_D[0]*a_D*c_plus*ht_D + V_U[0]*a_U*c_minus*ht_U;
    
    
    _F_eta[0][i][j]=_F_C_eta[0]+_F_P_eta[0];
    _F_eta[1][i][j]=_F_C_eta[1]+_F_P_eta[1];
    _F_eta[2][i][j]=_F_C_eta[2]+_F_P_eta[2];
    _F_eta[3][i][j]=_F_C_eta[3]+_F_P_eta[3];

         
    
    
}
void Dv_nozzle::print_res(){
    rho_vs_x<<"######### x ######   "<<"######### y ######   "<<"#########rho######   "<<std::endl;
    u_vs_x<<"######### x ######   "<<"######### y ######   "<<"#########ux######   "<<std::endl;
    v_vs_x<<"######### x ######   "<<"######### y ######   "<<"#########vy######   "<<std::endl;
    p_vs_x<<"######### x ######   "<<"######### y ######   "<<"#########P######   "<<std::endl;
    
    for (int i= 0;i<=im_c;i++)
    {
        for(int j=0;j<=im_c;j++)
        {
            rho_vs_x<<std::setprecision(5)<<x_c[i][j]<<","<<std::setprecision(5)<<y_c[i][j]<<","<<std::setprecision(5)<< _V[0][i][j]<<std::endl;
            
            
            
            u_vs_x<<std::setprecision(5)<<x_c[i][j]<<","<<std::setprecision(5)<<y_c[i][j]<<","<<std::setprecision(14)<< _V[1][i][j]<<std::endl;
            
            v_vs_x<<std::setprecision(5)<<x_c[i][j]<<","<<std::setprecision(5)<<y_c[i][j]<<","<<std::setprecision(14)<< _V[2][i][j]<<std::endl;
            
            
            p_vs_x<<std::setprecision(5)<<x_c[i][j]<<","<<std::setprecision(5)<<y_c[i][j]<<","<<std::setprecision(5)<< _V[3][i][j]<<std::endl;
            
        }
        
    }
    
}
    
