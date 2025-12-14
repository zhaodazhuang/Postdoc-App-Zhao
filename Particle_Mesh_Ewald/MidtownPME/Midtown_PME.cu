#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<sstream>
#include<complex>
#include<algorithm> 
#include<cufftXt.h>
#include<cuComplex.h>
#include <time.h>
#define a_number 128

__device__ double round(double x){ return (x>0.0)? floor(x+0.5):ceil(x-0.5); }



void read_xyzfile(int& npar,
                  std::string& file_name,
                  std::vector<std::string>& element,
                  std::vector<double>& rX,
                  std::vector<double>& rY,
                  std::vector<double>& rZ ){

    std::ifstream coor_I(file_name, std::ios::in);
        if ( !coor_I.is_open() ){ 
            std::cout << file_name <<" is not exist"<< std::endl; 
            exit(0); }
        std::string line;

        getline(coor_I, line);
        getline(coor_I, line);

        for (int i = 0; i < npar; ++i){

            getline(coor_I, line);

            if ( coor_I.eof() ){
                std::cout << "end of file" << std::endl;
                exit(0); 
            }

            std::istringstream is(line);
            is >> element[i] >> rX[i] >> rY[i] >> rZ[i];
        }   

    coor_I.close();

    return;
}



void calc_unit_conversion_coeffs( double& r4pie ){

    double pi = acos(-1.e0);
    double ele_to_clbp19 = 1.602176634e0;
    double epslonp12 = 8.8541878128e0;
    double nan23=6.02214076e0;

    r4pie = pow(ele_to_clbp19,2.e0)*nan23/(4.e0*pi*epslonp12)*1.0e6;

    return;
}


void record_Midt_grid(int& n_order,
                      int& points_n,
                      std::vector<int>& Midt_grid ){

    if ( n_order == 4 ){

        Midt_grid[0] = 1;    Midt_grid[1] = 1;    Midt_grid[2] = 1;
        Midt_grid[3] = 2;    Midt_grid[4] = 1;    Midt_grid[5] = 1;
        Midt_grid[6] = 1;    Midt_grid[7] = 2;    Midt_grid[8] = 1;
        Midt_grid[9] = 1;    Midt_grid[10] = 1;   Midt_grid[11] = 2;
    }
    else if ( n_order == 6 ){

        Midt_grid[0] = 1;    Midt_grid[1] = 1;    Midt_grid[2] = 1;
        Midt_grid[3] = 2;    Midt_grid[4] = 1;    Midt_grid[5] = 1;
        Midt_grid[6] = 1;    Midt_grid[7] = 2;    Midt_grid[8] = 1;
        Midt_grid[9] = 1;    Midt_grid[10] = 1;   Midt_grid[11] = 2;
        Midt_grid[12] = 3;   Midt_grid[13] = 1;   Midt_grid[14] = 1;
        Midt_grid[15] = 1;   Midt_grid[16] = 3;   Midt_grid[17] = 1;
        Midt_grid[18] = 1;   Midt_grid[19] = 1;   Midt_grid[20] = 3;
        Midt_grid[21] = 2;   Midt_grid[22] = 2;   Midt_grid[23] = 1;
        Midt_grid[24] = 2;   Midt_grid[25] = 1;   Midt_grid[26] = 2;
        Midt_grid[27] = 1;   Midt_grid[28] = 2;   Midt_grid[29] = 2;
        Midt_grid[30] = 2;   Midt_grid[31] = 2;   Midt_grid[32] = 2;
    }

    int num = points_n*3;
    int numd8 = num/8;
    int numd4 = num/4;
    int numd2 = num/2;

    for ( int i = numd8; i < numd4; i=i+3 ){

        Midt_grid[i] = 1 - Midt_grid[i-numd8];
        Midt_grid[i+1] = Midt_grid[i+1-numd8];
        Midt_grid[i+2] = Midt_grid[i+2-numd8];
    }


    for ( int i = numd4; i < numd2; i=i+3 ){

        Midt_grid[i] = Midt_grid[i-numd4];
        Midt_grid[i+1] = 1 - Midt_grid[i+1-numd4] ;
        Midt_grid[i+2] = Midt_grid[i+2-numd4];
    }


   for ( int i = numd2; i < num; i=i+3 ){

        Midt_grid[i] = Midt_grid[i-numd2];
        Midt_grid[i+1] = Midt_grid[i+1-numd2] ;
        Midt_grid[i+2] = 1 - Midt_grid[i+2-numd2];
    }

    return;
}


__global__ void cuda_round(double* d_rX,
                           double* d_rY,
                           double* d_rZ,
                           double* d_uX,
                           double* d_uY,
                           double* d_uZ,
                           int npar,
                           double lX,
                           double lY,
                           double lZ,
                           double invlX,
                           double invlY,
                           double invlZ,
                           int Kmax1,
                           int Kmax2,
                           int Kmax3)
{

    int n = blockDim.x*blockIdx.x + threadIdx.x;
    double uu;

    if ( n < npar ){
        uu = d_rX[n] - round(d_rX[n]/lX)*lX; d_uX[n] = double(Kmax1)*(invlX*uu+0.5e0);
        uu = d_rY[n] - round(d_rY[n]/lY)*lY; d_uY[n] = double(Kmax2)*(invlY*uu+0.5e0);
        uu = d_rZ[n] - round(d_rZ[n]/lZ)*lZ; d_uZ[n] = double(Kmax3)*(invlZ*uu+0.5e0);
    }

}


__device__ double Midt_L(double theta,
                         double zeta){

    double L;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    L = -0.5e0*thetasq*theta + 0.5e0*thetasq
        -(9.e0*zetasq-2.e0)*theta/6.e0 + 0.5e0*zetasq;

    return L;
}


__device__ double Midt_dL(double theta,
                          double zeta){

    double dL;

    dL = -1.5e0*theta*theta + theta
         -(9.e0*zeta*zeta-2.e0)/6.e0;

    return dL;
}


__device__ double Midt_R(double theta,
                         double zeta){

    double R;

    R = theta*theta*theta/6.e0 + (3.e0*zeta*zeta-1.e0)*theta/6.e0;

    return R;
}



__device__ double Midt_dR(double theta,
                          double zeta){

    double dR;

    dR = 0.5e0*theta*theta + (3.e0*zeta*zeta-1.e0)/6.e0;

    return dR;

}



__device__ double Midt_L111(double theta,
                            double zeta){

    double L;
    double thetasq = theta*theta;
    double thetasqsq = thetasq*thetasq;
    double zetasq = zeta*zeta;
    double zetasqsq = zetasq*zetasq;
  
    L = thetasqsq*theta/12.e0 - thetasqsq/6.e0
      + (10.e0*zetasq-1.e0)/12.e0*thetasq*theta
      - (6.e0*zetasq-1.e0)/6.e0*thetasq
      + (5.e0*zetasqsq-zetasq)/4.e0*theta
      - (3.e0*zetasqsq-zetasq)/6.e0;
  
    return L;
}


__device__ double Midt_dL111(double theta,
                             double zeta){

    double dL;
    double thetasq = theta*theta;
    double thetasqsq = thetasq*thetasq;
    double zetasq = zeta*zeta;
    double zetasqsq = zetasq*zetasq;

    dL = 5.e0*thetasqsq/12.e0 - 2.e0*thetasq*theta/3.e0
       + (10.e0*zetasq-1.e0)/4.e0*thetasq
       - (6.e0*zetasq-1.e0)/3.e0*theta
       + (5.e0*zetasqsq-zetasq)/4.e0;

    return dL;
}


__device__ double Midt_L311(double theta,
                            double zeta){

    double L;
    double thetasq = theta*theta;
    double thetasqsq = thetasq*thetasq;
    double zetasq = zeta*zeta;
    double zetasqsq = zetasq*zetasq;
  
    L = thetasqsq*theta/120.e0
      + (2.e0*zetasq-1.e0)/24.e0*thetasq*theta
      + (15.e0*(zetasqsq-zetasq) + 4.e0)/120.e0*theta;
  
    return L;
}


__device__ double Midt_dL311(double theta,
                             double zeta){

    double dL;
    double thetasq = theta*theta;
    double thetasqsq = thetasq*thetasq;
    double zetasq = zeta*zeta;
    double zetasqsq = zetasq*zetasq;

    dL = thetasqsq/24.e0
       + (2.e0*zetasq-1.e0)/8.e0*thetasq
       + (15.e0*(zetasqsq-zetasq) + 4.e0)/120.e0;

    return dL;
}


__device__ double Midt_L211(double theta,
                            double zeta){

    double L;

    L = -0.25e0*Midt_L111(theta,zeta)
        - 2.5e0*Midt_L311(theta,zeta);

    return L;
}


__device__ double Midt_dL211(double theta,
                             double zeta){

    double dL;

    dL = -0.25e0*Midt_dL111(theta,zeta)
         - 2.5e0*Midt_dL311(theta,zeta);

    return dL;
}


__device__ double Midt_S(double theta,
                         double zeta){

    double S;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    S = 0.5e0*( -thetasq*theta + thetasq
      -(3.e0*zetasq-2.e0)*theta + zetasq );

    return S;
}


__device__ double Midt_dS(double theta,
                          double zeta){

    double dS;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    dS = 0.5e0*( -3.e0*thetasq + 2.e0*theta
       -(3.e0*zetasq-2.e0) );

    return dS;
}


__global__ void  cuda_calc_order4_coeffs(int npar,
                                         double *d_uX,
                                         double *d_uY,
                                         double *d_uZ,
                                         int n_order,
                                         int points_n,
                                         double zeta,
                                         int *d_Midt_grid,
                                         double *d_Midt_coeffs,
                                         double *d_Midt_diff1,
                                         double *d_Midt_diff2,
                                         double *d_Midt_diff3)
{

    int par_i = blockIdx.x; // partical Index
    int point_i = threadIdx.x; // grid point Index
    int point_coor = point_i*3; // grid point coor

    int n = par_i*points_n + point_i; // coeffs Index

    if ( n < points_n*npar ){

        double thetaX = d_uX[par_i] - int(d_uX[par_i]);
        double thetaY = d_uY[par_i] - int(d_uY[par_i]);
        double thetaZ = d_uZ[par_i] - int(d_uZ[par_i]);

        double theta1, theta2, theta3;
        double sgn1, sgn2, sgn3;


        if ( 0.5e0 - d_Midt_grid[point_coor] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
        else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }

        if ( 0.5e0 - d_Midt_grid[point_coor+1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
        else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }

        if ( 0.5e0 - d_Midt_grid[point_coor+2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
        else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }


        if ( point_i%4 == 0 ){

            d_Midt_coeffs[n] = Midt_L(theta1,zeta)*theta2*theta3
                             + Midt_L(theta2,zeta)*theta1*theta3
                             + Midt_L(theta3,zeta)*theta2*theta1;

            /////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dL(theta1,zeta)*theta2*theta3
                             + Midt_L(theta2,zeta)*theta3
                             + Midt_L(theta3,zeta)*theta2 );

            d_Midt_diff2[n] = sgn2*(  Midt_L(theta1,zeta)*theta3
                             + Midt_dL(theta2,zeta)*theta1*theta3
                             + Midt_L(theta3,zeta)*theta1 );

            d_Midt_diff3[n] = sgn3*( Midt_L(theta1,zeta)*theta2
                             + Midt_L(theta2,zeta)*theta1
                             + Midt_dL(theta3,zeta)*theta1*theta2 );
        }
        else if ( point_i%4 == 1 ){

            d_Midt_coeffs[n] = Midt_R(theta1,zeta)*theta2*theta3;
            //////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*Midt_dR(theta1,zeta)*theta2*theta3;
            d_Midt_diff2[n] = sgn2*Midt_R(theta1,zeta)*theta3;
            d_Midt_diff3[n] = sgn3*Midt_R(theta1,zeta)*theta2;
        }
        else if ( point_i%4 == 2 ){

            d_Midt_coeffs[n] = Midt_R(theta2,zeta)*theta1*theta3;
            ///////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*Midt_R(theta2,zeta)*theta3;
            d_Midt_diff2[n] = sgn2*Midt_dR(theta2,zeta)*theta1*theta3;
            d_Midt_diff3[n] = sgn3*Midt_R(theta2,zeta)*theta1;
        }
        else if ( point_i%4 == 3 ){

            d_Midt_coeffs[n] = Midt_R(theta3,zeta)*theta2*theta1;
            ///////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*Midt_R(theta3,zeta)*theta2;
            d_Midt_diff2[n] = sgn2*Midt_R(theta3,zeta)*theta1;
            d_Midt_diff3[n] = sgn3*Midt_dR(theta3,zeta)*theta1*theta2;
        }

    }


}


__global__ void  cuda_calc_order6_coeffs(int npar,
                                         double *d_uX,
                                         double *d_uY,
                                         double *d_uZ,
                                         int n_order,
                                         int points_n,
                                         double zeta,
                                         int *d_Midt_grid,
                                         double *d_Midt_coeffs,
                                         double *d_Midt_diff1,
                                         double *d_Midt_diff2,
                                         double *d_Midt_diff3)
{

    int par_i = blockIdx.x; // partical Index
    int point_i = threadIdx.x;// grid point Index
    int point_coor = point_i*3; // grid point coor

    int n = par_i*points_n + point_i;

    if ( n < points_n*npar ){

        double thetaX = d_uX[par_i] - int(d_uX[par_i]);
        double thetaY = d_uY[par_i] - int(d_uY[par_i]);
        double thetaZ = d_uZ[par_i] - int(d_uZ[par_i]);

        double theta1, theta2, theta3;
        double sgn1, sgn2, sgn3;

        if ( 0.5e0 - d_Midt_grid[point_coor] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
        else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }

        if ( 0.5e0 - d_Midt_grid[point_coor+1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
        else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }

        if ( 0.5e0 - d_Midt_grid[point_coor+2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
        else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }


        if ( point_i%11 == 0 ){

            d_Midt_coeffs[n] = Midt_L111(theta1,zeta)*theta2*theta3
                             + Midt_L111(theta2,zeta)*theta1*theta3
                             + Midt_L111(theta3,zeta)*theta1*theta2
                             + Midt_S(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta);
            ////////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dL111(theta1,zeta)*theta2*theta3
                             + Midt_L111(theta2,zeta)*theta3
                             + Midt_L111(theta3,zeta)*theta2
            + Midt_dS(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff2[n] = sgn2*( Midt_L111(theta1,zeta)*theta3
                             + Midt_dL111(theta2,zeta)*theta1*theta3
                             + Midt_L111(theta3,zeta)*theta1
            + Midt_S(theta1,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff3[n] = sgn3*( Midt_L111(theta1,zeta)*theta2
                             + Midt_L111(theta2,zeta)*theta1
                             + Midt_dL111(theta3,zeta)*theta1*theta2
            + Midt_S(theta1,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta3,zeta) );
        }
        else if ( point_i%11 == 1 ){

            d_Midt_coeffs[n] = Midt_L211(theta1,zeta)*theta2*theta3
                             + Midt_R(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta);
            //////////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dL211(theta1,zeta)*theta2*theta3
            + Midt_dR(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff2[n] = sgn2*( Midt_L211(theta1,zeta)*theta3
            + Midt_R(theta1,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff3[n] = sgn3*( Midt_L211(theta1,zeta)*theta2
            + Midt_R(theta1,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta3,zeta) );
        }
        else if ( point_i%11 == 2 ){

            d_Midt_coeffs[n] = Midt_L211(theta2,zeta)*theta1*theta3
                             + Midt_R(theta2,zeta)*Midt_S(theta1,zeta)*Midt_S(theta3,zeta);
            ///////////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_L211(theta2,zeta)*theta3
            + Midt_R(theta2,zeta)*Midt_dS(theta1,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff2[n] = sgn2*( Midt_dL211(theta2,zeta)*theta1*theta3
            + Midt_dR(theta2,zeta)*Midt_S(theta1,zeta)*Midt_S(theta3,zeta) );

            d_Midt_diff3[n] = sgn3*( Midt_L211(theta2,zeta)*theta1
            + Midt_R(theta2,zeta)*Midt_S(theta1,zeta)*Midt_dS(theta3,zeta) );
        }
        else if ( point_i%11 == 3 ){

            d_Midt_coeffs[n] = Midt_L211(theta3,zeta)*theta2*theta1
                             + Midt_R(theta3,zeta)*Midt_S(theta2,zeta)*Midt_S(theta1,zeta);

            //////////////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_L211(theta3,zeta)*theta2
            + Midt_R(theta3,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta1,zeta) );

            d_Midt_diff2[n] = sgn2*( Midt_L211(theta3,zeta)*theta1
            + Midt_R(theta3,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta1,zeta) );

            d_Midt_diff3[n] = sgn3*( Midt_dL211(theta3,zeta)*theta1*theta2
            + Midt_dR(theta3,zeta)*Midt_S(theta2,zeta)*Midt_S(theta1,zeta) );
        }
        else if ( point_i%11 == 4 ){

            d_Midt_coeffs[n] = Midt_L311(theta1,zeta)*theta2*theta3;
            ///////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dL311(theta1,zeta)*theta2*theta3 );
            d_Midt_diff2[n] = sgn2*( Midt_L311(theta1,zeta)*theta3 );
            d_Midt_diff3[n] = sgn3*( Midt_L311(theta1,zeta)*theta2 );
        }
        else if ( point_i%11 == 5 ){

            d_Midt_coeffs[n] = Midt_L311(theta2,zeta)*theta1*theta3;
            ///////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_L311(theta2,zeta)*theta3 );
            d_Midt_diff2[n] = sgn2*( Midt_dL311(theta2,zeta)*theta1*theta3 );
            d_Midt_diff3[n] = sgn3*( Midt_L311(theta2,zeta)*theta1 );
        }
        else if ( point_i%11 == 6 ){

            d_Midt_coeffs[n] = Midt_L311(theta3,zeta)*theta1*theta2;
            //////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_L311(theta3,zeta)*theta2 );
            d_Midt_diff2[n] = sgn2*( Midt_L311(theta3,zeta)*theta1 );
            d_Midt_diff3[n] = sgn3*( Midt_dL311(theta3,zeta)*theta1*theta2 );
        }
        else if ( point_i%11 == 7 ){

            d_Midt_coeffs[n] = Midt_R(theta1,zeta)*Midt_R(theta2,zeta)*Midt_S(theta3,zeta);
            ////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dR(theta1,zeta)*
                               Midt_R(theta2,zeta)*Midt_S(theta3,zeta) );
            d_Midt_diff2[n] = sgn2*( Midt_R(theta1,zeta)*
                               Midt_dR(theta2,zeta)*Midt_S(theta3,zeta) );
            d_Midt_diff3[n] = sgn3*( Midt_R(theta1,zeta)*
                               Midt_R(theta2,zeta)*Midt_dS(theta3,zeta) );
        }
        else if ( point_i%11 == 8 ){

            d_Midt_coeffs[n] = Midt_R(theta1,zeta)*Midt_R(theta3,zeta)*Midt_S(theta2,zeta);
            ////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dR(theta1,zeta)*
                                 Midt_R(theta3,zeta)*Midt_S(theta2,zeta) );
            d_Midt_diff2[n] = sgn2*( Midt_R(theta1,zeta)*
                                 Midt_R(theta3,zeta)*Midt_dS(theta2,zeta) );
            d_Midt_diff3[n] = sgn3*( Midt_R(theta1,zeta)*
                                 Midt_dR(theta3,zeta)*Midt_S(theta2,zeta) );
        }
        else if ( point_i%11 == 9 ){

            d_Midt_coeffs[n] = Midt_R(theta3,zeta)*Midt_R(theta2,zeta)*Midt_S(theta1,zeta);
            //////////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_R(theta3,zeta)*
                                  Midt_R(theta2,zeta)*Midt_dS(theta1,zeta) );
            d_Midt_diff2[n] = sgn2*( Midt_R(theta3,zeta)*
                                  Midt_dR(theta2,zeta)*Midt_S(theta1,zeta) );
            d_Midt_diff3[n] = sgn3*( Midt_dR(theta3,zeta)*
                                  Midt_R(theta2,zeta)*Midt_S(theta1,zeta) );
        }
        else if ( point_i%11 == 10 ){

            d_Midt_coeffs[n] = Midt_R(theta1,zeta)*Midt_R(theta2,zeta)*Midt_R(theta3,zeta);
            //////////////////////////////////////////////////////////
            d_Midt_diff1[n] = sgn1*( Midt_dR(theta1,zeta)*
                               Midt_R(theta2,zeta)*Midt_R(theta3,zeta) );
            d_Midt_diff2[n] = sgn2*( Midt_R(theta1,zeta)*
                               Midt_dR(theta2,zeta)*Midt_R(theta3,zeta) );
            d_Midt_diff3[n] = sgn3*( Midt_R(theta1,zeta)*
                               Midt_R(theta2,zeta)*Midt_dR(theta3,zeta) );
        }
    }
}


__global__ void  cuda_assign_charge_grid(int npar,
                                         double *d_uX,
                                         double *d_uY,
                                         double *d_uZ,
                                         double *d_q,
                                         int Kmax1,
                                         int Kmax2,
                                         int Kmax3,
                                         int points_n,
                                         int *d_Midt_grid,
                                         double *d_Midt_coeffs,
                                         double *d_charge_grid)
{

    int par_i = blockIdx.x; // partical Index
    int point_i = threadIdx.x; // grid point Index
    int point_coor = point_i*3; // grid point coor
    int coeffs_i = par_i*points_n + point_i;

    // (k1, k2, k3 ) charge point coor
    int k1 = int(d_uX[par_i])+d_Midt_grid[point_coor];
    if ( k1 < 0 ){ k1 = k1 + Kmax1; }
    if ( k1 >= Kmax1 ){ k1 = k1 - Kmax1; }

    int k2 = int(d_uY[par_i])+d_Midt_grid[point_coor+1];
    if ( k2 < 0 ){ k2 = k2 + Kmax2; }
    if ( k2 >= Kmax2 ){ k2 = k2 - Kmax2; }

    int k3 = int(d_uZ[par_i])+d_Midt_grid[point_coor+2];
    if ( k3 < 0 ){ k3 = k3 + Kmax3; }
    if ( k3 >= Kmax3 ){ k3 = k3 - Kmax3; }


    int n = k1 + k2*Kmax1 + k3*Kmax1*Kmax2; // charge point Index

    if ( n < Kmax1*Kmax2*Kmax3 ){

        double val = d_q[par_i]*d_Midt_coeffs[coeffs_i];
        atomicAdd(&d_charge_grid[n], val);
    }

}

__global__ void cuda_initialization_ZZU(int Kmax1,
                                        int Kmax2,
                                        int Kmax3,
                                        double *d_charge_grid,
                                        cuDoubleComplex *d_ZZU)
{
    int m1 = blockIdx.x*blockDim.x + threadIdx.x;
    int m2 = blockIdx.y;
    int m3 = blockIdx.z;

    if ( m1 < Kmax1 && m2 < Kmax2 && m3 < Kmax3 ){

       int n = m1 + m2*Kmax1 + m3*Kmax1*Kmax2;

       d_ZZU[n] = make_cuDoubleComplex(d_charge_grid[n],0);
    }

}


__global__ void cuda_add_coeffs_ZZU(double r4alpsq,
                                    double tpi,
                                    double invlX,
                                    double invlY, 
                                    double invlZ,
                                    int Kmax1,
                                    int Kmax2,
                                    int Kmax3,
                                    double rec_cutsq,
                                    double zetasq,
                                    cuDoubleComplex *d_ZZU)
{

    int m1 = blockIdx.x*blockDim.x + threadIdx.x;
    int m2 = blockIdx.y;
    int m3 = blockIdx.z;

    if ( m1 < Kmax1 && m2 < Kmax2 && m3 < Kmax3 ){

        int n = m1 + m2*Kmax1 + m3*Kmax1*Kmax2;
 
        if ( m1 > Kmax1/2 ){ m1 = m1-Kmax1; };
        if ( m2 > Kmax2/2 ){ m2 = m2-Kmax2; };
        if ( m3 > Kmax3/2 ){ m3 = m3-Kmax3; };    

        double rkX=tpi*double(m1)*invlX;
        double rkY=tpi*double(m2)*invlY;
        double rkZ=tpi*double(m3)*invlZ;

        double rksq = rkX*rkX + rkY*rkY + rkZ*rkZ;
               
        if ( rksq > 1.e-6 && rksq <= rec_cutsq ){

            double psi = pow(double(m1)/double(Kmax1),2.e0)
                       + pow(double(m2)/double(Kmax2),2.e0)
                       + pow(double(m3)/double(Kmax3),2.e0);

            d_ZZU[n] = cuCmul(d_ZZU[n], make_cuDoubleComplex(
                       exp(tpi*tpi*zetasq*psi)*exp(r4alpsq*rksq)/rksq, 0));
        }
        else{

            d_ZZU[n] = make_cuDoubleComplex(0,0);
        } 
        
    }  

}


__global__ void product_kernel(int Kmax1,
                           int Kmax2,
                           int Kmax3,
                           cuDoubleComplex *d_E_complex,
                           double *d_charge_grid,
                           cuDoubleComplex *d_ZZU)
{

    int m1 = blockIdx.x*blockDim.x + threadIdx.x;
    int m2 = blockIdx.y;
    int m3 = blockIdx.z;

    if ( m1 < Kmax1 && m2 < Kmax2 && m3 < Kmax3 ){

       int n = m1 + m2*Kmax1 + m3*Kmax1*Kmax2;

       d_E_complex[n] = cuCmul(d_ZZU[n], make_cuDoubleComplex(d_charge_grid[n],0));
 
    }

}



__global__ void reduction_complex_sum(cuDoubleComplex *X,
                                      size_t input_size)
{

    __shared__ cuDoubleComplex partialSum[2 * a_number];

    int i = 2 * blockIdx.x * blockDim.x + threadIdx.x;

    if(i < input_size) partialSum[threadIdx.x] = X[i];
    else partialSum[threadIdx.x] = make_cuDoubleComplex(0,0);

    if(i + blockDim.x < input_size) partialSum[threadIdx.x + blockDim.x] = X[i + blockDim.x];
    else partialSum[threadIdx.x + blockDim.x] = make_cuDoubleComplex(0,0);

    __syncthreads();


    int t = threadIdx.x;
    for(int stride = blockDim.x; stride >= 1; stride /= 2){
        if(t < stride)
            partialSum[t] = cuCadd(partialSum[t], partialSum[t + stride]);
            __syncthreads();
    }

    if(t == 0){
        X[blockIdx.x] = partialSum[t];
    }

}



__global__  void  cuda_calc_Force(int npar,
                                  double *d_uX,
                                  double *d_uY,
                                  double *d_uZ,
                                  double *d_FX,
                                  double *d_FY,
                                  double *d_FZ,
                                  int Kmax1,
                                  int Kmax2,
                                  int Kmax3,
                                  int points_n,
                                  int *d_Midt_grid,
                                  cuDoubleComplex *d_ZZU,
                                  double *d_Midt_diff1,
                                  double *d_Midt_diff2,
                                  double *d_Midt_diff3)
{

    int par_i = blockIdx.x; // partical Index
    int point_i = threadIdx.x; // grid point Index
    int point_coor = point_i*3; // grid point coor
    int coeffs_i = par_i*points_n + point_i;
    
    // (k1, k2, k3 ) charge point coor
    int k1 = int(d_uX[par_i])+d_Midt_grid[point_coor];
    if ( k1 < 0 ){ k1 = k1 + Kmax1; }
    if ( k1 >= Kmax1 ){ k1 = k1 - Kmax1; }

    int k2 = int(d_uY[par_i])+d_Midt_grid[point_coor+1];
    if ( k2 < 0 ){ k2 = k2 + Kmax2; }
    if ( k2 >= Kmax2 ){ k2 = k2 - Kmax2; }

    int k3 = int(d_uZ[par_i])+d_Midt_grid[point_coor+2];
    if ( k3 < 0 ){ k3 = k3 + Kmax3; }
    if ( k3 >= Kmax3 ){ k3 = k3 - Kmax3; }

    int n = k1 + k2*Kmax1 + k3*Kmax1*Kmax2; // charge point Index

    if ( n < Kmax1*Kmax2*Kmax3 ){

        double ZZUR = cuCreal(d_ZZU[n]);
        atomicAdd(&d_FX[par_i], ZZUR*d_Midt_diff1[coeffs_i]);
        atomicAdd(&d_FY[par_i], ZZUR*d_Midt_diff2[coeffs_i]);
        atomicAdd(&d_FZ[par_i], ZZUR*d_Midt_diff3[coeffs_i]);
    }

}


__global__  void  cuda_multi_coeffs_for_force(int npar,
                                              double invlX,
                                              double invlY,
                                              double invlZ,
                                              double *d_q,
                                              int Kmax1,
                                              int Kmax2,
                                              int Kmax3,
                                              double F_coeffs,
                                              double *d_FX, 
                                              double *d_FY, 
                                              double *d_FZ)
{

    int par_i = blockIdx.x*blockDim.x + threadIdx.x;

    if ( par_i < npar ){

        d_FX[par_i] = F_coeffs*d_q[par_i]*d_FX[par_i]*double(Kmax1)*invlX;
        d_FY[par_i] = F_coeffs*d_q[par_i]*d_FY[par_i]*double(Kmax2)*invlY;
        d_FZ[par_i] = F_coeffs*d_q[par_i]*d_FZ[par_i]*double(Kmax3)*invlZ;
    }

}


__global__ void reduction_sum(double *X, 
                              size_t input_size)
{

    __shared__ double partialSum[2 * a_number];

    int i = 2 * blockIdx.x * blockDim.x + threadIdx.x;

    if(i < input_size) partialSum[threadIdx.x] = X[i];
    else partialSum[threadIdx.x] = 0.0;

    if(i + blockDim.x < input_size) partialSum[threadIdx.x + blockDim.x] = X[i + blockDim.x];
    else partialSum[threadIdx.x + blockDim.x] = 0.0;

    __syncthreads();


    int t = threadIdx.x;
    for(int stride = blockDim.x; stride >= 1; stride /= 2){
        if(t < stride)
            partialSum[t] += partialSum[t + stride];
            __syncthreads();
    }

    if(t == 0){
        X[blockIdx.x] = partialSum[t];
    }

}


__global__ void  Force_sum_to_0(int npar,
                                double *d_FX, 
                                double *d_FY,
                                double *d_FZ,
                                double FXToT,
                                double FYToT, 
                                double FZToT)
{

    int par_i = blockIdx.x*blockDim.x + threadIdx.x;

    if ( par_i < npar ){

        d_FX[par_i] = d_FX[par_i] - FXToT;
        d_FY[par_i] = d_FY[par_i] - FYToT;
        d_FZ[par_i] = d_FZ[par_i] - FZToT;  
    }

}


void calc_recSpace_Energy_Force(int& npar,
                                double& lX,
                                double& lY,
                                double& lZ,
                                double& r4pie,
                                std::vector<double>& q,
                                std::vector<double>& rX,
                                std::vector<double>& rY,
                                std::vector<double>& rZ,
                                double& alpha,
                                int& n_order,
                                int& Kmax1,
                                int& Kmax2,
                                int& Kmax3,
                                cufftHandle& plan,
                                int& points_n,
                                double& zeta,
                                std::vector<int>& Midt_grid,
                                double& E,
                                std::vector<double>& FX,
                                std::vector<double>& FY,
                                std::vector<double>& FZ)
{

    /*----------------------------------------------
    *
    *  Calculate energy and forces in reciprocal space 
    *  through smooth partical mesh Ewald (SPME)
    *
    *  Read this article :
    *  " A smooth particle mesh Ewald method "
    *-------------------------------------------------
    */ 
 
    double invlX = 1.e0/lX;
    double invlY = 1.e0/lY;
    double invlZ = 1.e0/lZ;

    double tpi = 2.e0*acos(-1.e0); 
  
    /*---------------------------------------------
    *  STEP1 : Obtain scaled fraction coordinates
    *  
    *      uX = rX/lX*Kamx1     see page 8579
    *----------------------------------------------
    */
   
    int M_r = npar*sizeof(double);
    double *d_rX, *d_rY, *d_rZ;
    double *d_uX, *d_uY, *d_uZ;

    cudaMalloc(&d_rX, M_r);
    cudaMalloc(&d_rY, M_r);
    cudaMalloc(&d_rZ, M_r);
    cudaMalloc(&d_uX, M_r);
    cudaMalloc(&d_uY, M_r);
    cudaMalloc(&d_uZ, M_r);

    cudaMemcpy(d_rX, rX.data(), M_r, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rY, rY.data(), M_r, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rZ, rZ.data(), M_r, cudaMemcpyHostToDevice);


    int N_block = a_number;
    int N_grid = (npar+N_block-1)/N_block;

    cuda_round<<<N_grid,N_block>>>(d_rX, d_rY, d_rZ,
                                   d_uX, d_uY, d_uZ,
                                   npar,
                                   lX, lY, lZ,
                                   invlX, invlY, invlZ,
                                   Kmax1, Kmax2, Kmax3);
      
    /*-------------------------------------------------------------------
    * STEP 2: Obtain Midtown-splines interpolation coefficient 
    *        
    *
    *         Obtain charge_grid(k1,k2,k3) in real space see Eq.(4.6)
    *----------------------------------------------------------------
    */ 


    int *d_Midt_grid;
    double *d_Midt_coeffs;
    double *d_Midt_diff1;
    double *d_Midt_diff2;
    double *d_Midt_diff3;


    int M_Midt_coeffs = points_n*npar*sizeof(double);
    int M_Midt_grid = points_n*3*sizeof(int);

    cudaMalloc(&d_Midt_grid, M_Midt_grid);
    cudaMalloc(&d_Midt_coeffs, M_Midt_coeffs);
    cudaMalloc(&d_Midt_diff1, M_Midt_coeffs);
    cudaMalloc(&d_Midt_diff2, M_Midt_coeffs);
    cudaMalloc(&d_Midt_diff3, M_Midt_coeffs);


    cudaMemcpy(d_Midt_grid, Midt_grid.data(), M_Midt_grid, cudaMemcpyHostToDevice);


    N_block = points_n;
    N_grid = (points_n*npar+N_block-1)/N_block;

    if ( n_order == 4 ){

        cuda_calc_order4_coeffs<<<N_grid, N_block>>>(npar,
                                                     d_uX, d_uY, d_uZ,
                                                     n_order,
                                                     points_n,
                                                     zeta,
                                                     d_Midt_grid,
                                                     d_Midt_coeffs,
                                                     d_Midt_diff1,
                                                     d_Midt_diff2,
                                                     d_Midt_diff3);
    }
    else if ( n_order == 6 ){

        cuda_calc_order6_coeffs<<<N_grid, N_block>>>(npar,
                                                     d_uX, d_uY, d_uZ,
                                                     n_order,
                                                     points_n,
                                                     zeta,
                                                     d_Midt_grid,
                                                     d_Midt_coeffs,
                                                     d_Midt_diff1,
                                                     d_Midt_diff2,
                                                     d_Midt_diff3);
    }
 

    double *d_q;
    double *d_charge_grid;

    int M_q = npar*sizeof(double);
    int M_charge_grid = Kmax1*Kmax2*Kmax3*sizeof(double);

    cudaMalloc(&d_q, M_q);
    cudaMalloc(&d_charge_grid, M_charge_grid );

    cudaMemcpy(d_q, q.data(), M_q, cudaMemcpyHostToDevice);

    cuda_assign_charge_grid<<<N_grid, N_block>>>(npar,
                                                 d_uX, d_uY, d_uZ,
                                                 d_q,
                                                 Kmax1, Kmax2, Kmax3,
                                                 points_n,
                                                 d_Midt_grid,
                                                 d_Midt_coeffs,
                                                 d_charge_grid);

    cudaFree(d_Midt_grid);

    /*--------------------------------------------------------
    * STEP 3: Obtain F[Q](-m1, -m2, -m3) form FFT
    *   
    *   F[Q](-m1, -m2, -m3) =
    *
    *   \sum_{k1=0}^{Kmax1-1} \sum_{k2=0}^{Kmax2-1}
    *   \sum_{k3=0}^{Kmax3-1} Q(k1,k2,k3)*exp(-i2\pi m1*k1/Kmax1)*
    *                                     exp(-i2\pi m2*k2/Kmax2)*
    *                                     exp(-i2\pi m1*k3/Kmax3)
    *
    *   now in the code ZZU  <==> F[charge_grid](-m1, -m2, -m3) 
    *---------------------------------------------------------
    */ 


    cufftDoubleComplex *d_ZZU;
    int M_ZZU = Kmax1*Kmax2*Kmax3*sizeof(std::complex<double>);

    cudaMalloc(&d_ZZU, M_ZZU);

    
    N_block = a_number;
 
    dim3 grid_size((Kmax1+N_block-1)/N_block, Kmax2, Kmax3);

    cuda_initialization_ZZU<<<grid_size, N_block>>>(Kmax1, Kmax2, Kmax3,
                                                    d_charge_grid,
                                                    d_ZZU);
 
    cufftExecZ2Z(plan, d_ZZU, d_ZZU, CUFFT_FORWARD);


    /*-------------------------------------------------------
    *  STEP 4: Obtain 
    *    
    *    exp(-\pi^2 m^2/ beta^2)  
    *    ---------------------- * psi*F[Q](-m1, -m2, -m3)
    *           m^2
    *
    *   psi = exp(4\pi^2 \zeta^2 [(p1/Kmax1)^2 + (p2/Kmax2)^2 + (p3/Kmax3)^3])
    *  
    *   m = p/L
    *
    *   now in the code ZZU <==>
    *    
    *    exp(-\pi^2 m^2/ beta^2)
    *   ---------------------- * psi*F[Q](-m1, -m2, -m3)
    *           m^2
    *
    *------------------------------------------------------
    */ 
  
  
    double rec_cut, rec_cutsq;
    rec_cut = std::min( {double(Kmax1)*invlX, double(Kmax2)*invlY, double(Kmax3)*invlZ} );
    rec_cut = 0.5e0*rec_cut; rec_cut = 1.05e0*rec_cut*tpi;
    rec_cutsq = rec_cut*rec_cut; 
  
    double r4alpsq = -0.25e0/(alpha*alpha);
    double zetasq = zeta*zeta; 
 
    /*-------------------------------------------------------------
    *  K space needs to be centered on 0 to maintain central symmetry, 
    *  so m needs to be translated (m = m-Kmax/2),
    *  although Kmax is an even number resulting in asymmetry, 
    *  but there is rec_cut to make it ultimately centrally symmetric
    *------------------------------------------------------------------
    */
 
    cuda_add_coeffs_ZZU<<<grid_size,N_block>>>(r4alpsq,
                                               tpi,
                                               invlX, invlY, invlZ,
                                               Kmax1, Kmax2, Kmax3,
                                               rec_cutsq,
                                               zetasq,
                                               d_ZZU);

    /*---------------------------------------------------------
    *
    *  STEP 5): Complete transformation  F[f]*g = f*F[g]  see Eq.(B3)
    *
    *  \sum_{m1,m2,m3} F[Q](m1,m2,m3)*
    *
    *    exp(-\pi^2 m^2/ beta^2)
    * {  ---------------------- * psi*F[Q](-m1, -m2, -m3) }
    *            m^2
    *
    *  =
    *
    *  \sum_{m1,m2,m3} Q(k1,k2,k3)*
    *
    *      exp(-\pi^2 m^2/ beta^2)
    *  F[ ------------------------ * psi*F[Q](-m1, -m2, -m3) ]
    *              m^2
    * 
    *  
    *  see Eq.(4.7)
    *----------------------------------------------------------
    */ 
  

    cufftExecZ2Z(plan, d_ZZU, d_ZZU, CUFFT_INVERSE);
 
    cuDoubleComplex *d_E_complex;
    int M_E_complex = Kmax1*Kmax2*Kmax3*sizeof(std::complex<double>);

    cudaMalloc(&d_E_complex, M_E_complex);

    product_kernel<<<grid_size,N_block>>>(Kmax1,
                                          Kmax2,
                                          Kmax3,
                                          d_E_complex,
                                          d_charge_grid,
                                          d_ZZU);

    cudaFree(d_charge_grid);

    int len = Kmax1*Kmax2*Kmax3;
    N_block = a_number;
    N_grid = (len - 1)/N_block/2 + 1;
    while(N_grid >= 1){
        reduction_complex_sum<<<N_grid, N_block>>>(d_E_complex, len);
        if(N_grid == 1)
            break;
        len = N_grid;
        N_grid = (N_grid - 1)/N_block/2 + 1;
    }
   
    std::complex<double> E_complex;
    cudaMemcpy(&E_complex, d_E_complex, sizeof(std::complex<double>), cudaMemcpyDeviceToHost);
    cudaFree(d_E_complex);


    double E_coeffs = tpi*invlX*invlY*invlZ*r4pie;
    E = E_complex.real()*E_coeffs;
  
 
    /*----------------------------------------------------
    *
    *  STEP 6): Calculate force
    *
    *  \partial Q(k1,k2,k3) \partial M_n(u_i)   d u_i
    *  ------------------- = ---------------- * ------
    *   \partial r_i          \partial u_i      d r_i
    *
    *   u_i = Kmax/L*r_i
    *-------------------------------------------------------
    */

    double *d_FX, *d_FY, *d_FZ;
    int M_F = npar*sizeof(double);

    cudaMalloc(&d_FX, M_F);
    cudaMalloc(&d_FY, M_F);
    cudaMalloc(&d_FZ, M_F);
    
    cudaMemset(d_FX, 0, M_F);
    cudaMemset(d_FY, 0, M_F);
    cudaMemset(d_FZ, 0, M_F);

    N_block = points_n;
    N_grid = (points_n*npar+N_block-1)/N_block;

    cuda_calc_Force<<<N_grid, N_block>>>(npar,
                                         d_uX, d_uY, d_uZ,
                                         d_FX, d_FY, d_FZ,
                                         Kmax1, Kmax2, Kmax3,
                                         points_n,
                                         d_Midt_grid,
                                         d_ZZU,
                                         d_Midt_diff1,
                                         d_Midt_diff2,
                                         d_Midt_diff3);   

    cudaFree(d_Midt_grid);
    cudaFree(d_ZZU);
    cudaFree(d_Midt_diff1);
    cudaFree(d_Midt_diff2);
    cudaFree(d_Midt_diff3);

    double F_coeffs = -2.e0*E_coeffs;

    N_block = a_number;
    N_grid = (npar + N_block -1)/N_block;

    cuda_multi_coeffs_for_force<<<N_grid, N_block>>>(npar,
                                                     invlX, invlY, invlZ,
                                                     d_q,
                                                     Kmax1, Kmax2, Kmax3,
                                                     F_coeffs,
                                                     d_FX, d_FY, d_FZ); 

    double *dd_FX, *dd_FY, *dd_FZ;

    cudaMalloc(&dd_FX, M_F);
    cudaMalloc(&dd_FY, M_F);
    cudaMalloc(&dd_FZ, M_F);
 
    cudaMemcpy(dd_FX, d_FX, M_F, cudaMemcpyDeviceToDevice);
    cudaMemcpy(dd_FY, d_FY, M_F, cudaMemcpyDeviceToDevice);
    cudaMemcpy(dd_FZ, d_FZ, M_F, cudaMemcpyDeviceToDevice);


    len = npar;
    N_block = a_number;
    N_grid = (len - 1)/N_block/2 + 1;
    while(N_grid >= 1){
        reduction_sum<<<N_grid, N_block>>>(dd_FX, len);
        reduction_sum<<<N_grid, N_block>>>(dd_FY, len);
        reduction_sum<<<N_grid, N_block>>>(dd_FZ, len);
        if(N_grid == 1)
            break;
        len = N_grid;
        N_grid = (N_grid - 1)/N_block/2 + 1;
    }


    double FXToT, FYToT, FZToT;
 
    cudaMemcpy(&FXToT, dd_FX, sizeof(double), cudaMemcpyDeviceToHost);    
    cudaMemcpy(&FYToT, dd_FY, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(&FZToT, dd_FZ, sizeof(double), cudaMemcpyDeviceToHost);

    FXToT = FXToT/double(npar);
    FYToT = FYToT/double(npar);
    FZToT = FZToT/double(npar);

    N_block = a_number;
    N_grid = (npar+N_block-1)/N_block;

    // Force the sum of forces to be 0

    Force_sum_to_0<<<N_grid, N_block>>>(npar,
                                        d_FX, d_FY, d_FZ,
                                        FXToT, FYToT, FZToT);

    cudaMemcpy(FX.data(), d_FX, M_F, cudaMemcpyDeviceToHost);
    cudaMemcpy(FY.data(), d_FY, M_F, cudaMemcpyDeviceToHost);
    cudaMemcpy(FZ.data(), d_FZ, M_F, cudaMemcpyDeviceToHost);

    return;

}



int main(){


    // number of partical
    int npar = 12288;

    // Box size
    double lX = 49.6685068176563650e0;
    double lY = 49.6685068176563650e0;
    double lZ = 49.6685068176563650e0;

    // Element attributes and coordinates of particles    
    std::vector<std::string> element(npar);
    std::vector<double> rX(npar);
    std::vector<double> rY(npar);
    std::vector<double> rZ(npar);

 
    std::string file_name = "H2O.xyz";
    read_xyzfile( npar, file_name, element, rX, rY, rZ );

    double qO = -0.8476e0;
    double qH =  0.4238e0;
    std::vector<double> q(npar);

    for ( int i = 0; i < npar; ++i ){

        if ( element[i] == "O" ){
            q[i] = qO;
        }
        else if ( element[i] == "H" ){  
            q[i] = qH;
        }
        else{
            std::cout <<"Unrecognized chemical element appears"<< std::endl;
            return 0;
        }
    }               
 

    // Conversion of reduced units
    double r4pie;
    calc_unit_conversion_coeffs( r4pie );


    // Ewald parameter 
    double alpha = 0.55e0;

    // PME parameter
    int n_order = 6;
    int Kmax1 = 128;
    int Kmax2 = 128;
    int Kmax3 = 128;


    //A point charge is distributed to the number of surrounding grid points
    int points_n;
    //Parameters specific to Midtown-splines
    double zeta;
    if ( n_order == 4 ){ points_n = 32; zeta = 2.e0/sqrt(15.e0); }
    else if ( n_order == 6 ){ points_n = 88; zeta = 0.6503998e0; }
    
    //Record the relative coordinates of the grid points around the point charge
    std::vector<int> Midt_grid(points_n*3);
    record_Midt_grid( n_order, points_n, Midt_grid );

    cufftHandle plan;
    cufftCreate(&plan);
    cufftPlan3d(&plan, Kmax3, Kmax2, Kmax1, CUFFT_Z2Z);

    // reciprocal space energy
    double recE;
    // reciprocal space force
    std::vector<double> recFX(npar);
    std::vector<double> recFY(npar);
    std::vector<double> recFZ(npar);

clock_t start,end;
start = clock();
for ( int step = 0; step < 100; ++step ){ 

    calc_recSpace_Energy_Force(npar,
                               lX, lY, lZ,
                               r4pie, q, 
                               rX, rY, rZ, 
                               alpha,
                               n_order,
                               Kmax1, Kmax2, Kmax3,
                               plan,
                               points_n, 
                               zeta,
                               Midt_grid,
                               recE,
                               recFX, recFY, recFZ);
}
end = clock();
std::cout<<"time = "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;

    cufftDestroy(plan);

    std::ofstream E_O("CUDA_Energy.dat");
        E_O.precision(16);
        E_O << recE <<" 10J/mol "<< std::endl;
    E_O.close();


    std::ofstream dat_O("CUDA_Force.dat");
        dat_O.precision(16);
        dat_O <<" 10J/mol/angstrom "<< std::endl;

        for (int i= 0; i < npar; ++i ){

            dat_O << recFX[i] <<"  "
                  << recFY[i] <<"   "
                  << recFZ[i] << std::endl;
        }

    dat_O.close();

    return 0;
}

