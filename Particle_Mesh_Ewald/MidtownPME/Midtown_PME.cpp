#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<sstream>
#include<complex>
#include<algorithm> 
#include <time.h>

double round(double& x){ return (x>0.0)? floor(x+0.5):ceil(x-0.5); }



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
    double kbp23 = 1.38064852e0; // Kg*m^2/s^2/k
    double nan23=6.02214076e0;

    r4pie = pow(ele_to_clbp19,2.e0)*nan23/(4.e0*pi*epslonp12)*1.0e6;

    return;
}



void calc_root_unity__bit_address( int& Kmax1,
                                   int& Kmax2,
                                   int& Kmax3,
                                   std::vector<int>& key1,
                                   std::vector<int>& key2,
                                   std::vector<int>& key3,
                                   std::vector<std::complex<double>>& ww1,
                                   std::vector<std::complex<double>>& ww2,
                                   std::vector<std::complex<double>>& ww3){

    /*---------------------------------------------
    *
    * Calculate the primitive root of unity
    *   and  bit address which FFT requires
    *------------------------------------------
    */

    /*-------------------------------------
    *
    *   set reverse bit address arrays
    *
    *          0 1 2 3 4 5 6 7
    *                / \
    *               /   \
    *        0 2 4 6     1 3 5 7
    *          /\           /\
    *         /  \         /  \
    *      0 4    2 6    1 5  3 7
    *------------------------------------
    */

    // ceng = log_{2}(Kmax)
    // Kmax must be 2^{integer}
    int ceng1 = int( log(double(Kmax1)+1.e-10)/log(2) );
    int ceng2 = int( log(double(Kmax2)+2.e-20)/log(2) );
    int ceng3 = int( log(double(Kmax3)+3.e-30)/log(2) );
  
    bool is_ok = ( Kmax1 != pow(2,ceng1) || Kmax2 != pow(2,ceng2) || 
                   Kmax3 != pow(2,ceng3) );
 
    if ( is_ok ){ 
        std::cout <<" fft array not 2^integer; stop fft_efb_3"<< std::endl;
        exit(0);
    }


    int iii, jjj, jj2;


    for (int k = 0; k < Kmax1; ++k ){
        iii = 0; jjj = k;

        for (int j = 0; j < ceng1; j++ ){
            jj2 = jjj/2; iii = 2*(iii-jj2)+jjj; jjj = jj2;  
        }

        key1[k] = iii;
    }
   
 
    for (int k = 0; k < Kmax2; ++k ){
        iii = 0; jjj = k;

        for (int j = 0; j < ceng2; j++ ){
            jj2 = jjj/2; iii = 2*(iii-jj2)+jjj; jjj = jj2;
        }

        key2[k] = iii;
    }



    for (int k = 0; k < Kmax3; ++k ){
        iii = 0; jjj = k;

        for (int j = 0; j < ceng3; j++ ){
            jj2 = jjj/2; iii = 2*(iii-jj2)+jjj; jjj = jj2;
        }

        key3[k] = iii;
    }



    /*----------------------------------------------
    *  ww1[int] = exp(i*2pi*int/Kmax1)
    *  ww2[int] = exp(i*2pi*int/Kmax2)
    *  ww3[int] = exp(i*2pi*int/Kmax3)
    *
    *  exp(i*2pi*(Kmax-int)/Kmax) = exp(i*2pi)*exp(i*2pi*(-int)/Kmax3)
    *                             = exp(i*2pi*(-int)/Kmax3)
    *  So  ww[int] = conj( ww[Kmax-int] )
    *----------------------------------------------------
    */

    double tpi = 2.e0*acos(-1.e0);



    ww1[0] = {1.e0, 0.e0}; ww1[Kmax1/2] = {-1.e0, 0.e0};
   
    for (int i = 1; i < Kmax1/2; ++i ){

     double arg = tpi/double(Kmax1)*double(i);
     ww1[i] = {cos(arg), sin(arg)}; ww1[Kmax1-i] = conj(ww1[i]);
    }



    ww2[0] = {1.e0, 0.e0}; ww2[Kmax2/2] = {-1.e0, 0.e0};

    for (int i = 1; i < Kmax2/2; ++i ){

        double arg = tpi/double(Kmax2)*double(i);
        ww2[i] = {cos(arg), sin(arg)}; ww2[Kmax2-i] = conj(ww2[i]);
    }




    ww3[0] = {1.e0, 0.e0}; ww3[Kmax3/2] = {-1.e0, 0.e0};

    for ( int i = 1; i < Kmax3/2; ++i ){

        double arg = tpi/double(Kmax3)*double(i);
        ww3[i] = {cos(arg), sin(arg)}; ww3[Kmax3-i] = conj(ww3[i]);
    }

  
    return;
}


void record_Midt_grid(int& n_order,
                      int& points_n,
                      std::vector<std::vector<int>>& Midt_grid ){

    if ( n_order == 4 ){
  
        Midt_grid[0] = {1,1,1};
        Midt_grid[1] = {2,1,1};
        Midt_grid[2] = {1,2,1};
        Midt_grid[3] = {1,1,2};
    }
    else if ( n_order == 6 ){
  
  
        Midt_grid[0] = {1,1,1};
        Midt_grid[1] = {2,1,1};
        Midt_grid[2] = {1,2,1};
        Midt_grid[3] = {1,1,2};
        Midt_grid[4] = {3,1,1};
        Midt_grid[5] = {1,3,1};
        Midt_grid[6] = {1,1,3};
        Midt_grid[7] = {2,2,1};
        Midt_grid[8] = {2,1,2};
        Midt_grid[9] = {1,2,2};
        Midt_grid[10] = {2,2,2};
    }
  
    int points_nd8 = points_n/8;
    int points_nd4 = points_n/4;
    int points_nd2 = points_n/2;
    int a;
  
    for ( int i = points_nd8; i < points_nd4; ++i ){
  
        if ( 0.5e0 - Midt_grid[i-points_nd8][0] < 0 ){
            a = 1-Midt_grid[i-points_nd8][0]; 
        }

        Midt_grid[i] = { a, Midt_grid[i-points_nd8][1], Midt_grid[i-points_nd8][2] };
    }
  
  
    for ( int i = points_nd4; i < points_nd2; ++i ){
  
        if ( 0.5e0 - Midt_grid[i-points_nd4][1] < 0 ){
            a = 1-Midt_grid[i-points_nd4][1]; 
        }

        Midt_grid[i] = { Midt_grid[i-points_nd4][0], a, Midt_grid[i-points_nd4][2] };
    }
  
  
    for ( int i = points_nd2; i < points_n; ++i ){
  
      if ( 0.5e0 - Midt_grid[i-points_nd2][2] < 0 ){
          a = 1-Midt_grid[i-points_nd2][2]; 
      }

      Midt_grid[i] = { Midt_grid[i-points_nd2][0], Midt_grid[i-points_nd2][1], a };
    }
  
    return;
}


double Midt_L(double theta,
              double zeta){

    double L;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    L = -0.5e0*thetasq*theta + 0.5e0*thetasq
        -(9.e0*zetasq-2.e0)*theta/6.e0 + 0.5e0*zetasq;

    return L;
}



double Midt_dL(double theta,
               double zeta){

  double dL;

  dL = -1.5e0*theta*theta + theta
      -(9.e0*zeta*zeta-2.e0)/6.e0;

  return dL;
}



double Midt_R(double theta,
              double zeta){

    double R;

    R = theta*theta*theta/6.e0 + (3.e0*zeta*zeta-1.e0)*theta/6.e0;

    return R;
}


double Midt_dR(double theta,
               double zeta){

  double dR;

  dR = 0.5e0*theta*theta + (3.e0*zeta*zeta-1.e0)/6.e0;

  return dR;

}


double Midt_L111(double theta,
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



double Midt_dL111(double theta,
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



double Midt_L311(double theta,
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


double Midt_dL311(double theta,
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


double Midt_L211(double theta,
                 double zeta){

    double L;

    L = -0.25e0*Midt_L111(theta,zeta)
        - 2.5e0*Midt_L311(theta,zeta);

    return L;
}


double Midt_dL211(double theta,
                  double zeta){

    double dL;

    dL = -0.25e0*Midt_dL111(theta,zeta)
         - 2.5e0*Midt_dL311(theta,zeta);

    return dL;
}


double Midt_S(double theta,
              double zeta){

    double S;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    S = 0.5e0*( -thetasq*theta + thetasq
      -(3.e0*zetasq-2.e0)*theta + zetasq );

    return S;
}


double Midt_dS(double theta,
               double zeta){

    double dS;
    double thetasq = theta*theta;
    double zetasq = zeta*zeta;

    dS = 0.5e0*( -3.e0*thetasq + 2.e0*theta
       -(3.e0*zetasq-2.e0) );

    return dS;
}


void calc_Midt_coeffs(int& npar,
                      std::vector<double>& uX,
                      std::vector<double>& uY,
                      std::vector<double>& uZ,
                      int& n_order,
                      int& points_n,
                      double& zeta,
                      std::vector<std::vector<int>>& Midt_grid,
                      std::vector<std::vector<double>>& Midt_coeffs,
                      std::vector<std::vector<double>>& Midt_diff1,
                      std::vector<std::vector<double>>& Midt_diff2,
                      std::vector<std::vector<double>>& Midt_diff3){

    int points_nd8 = points_n/8;
    int points_nd4 = points_n/4;
    int points_nd2 = points_n/2;
  
    double thetaX, thetaY, thetaZ;
    double theta1, theta2, theta3;
    double sgn1, sgn2, sgn3; 
 
    if ( n_order == 4 ){
  
        for (int j = 0; j < npar; ++j ){
  
            thetaX = uX[j] - int(uX[j]);
            thetaY = uY[j] - int(uY[j]);
            thetaZ = uZ[j] - int(uZ[j]);
  
            for (int i = 0; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_L(theta1,zeta)*theta2*theta3
                                  + Midt_L(theta2,zeta)*theta1*theta3
                                  + Midt_L(theta3,zeta)*theta2*theta1;

                /////////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dL(theta1,zeta)*theta2*theta3
                                 + Midt_L(theta2,zeta)*theta3
                                 + Midt_L(theta3,zeta)*theta2 );

                Midt_diff2[i][j] = sgn2*(  Midt_L(theta1,zeta)*theta3
                                 + Midt_dL(theta2,zeta)*theta1*theta3
                                 + Midt_L(theta3,zeta)*theta1 );

                Midt_diff3[i][j] = sgn3*( Midt_L(theta1,zeta)*theta2
                                 + Midt_L(theta2,zeta)*theta1
                                 + Midt_dL(theta3,zeta)*theta1*theta2 );

            }

            for (int i = 1; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_R(theta1,zeta)*theta2*theta3;
                Midt_coeffs[i+1][j] = Midt_R(theta2,zeta)*theta1*theta3;
                Midt_coeffs[i+2][j] = Midt_R(theta3,zeta)*theta2*theta1;

                //////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*Midt_dR(theta1,zeta)*theta2*theta3;
                Midt_diff2[i][j] = sgn2*Midt_R(theta1,zeta)*theta3;
                Midt_diff3[i][j] = sgn3*Midt_R(theta1,zeta)*theta2;
        
                ///////////////////////////////////////////////////////////
                Midt_diff1[i+1][j] = sgn1*Midt_R(theta2,zeta)*theta3;
                Midt_diff2[i+1][j] = sgn2*Midt_dR(theta2,zeta)*theta1*theta3;
                Midt_diff3[i+1][j] = sgn3*Midt_R(theta2,zeta)*theta1;
        
                ///////////////////////////////////////////////////////////
                Midt_diff1[i+2][j] = sgn1*Midt_R(theta3,zeta)*theta2;
                Midt_diff2[i+2][j] = sgn2*Midt_R(theta3,zeta)*theta1;
                Midt_diff3[i+2][j] = sgn3*Midt_dR(theta3,zeta)*theta1*theta2;

            }
  
        } // loop for j
  
    } // n_order == 4
    else if ( n_order == 6 ){
  
        for (int j = 0; j < npar; ++j ){
  
            thetaX = uX[j] - int(uX[j]);
            thetaY = uY[j] - int(uY[j]);
            thetaZ = uZ[j] - int(uZ[j]);
  
            for (int i = 0; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_L111(theta1,zeta)*theta2*theta3
                                  + Midt_L111(theta2,zeta)*theta1*theta3
                                  + Midt_L111(theta3,zeta)*theta1*theta2
                + Midt_S(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta);

                ////////////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dL111(theta1,zeta)*theta2*theta3 
                                 + Midt_L111(theta2,zeta)*theta3
                                 + Midt_L111(theta3,zeta)*theta2
                + Midt_dS(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta) );

                Midt_diff2[i][j] = sgn2*( Midt_L111(theta1,zeta)*theta3
                                 + Midt_dL111(theta2,zeta)*theta1*theta3
                                 + Midt_L111(theta3,zeta)*theta1
                + Midt_S(theta1,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta3,zeta) );

                Midt_diff3[i][j] = sgn3*( Midt_L111(theta1,zeta)*theta2
                                 + Midt_L111(theta2,zeta)*theta1
                                 + Midt_dL111(theta3,zeta)*theta1*theta2
                + Midt_S(theta1,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta3,zeta) );
            }
  
            for (int i = 1; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_L211(theta1,zeta)*theta2*theta3
                + Midt_R(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta);

                Midt_coeffs[i+1][j] = Midt_L211(theta2,zeta)*theta1*theta3
                + Midt_R(theta2,zeta)*Midt_S(theta1,zeta)*Midt_S(theta3,zeta);

                Midt_coeffs[i+2][j] = Midt_L211(theta3,zeta)*theta2*theta1
                + Midt_R(theta3,zeta)*Midt_S(theta2,zeta)*Midt_S(theta1,zeta);

                //////////////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dL211(theta1,zeta)*theta2*theta3
                + Midt_dR(theta1,zeta)*Midt_S(theta2,zeta)*Midt_S(theta3,zeta) );

                Midt_diff2[i][j] = sgn2*( Midt_L211(theta1,zeta)*theta3
                + Midt_R(theta1,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta3,zeta) );

                Midt_diff3[i][j] = sgn3*( Midt_L211(theta1,zeta)*theta2
                + Midt_R(theta1,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta3,zeta) );

                ///////////////////////////////////////////////////////////////////
                Midt_diff1[i+1][j] = sgn1*( Midt_L211(theta2,zeta)*theta3
                + Midt_R(theta2,zeta)*Midt_dS(theta1,zeta)*Midt_S(theta3,zeta) );

                Midt_diff2[i+1][j] = sgn2*( Midt_dL211(theta2,zeta)*theta1*theta3
                + Midt_dR(theta2,zeta)*Midt_S(theta1,zeta)*Midt_S(theta3,zeta) );

                Midt_diff3[i+1][j] = sgn3*( Midt_L211(theta2,zeta)*theta1
                + Midt_R(theta2,zeta)*Midt_S(theta1,zeta)*Midt_dS(theta3,zeta) );

                //////////////////////////////////////////////////////////////////////
                Midt_diff1[i+2][j] = sgn1*( Midt_L211(theta3,zeta)*theta2
                + Midt_R(theta3,zeta)*Midt_S(theta2,zeta)*Midt_dS(theta1,zeta) );

                Midt_diff2[i+2][j] = sgn2*( Midt_L211(theta3,zeta)*theta1
                + Midt_R(theta3,zeta)*Midt_dS(theta2,zeta)*Midt_S(theta1,zeta) );

                Midt_diff3[i+2][j] = sgn3*( Midt_dL211(theta3,zeta)*theta1*theta2
                + Midt_dR(theta3,zeta)*Midt_S(theta2,zeta)*Midt_S(theta1,zeta) );
            }
  
            for (int i = 4; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_L311(theta1,zeta)*theta2*theta3;
                Midt_coeffs[i+1][j] = Midt_L311(theta2,zeta)*theta1*theta3;
                Midt_coeffs[i+2][j] = Midt_L311(theta3,zeta)*theta1*theta2;

                ///////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dL311(theta1,zeta)*theta2*theta3 );
                Midt_diff2[i][j] = sgn2*( Midt_L311(theta1,zeta)*theta3 );
                Midt_diff3[i][j] = sgn3*( Midt_L311(theta1,zeta)*theta2 );

                ///////////////////////////////////////////////////////////
                Midt_diff1[i+1][j] = sgn1*( Midt_L311(theta2,zeta)*theta3 );
                Midt_diff2[i+1][j] = sgn2*( Midt_dL311(theta2,zeta)*theta1*theta3 );
                Midt_diff3[i+1][j] = sgn3*( Midt_L311(theta2,zeta)*theta1 );

                //////////////////////////////////////////////////////////
                Midt_diff1[i+2][j] = sgn1*( Midt_L311(theta3,zeta)*theta2 );
                Midt_diff2[i+2][j] = sgn2*( Midt_L311(theta3,zeta)*theta1 );
                Midt_diff3[i+2][j] = sgn3*( Midt_dL311(theta3,zeta)*theta1*theta2 );
            }
  
            for (int i = 7; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_R(theta1,zeta)*
                                    Midt_R(theta2,zeta)*Midt_S(theta3,zeta);
                Midt_coeffs[i+1][j] = Midt_R(theta1,zeta)*
                                      Midt_R(theta3,zeta)*Midt_S(theta2,zeta);
                Midt_coeffs[i+2][j] = Midt_R(theta3,zeta)*
                                      Midt_R(theta2,zeta)*Midt_S(theta1,zeta);

                ////////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dR(theta1,zeta)*
                                   Midt_R(theta2,zeta)*Midt_S(theta3,zeta) );
                Midt_diff2[i][j] = sgn2*( Midt_R(theta1,zeta)*
                                   Midt_dR(theta2,zeta)*Midt_S(theta3,zeta) );
                Midt_diff3[i][j] = sgn3*( Midt_R(theta1,zeta)*
                                   Midt_R(theta2,zeta)*Midt_dS(theta3,zeta) );

                ////////////////////////////////////////////////////////////
                Midt_diff1[i+1][j] = sgn1*( Midt_dR(theta1,zeta)*
                                     Midt_R(theta3,zeta)*Midt_S(theta2,zeta) );
                Midt_diff2[i+1][j] = sgn2*( Midt_R(theta1,zeta)*
                                     Midt_R(theta3,zeta)*Midt_dS(theta2,zeta) );
                Midt_diff3[i+1][j] = sgn3*( Midt_R(theta1,zeta)*
                                     Midt_dR(theta3,zeta)*Midt_S(theta2,zeta) );

                //////////////////////////////////////////////////////////////
                Midt_diff1[i+2][j] = sgn1*( Midt_R(theta3,zeta)*
                                      Midt_R(theta2,zeta)*Midt_dS(theta1,zeta) );
                Midt_diff2[i+2][j] = sgn2*( Midt_R(theta3,zeta)*
                                      Midt_dR(theta2,zeta)*Midt_S(theta1,zeta) );
                Midt_diff3[i+2][j] = sgn3*( Midt_dR(theta3,zeta)*
                                      Midt_R(theta2,zeta)*Midt_S(theta1,zeta) );
            }
  
            for (int i = 10; i < points_n; i=i+points_nd8 ){
  
                if ( 0.5e0 - Midt_grid[i][0] < 0 ){ theta1 = thetaX; sgn1 = 1.e0; }
                else {  theta1 = 1.e0 - thetaX; sgn1 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][1] < 0 ){ theta2 = thetaY; sgn2 = 1.e0; }
                else {  theta2 = 1.e0 - thetaY; sgn2 = -1.e0; }
  
                if ( 0.5e0 - Midt_grid[i][2] < 0 ){ theta3 = thetaZ; sgn3 = 1.e0; }
                else {  theta3 = 1.e0 - thetaZ; sgn3 = -1.e0; }
  
                Midt_coeffs[i][j] = Midt_R(theta1,zeta)*
                                    Midt_R(theta2,zeta)*Midt_R(theta3,zeta);

                //////////////////////////////////////////////////////////
                Midt_diff1[i][j] = sgn1*( Midt_dR(theta1,zeta)*
                                   Midt_R(theta2,zeta)*Midt_R(theta3,zeta) );
                Midt_diff2[i][j] = sgn2*( Midt_R(theta1,zeta)*
                                   Midt_dR(theta2,zeta)*Midt_R(theta3,zeta) );
                Midt_diff3[i][j] = sgn3*( Midt_R(theta1,zeta)*
                                   Midt_R(theta2,zeta)*Midt_dR(theta3,zeta) );
            }
        } // loop for j
    }
  
    return;
}


void Fast_Fourier_Transform_3D(int& Kmax1,
                               int& Kmax2,
                               int& Kmax3,
                               std::vector<int>& key1,
                               std::vector<int>& key2,
                               std::vector<int>& key3,
                               std::vector<std::complex<double>>& ww1,
                               std::vector<std::complex<double>>& ww2,
                               std::vector<std::complex<double>>& ww3,
                               std::vector<std::vector<
                               std::vector<std::complex<double>>>>& f,
                               bool& is_negetive){

    /*-------------------------------------------------
    * Note that this FFT program is not strictly a Fast Fourier transform:
    *  
    *  is_negetive = False
    *
    *  F(m1,m2,m3) = sum_{k1,k2,k3} f(k1,k2,k3)exp(i*2\pi k1*m1/K1)
    *                    *exp(i*2\pi k2*m2/K2)*exp(i*2\pi k3*m3/K3)
    *
    *  is_negetive = True
    *
    *  f(k1,k2,k3) = sum_{m1,m2,m3} F(m1,m2,m3)exp(-i*2\pi k1*m1/K1)
    *                   *exp(-i*2\pi k2*m2/K2)*exp(-i*2\pi k3*m3/K3)
    *
    *---------------------------------------------------------------
    */


    /*--------------------------------------------------------------
    *
    * The logical structure of the code here must look at this picture: 
    *
    * https://www.researchgate.net/figure/Traffic-
    * patterns-in-the-Cooley-Tukey-FFT-algorithm-with-16-nodes_fig9_265591873
    *---------------------------------------------------------------------
    */ 


    // ceng = log_{2}(Kmax)
    // Kmax must be 2^{integer}
    int ceng1 = int( log(double(Kmax1)+1.e-10)/log(2) );
    int ceng2 = int( log(double(Kmax2)+1.e-10)/log(2) );
    int ceng3 = int( log(double(Kmax3)+1.e-10)/log(2) );
   
    // take conjugate of exponentials if required
    if (is_negetive ){ 
        for (int i = 0; i < Kmax1; ++i ){ ww1[i] = conj(ww1[i]); }
        for (int i = 0; i < Kmax2; ++i ){ ww2[i] = conj(ww2[i]); }
        for (int i = 0; i < Kmax3; ++i ){ ww3[i] = conj(ww3[i]); }
    }




    int m1, m2, m3;
    int hN, nbit, m_add_hN;
    std::complex<double> Bm;   


    // perform fourier transform in X direction
    m1 = 0; 
    hN = Kmax1/2; // halfN

    for (int ceng = 0; ceng < ceng1; ++ceng ){
        while( m1 < Kmax1 ){
            for (int i = 0; i < hN; ++i ){
                nbit = key1[m1/hN]; m_add_hN = m1+hN;

                for (m3 = 0; m3 < Kmax3; ++m3 ){
                    for (m2 = 0; m2 < Kmax2; ++m2){

                        Bm = f[m_add_hN][m2][m3]*ww1[nbit];
                        f[m_add_hN][m2][m3] = f[m1][m2][m3] - Bm;
                        f[m1][m2][m3] = f[m1][m2][m3] + Bm;

                    } // loop for m2
                } // loop fot m3
                m1 = m1 + 1;
            } // loop for i
            m1 = m1 + hN;
        } // while
        m1 = 0; hN = hN/2;
    } // loop for ceng

    // unscramble the fft using bit address array
    for (int i = 0; i < Kmax1; ++i ){
        m1 = key1[i];

        if (m1 > i ){

           for (m3 = 0; m3 < Kmax3; ++m3 ){
               for (m2 = 0; m2 < Kmax2; ++m2){
                   Bm = f[i][m2][m3];
                   f[i][m2][m3] = f[m1][m2][m3];
                   f[m1][m2][m3] = Bm;
    } } } }



    // perform fourier transform in Y direction
    m2 = 0;
    hN = Kmax2/2; // halfN

    for (int ceng = 0; ceng < ceng2; ++ceng ){
        while( m2 < Kmax2 ){
            for (int i = 0; i < hN; ++i ){
                nbit = key2[m2/hN]; m_add_hN = m2+hN;

                for (m3 = 0; m3 < Kmax3; ++m3 ){
                    for (m1 = 0; m1 < Kmax1; ++m1){

                         Bm = f[m1][m_add_hN][m3]*ww2[nbit];
                         f[m1][m_add_hN][m3] = f[m1][m2][m3] - Bm;
                         f[m1][m2][m3] = f[m1][m2][m3] + Bm;

                    } // loop for m1
                } // loop fot m3

                m2 = m2 + 1;
            } // loop for i
            m2 = m2 + hN;
        } // while
        m2 = 0; hN = hN/2;
    } // loop for ceng


    // unscramble the fft using bit address array
    for (int i = 0; i < Kmax2; ++i ){
        m2 = key2[i];
        if ( m2 > i ){
            for (m3 = 0; m3 < Kmax3; ++m3 ){
                for (m1 = 0; m1 < Kmax1; ++m1){
                    Bm = f[m1][i][m3];
                    f[m1][i][m3] = f[m1][m2][m3];
                    f[m1][m2][m3] = Bm;
    } } } }


    // perform fourier transform in Z direction
    m3 = 0;
    hN = Kmax3/2; // halfN

    for (int ceng = 0; ceng < ceng3; ++ceng ){
        while( m3 < Kmax3 ){
            for (int i = 0; i < hN; ++i ){
                nbit = key3[m3/hN]; m_add_hN = m3+hN;

                for (m2 = 0; m2 < Kmax2; ++m2 ){
                    for (m1 = 0; m1 < Kmax1; ++m1){

                        Bm = f[m1][m2][m_add_hN]*ww3[nbit];
                        f[m1][m2][m_add_hN] = f[m1][m2][m3] - Bm;
                        f[m1][m2][m3] = f[m1][m2][m3] + Bm;

                    } // loop for m1
                } // loop fot m2

                m3 = m3 + 1;
            } // loop for i
            m3 = m3 + hN;
        } // while
        m3 = 0; hN = hN/2;
    } // loop for ceng

    // unscramble the fft using bit address array
    for (int i = 0; i < Kmax3; ++i ){
        m3 = key3[i];
        if ( m3 > i ){
            for (m2 = 0; m2 < Kmax2; ++m2 ){
                for (m1 = 0; m1 < Kmax1; ++m1){
                    Bm = f[m1][m2][i];
                    f[m1][m2][i] = f[m1][m2][m3];
                    f[m1][m2][m3] = Bm;
    } } } }




    // take conjugate of exponentials if required
    if ( is_negetive ){
        for (int i = 0; i < Kmax1; ++i ){ ww1[i] = conj(ww1[i]); }
        for (int i = 0; i < Kmax2; ++i ){ ww2[i] = conj(ww2[i]); }
        for (int i = 0; i < Kmax3; ++i ){ ww3[i] = conj(ww3[i]); }
    }


    return;
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
                                std::vector<int>& key1,
                                std::vector<int>& key2,
                                std::vector<int>& key3,
                                std::vector<std::complex<double>>& ww1,
                                std::vector<std::complex<double>>& ww2,
                                std::vector<std::complex<double>>& ww3,
                                int& points_n,
                                double& zeta,
                                std::vector<std::vector<int>>& Midt_grid,
                                double& E,
                                std::vector<double>& FX,
                                std::vector<double>& FY,
                                std::vector<double>& FZ){



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

    for (int i = 0; i < npar; ++i ){
        FX[i] = 0.e0;  FY[i] = 0.e0;  FZ[i] = 0.e0; 
    }

  
    /*---------------------------------------------
    *  STEP1 : Obtain scaled fraction coordinates
    *  
    *      uX = rX/lX*Kamx1     see page 8579
    *----------------------------------------------
    */
  
    std::vector<double> uX(npar);
    std::vector<double> uY(npar);
    std::vector<double> uZ(npar);
  
    for (int i = 0; i < npar; ++i ){
        double uu;

        // test rX[i] - floor(rX[i]*invlX)*lX; uX[i] = double(Kmax1)*invlX*uu;
        uu = rX[i] - round(rX[i]*invlX)*lX; uX[i] = double(Kmax1)*(invlX*uu+0.5e0);
        uu = rY[i] - round(rY[i]*invlY)*lY; uY[i] = double(Kmax2)*(invlY*uu+0.5e0);
        uu = rZ[i] - round(rZ[i]*invlZ)*lZ; uZ[i] = double(Kmax3)*(invlZ*uu+0.5e0);  
    }
      
 
    /*-------------------------------------------------------------------
    * STEP 2: Obtain Midtown-splines interpolation coefficient 
    *        
    *
    *         Obtain charge_grid(k1,k2,k3) in real space see Eq.(4.6)
    *----------------------------------------------------------------
    */ 
  
    std::vector<std::vector<double>> Midt_coeffs(points_n, std::vector<double>(npar));

    std::vector<std::vector<double>> Midt_diff1(points_n, std::vector<double>(npar));
    std::vector<std::vector<double>> Midt_diff2(points_n, std::vector<double>(npar));
    std::vector<std::vector<double>> Midt_diff3(points_n, std::vector<double>(npar));
 
    calc_Midt_coeffs(npar, 
                     uX, uY, uZ,
                     n_order, 
                     points_n,
                     zeta, 
                     Midt_grid, 
                     Midt_coeffs,
                     Midt_diff1, Midt_diff2, Midt_diff3);

  
    std::vector<std::vector<std::vector<double>>>
    charge_grid(Kmax1, std::vector<std::vector<double>>
               (Kmax2, std::vector<double>
               (Kmax3,0) ) );


    for (int j = 0; j < npar; ++j ){
  
        int uXint, uYint, uZint;
        int k1, k2, k3;
  
        uXint = int(uX[j]);
        uYint = int(uY[j]);
        uZint = int(uZ[j]);
  
        for (int i = 0; i < points_n; ++i ){
  
            k1 = uXint + Midt_grid[i][0];
            if ( k1 < 0 ){ k1 = k1 + Kmax1; }
            if ( k1 >= Kmax1 ){ k1 = k1 - Kmax1; }
  
            k2 = uYint + Midt_grid[i][1];
            if ( k2 < 0 ){ k2 = k2 + Kmax2; }
            if ( k2 >= Kmax2 ){ k2 = k2 - Kmax2; }
  
            k3 = uZint + Midt_grid[i][2];
            if ( k3 < 0 ){ k3 = k3 + Kmax3; }
            if ( k3 >= Kmax3 ){ k3 = k3 - Kmax3; }
  
            charge_grid[k1][k2][k3] = charge_grid[k1][k2][k3] 
                                    + q[j]*Midt_coeffs[i][j];
        }
    }


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
  
    std::vector<std::vector<std::vector<std::complex<double>>>>
    ZZU(Kmax1, std::vector<std::vector<std::complex<double>>>
       (Kmax2, std::vector<std::complex<double>>
       (Kmax3) ) );

  
    for (int k3 = 0; k3 < Kmax3; ++k3 ){
        for (int k2 = 0; k2 < Kmax2; ++k2 ){
            for (int k1 = 0; k1 < Kmax1; ++k1 ){
                ZZU[k1][k2][k3] = {charge_grid[k1][k2][k3],0.e0};
    } } }
  
  
    bool is_negetive = true;
    Fast_Fourier_Transform_3D(Kmax1, Kmax2, Kmax3,
                              key1,key2,key3, 
                              ww1,ww2,ww3,
                              ZZU,
                              is_negetive );
 
 

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
 
    int mm1, mm2, mm3;
    double rkX, rkY, rkZ, rksq;
    double psi3, psi2, psi;

    /*-------------------------------------------------------------
    *  K space needs to be centered on 0 to maintain central symmetry, 
    *  so m needs to be translated (m = m-Kmax/2),
    *  although Kmax is an even number resulting in asymmetry, 
    *  but there is rec_cut to make it ultimately centrally symmetric
    *------------------------------------------------------------------
    */
  
    for (int m3 = 0; m3 < Kmax3; ++m3 ){
        mm3 = m3; if ( m3 > Kmax3/2 ){ mm3 = m3-Kmax3; };
        rkZ=tpi*double(mm3)*invlZ;
        psi3 = pow(double(mm3)/double(Kmax3),2.e0);
 
        for (int m2 = 0; m2 < Kmax2; ++m2 ){
            mm2 = m2; if ( m2 > Kmax2/2 ){ mm2 = m2-Kmax2; };
            rkY=tpi*double(mm2)*invlY;
            psi2 = psi3 + pow(double(mm2)/double(Kmax2),2.e0);

            for (int m1 = 0; m1 < Kmax1; ++m1 ){
                mm1 = m1; if ( m1 > Kmax1/2 ){ mm1 = m1-Kmax1; };
                rkX=tpi*double(mm1)*invlX;
                psi = psi2 + pow(double(mm1)/double(Kmax1),2.e0);

                rksq = rkX*rkX + rkY*rkY + rkZ*rkZ;
   
                if ( rksq > 1.e-6 && rksq <= rec_cutsq ){
                    ZZU[m1][m2][m3] = exp(tpi*tpi*zetasq*psi)
                                    *exp(r4alpsq*rksq)/rksq*ZZU[m1][m2][m3]; 
                }
                else{
  
                    ZZU[m1][m2][m3] = std::complex<double>(0, 0);
                } // if
    }  }  }
   

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
  
    is_negetive = false;
    Fast_Fourier_Transform_3D(Kmax1, Kmax2, Kmax3,
                              key1,key2,key3,
                              ww1,ww2,ww3,
                              ZZU,
                              is_negetive );
 
  
    std::complex<double> E_complex = {0.e0,0.e0};
   
    for (int m1 = 0; m1 < Kmax1; ++m1 ){
        for (int m2 = 0; m2 < Kmax2; ++m2 ){
            for (int m3 = 0; m3 < Kmax3; ++m3 ){
                E_complex = E_complex + charge_grid[m1][m2][m3]*ZZU[m1][m2][m3];
    } } }


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
  
    double F_coeffs = -2.e0*E_coeffs;
    double ZZUR;

    for (int j = 0; j < npar; ++j ){
  
        int uXint, uYint, uZint;
        int k1, k2, k3;
  
        uXint = int(uX[j]);
        uYint = int(uY[j]);
        uZint = int(uZ[j]);
  
        for (int i = 0; i < points_n; ++i ){
  
            k1 = uXint + Midt_grid[i][0];
            if ( k1 < 0 ){ k1 = k1 + Kmax1; }
            if ( k1 >= Kmax1 ){ k1 = k1 - Kmax1; }
      
            k2 = uYint + Midt_grid[i][1];
            if ( k2 < 0 ){ k2 = k2 + Kmax2; }
            if ( k2 >= Kmax2 ){ k2 = k2 - Kmax2; }
      
            k3 = uZint + Midt_grid[i][2];
            if ( k3 < 0 ){ k3 = k3 + Kmax3; }
            if ( k3 >= Kmax3 ){ k3 = k3 - Kmax3; }
      
            ZZUR = ZZU[k1][k2][k3].real();
      
            FX[j] = FX[j] + ZZUR*Midt_diff1[i][j];
            FY[j] = FY[j] + ZZUR*Midt_diff2[i][j];
            FZ[j] = FZ[j] + ZZUR*Midt_diff3[i][j];
 
        }
  
        FX[j] = F_coeffs*q[j]*FX[j]*double(Kmax1)*invlX;
        FY[j] = F_coeffs*q[j]*FY[j]*double(Kmax2)*invlY;
        FZ[j] = F_coeffs*q[j]*FZ[j]*double(Kmax3)*invlZ;
  
    }


    // Force the sum of forces to be 0
    double FXtot = 0.e0, FYtot = 0.e0, FZtot = 0.e0;
  
    for (int i= 0; i < npar; ++i ){
  
        FXtot = FXtot + FX[i]; 
        FYtot = FYtot + FY[i];
        FZtot = FZtot + FZ[i];
    }

    FXtot = FXtot/double(npar);
    FYtot = FYtot/double(npar);
    FZtot = FZtot/double(npar);
 
    for (int i= 0; i < npar; ++i ){
  
        FX[i] = FX[i] - FXtot;
        FY[i] = FY[i] - FYtot;
        FZ[i] = FZ[i] - FZtot;
    }
   
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


    //Calculate the primitive root of unity and bit address which FFT requires
    std::vector<int> key1(Kmax1);
    std::vector<int> key2(Kmax2);
    std::vector<int> key3(Kmax3);
    std::vector<std::complex<double>> ww1(Kmax1);
    std::vector<std::complex<double>> ww2(Kmax2);
    std::vector<std::complex<double>> ww3(Kmax3);
    calc_root_unity__bit_address(Kmax1, Kmax2, Kmax3, 
                                 key1, key2, key3, 
                                 ww1, ww2, ww3);

    //A point charge is distributed to the number of surrounding grid points
    int points_n;
    //Parameters specific to Midtown-splines
    double zeta;
    if ( n_order == 4 ){ points_n = 32; zeta = 2.e0/sqrt(15.e0); }
    else if ( n_order == 6 ){ points_n = 88; zeta = 0.6503998e0; }
    
    //Record the relative coordinates of the grid points around the point charge
    std::vector<std::vector<int>> Midt_grid(points_n,std::vector<int>(3));
    record_Midt_grid( n_order, points_n, Midt_grid );

    // reciprocal space energy
    double recE;
    // reciprocal space force
    std::vector<double> recFX(npar);
    std::vector<double> recFY(npar);
    std::vector<double> recFZ(npar); 
 
    calc_recSpace_Energy_Force(npar,
                               lX, lY, lZ,
                               r4pie, q, 
                               rX, rY, rZ, 
                               alpha,
                               n_order,
                               Kmax1, Kmax2, Kmax3,
                               key1, key2, key3, 
                               ww1, ww2, ww3,
                               points_n, 
                               zeta,
                               Midt_grid,
                               recE, 
                               recFX, recFY, recFZ);

    std::ofstream E_O("Energy.dat");
        E_O.precision(16);
        E_O << recE <<" 10J/mol "<< std::endl;
    E_O.close();


    std::ofstream dat_O("Force.dat");
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

