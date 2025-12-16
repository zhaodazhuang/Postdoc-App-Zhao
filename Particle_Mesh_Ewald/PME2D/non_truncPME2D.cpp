#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<sstream>
#include<complex>
#include<algorithm> 


double round(double& x){ return (x>0.0)? floor(x+0.5):ceil(x-0.5); }


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


////////////////////////////////////////
//  Note: This FFT2D implements a standard forward FFT
//      but omits the  1/(N₁N₂) normalization in the inverse transform.
///////////////////////////////////////

void Fast_Fourier_Transform_2D(int& Kmax1,
                               int& Kmax2,
                               std::vector<int>& key1,
                               std::vector<int>& key2,
                               std::vector<std::complex<double>>& ww1,
                               std::vector<std::complex<double>>& ww2,
                               std::vector<std::vector<std::complex<double>>>& f,
                               bool& is_negetive){

    /*-------------------------------------------------
    *
    *  is_negetive = False
    *
    *  F(m1,m2,m3) = sum_{k1,k2} f(k1,k2)exp(i*2\pi k1*m1/K1)
    *                    *exp(i*2\pi k2*m2/K2)
    *
    *  is_negetive = True
    *
    *  f(k1,k2,k3) = sum_{m1,m2} F(m1,m2)exp(-i*2\pi k1*m1/K1)
    *                   *exp(-i*2\pi k2*m2/K2)
    *
    *---------------------------------------------------------------
    */

    // ceng = log_{2}(Kmax)
    // Kmax must be 2^{integer}
    int ceng1 = int( log(double(Kmax1)+1.e-10)/log(2) );
    int ceng2 = int( log(double(Kmax2)+1.e-10)/log(2) );

    // take conjugate of exponentials if required
    if (is_negetive ){
        for (int i = 0; i < Kmax1; ++i ){ ww1[i] = conj(ww1[i]); }
        for (int i = 0; i < Kmax2; ++i ){ ww2[i] = conj(ww2[i]); }
    }

    int m1, m2;
    int hN, nbit, m_add_hN;
    std::complex<double> Bm;


    // perform fourier transform in X direction
    m1 = 0;
    hN = Kmax1/2; // halfN

    for (int ceng = 0; ceng < ceng1; ++ceng ){
        while( m1 < Kmax1 ){
            for (int i = 0; i < hN; ++i ){
                nbit = key1[m1/hN]; m_add_hN = m1+hN;

                    for (m2 = 0; m2 < Kmax2; ++m2){

                        Bm = f[m_add_hN][m2]*ww1[nbit];
                        f[m_add_hN][m2] = f[m1][m2] - Bm;
                        f[m1][m2] = f[m1][m2] + Bm;

                    } // loop for m2
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

               for (m2 = 0; m2 < Kmax2; ++m2){
                   Bm = f[i][m2];
                   f[i][m2] = f[m1][m2];
                   f[m1][m2] = Bm;
    } } }



    // perform fourier transform in Y direction
    m2 = 0;
    hN = Kmax2/2; // halfN

    for (int ceng = 0; ceng < ceng2; ++ceng ){
        while( m2 < Kmax2 ){
            for (int i = 0; i < hN; ++i ){
                nbit = key2[m2/hN]; m_add_hN = m2+hN;

                    for (m1 = 0; m1 < Kmax1; ++m1){

                         Bm = f[m1][m_add_hN]*ww2[nbit];
                         f[m1][m_add_hN] = f[m1][m2] - Bm;
                         f[m1][m2] = f[m1][m2] + Bm;

                    } // loop for m1

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
                for (m1 = 0; m1 < Kmax1; ++m1){
                    Bm = f[m1][i];
                    f[m1][i] = f[m1][m2];
                    f[m1][m2] = Bm;
    } } }

    // take conjugate of exponentials if required
    if ( is_negetive ){
        for (int i = 0; i < Kmax1; ++i ){ ww1[i] = conj(ww1[i]); }
        for (int i = 0; i < Kmax2; ++i ){ ww2[i] = conj(ww2[i]); }
    }


    return;
}


////////////////////////////////////////
//  Note: This FFT3D implements a standard forward FFT
//      but omits the  1/(N₁N₂) normalization in the inverse transform.
///////////////////////////////////////

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



void calc_Euler_splines_coeffs(int& n_order,
                               int& Kmax1,
                               int& Kmax2,
                               int& Kmax3,
                               std::vector<std::complex<double>>& ww1,
                               std::vector<std::complex<double>>& ww2,
                               std::vector<std::complex<double>>& ww3,
                               std::vector<std::complex<double>>& Euler_b_1,
                               std::vector<std::complex<double>>& Euler_b_2,
                               std::vector<std::complex<double>>& Euler_b_3){            


    ////////////////////////////////////////////
    // see  eq.(4.1)
    ////////////////////////////////////////////
   
    /*----------------------------------------------
    *               0        0          0
    * 0  M_2(0) - M_3(0) - M_4(0) ... M_n(0)
    *           \        \
    * 1  M_2(1) - M_3(1) - M_4(1) ... M_n(1)
    *           \        \
    * 0  M_2(2) - M_3(2) - M_4(2) ... M_n(2)
    *           \        \
    * 0  M_2(3) - M_3(3) - M_4(3) ... M_n(3)      
    *           \        \
    * 0  M_2(4) - M_3(4) - M_4(4) ... M_n(4)
    *     .         .       .           .
    *     .         .       .           .
    *     .         .       .           .
    * 0  M_2(n) - M_3(n) - M_4(n) ... M_n(n)
    *
    *------------------------------------------------
    */
   
    /*-----------------------------------------------
    *
    *   M_k[u] = M_n(u);
    *
    *   For M_n of order n only need to obtain M_n(1), ... , M_n(n-1)
    *
    *   M_n(0) =  M_n(n) = M_n(n+1) = ... = 0
    *
    *   sum_{i=1}_{n-1} M_n(i) = 1
    *
    *   M_2[0] = 0, M_2[1] = 1, M_2[2] = 0
    *
    *   j cycles nust from large to small
    *
    * -----------------------------------------------
    */
   
    std::vector<double> M_k(n_order+1);

    // initialize : order k = 2
    M_k[0] = 0.e0; 
    M_k[1] = 1.e0; 
    M_k[2] = 0.e0;
   
    for (int k = 3; k <= n_order; ++k){  
        M_k[k] = 0.e0; // M_{k-1}(k) = 0 
   
        for (int j = k-1; j > 0; j = j-1 ){ // not need M_{n_order}(n_order)
            //           j               k-j
            // M_k(j) = --- M_{k-1}(j) + --- M_{k-1}(j-1) 
            //          k-1              k-1   
            M_k[j]=(double(j)*M_k[j]+double(k-j)*M_k[j-1])/double(k-1);
      
        } // loop for j
    } // loop for k
  


    /*----------------------------------------------- 
    *
    *  Calculate Euler Exponential Spline Coefficients  see eq.(4.4)
    *
    * 
    *  ww1[int] = exp(i*2pi*int/Kmax1)  
    *  ww2[int] = exp(i*2pi*int/Kmax2)
    *  ww3[int] = exp(i*2pi*int/Kmax3)
    *--------------------------------------------------------
    */

    std::complex<double> ccc;
   


    for (int mi = 0; mi < Kmax1; ++mi ){
        ccc = {0.e0, 0.e0};
    
        for (int k = 0; k < n_order-1; ++k ){
            ccc = ccc + M_k[k+1]*ww1[(mi*k)%Kmax1];  
        }
   
        Euler_b_1[mi] = ww1[( mi*(n_order-1) )%Kmax1]/ccc;
    } // loop for i
   


    for (int mi = 0; mi < Kmax2; ++mi ){
        ccc = {0.e0, 0.e0};
   
        for (int k = 0; k < n_order-1; ++k ){
            ccc = ccc + M_k[k+1]*ww2[(mi*k)%Kmax2];
        }
   
        Euler_b_2[mi] = ww2[( mi*(n_order-1) )%Kmax2]/ccc;
    } // loop for i
   



    for (int mi = 0; mi < Kmax3; ++mi ){
        ccc = {0.e0, 0.e0};
   
         for (int k = 0; k < n_order-1; ++k ){
             ccc = ccc + M_k[k+1]*ww3[(mi*k)%Kmax3];
         }
   
         Euler_b_3[mi] = ww3[( mi*(n_order-1) )%Kmax3]/ccc;
    } // loop for i

    return;     
}  




void calc_B_splines_coeffs(int& npar,
                           std::vector<double>& uX,
                           std::vector<double>& uY,
                           std::vector<double>& uZ,
                           int& n_order,
                           std::vector<std::vector<double>>& Mn_1,
                           std::vector<std::vector<double>>& Mn_2,
                           std::vector<std::vector<double>>& Mn_3,
                           std::vector<std::vector<double>>& dMn_1,
                           std::vector<std::vector<double>>& dMn_2,
                           std::vector<std::vector<double>>& dMn_3){

    /*---------------------------------------------
    *
    * Calculate B splines coefficients via the 
    * Cox-de Boor recursive formula
    *----------------------------------------------
    */

    /*----------------------------------------------
    *  uui = u_i - int(u_i)
    *
    * uui    M_2(uui)   - M_3(uui)   - M_4(uui)   ... M_n(uui)
    *                   \            \
    * 1-uui  M_2(uui+1) - M_3(uui+1) - M_4(uui+1) ... M_n(uui+1)
    *                   \            \
    * 0      M_2(uui+2) - M_3(uui+2) - M_4(uui+2) ... M_n(uui+2)
    *                   \            \
    * 0      M_2(uui+3) - M_3(uui+3) - M_4(uui+3) ... M_n(uui+3)
    *                   \            \
    * 0      M_2(uui+4) - M_3(uui+4) - M_4(uui+4) ... M_n(uui+4)
    *           .            .            .              .
    *           .            .            .              .
    *           .            .            .              .
    * 0      M_2(uui+n) - M_3(uui+n) - M_4(uui+n) ... M_n(uui+n)
    *
    * ------------------------------------------------------- 
    */
  
    /*----------------------------------------------------
    *
    *  M_2(uui) = uui, M_2(uui+1) = 1-uui
    *
    *  M_k(uui+k) = 0,  M_k(uui+k+1) = 0 ... 
    *
    *  M_n[i][j] = M_k(uui+j)
    *
    *  dM_n[i][j] = dM_k(uui+j)
    *------------------------------------------------------
    */
   
    for (int i = 0; i < npar; ++i ){
  
         // dM_2(uui+1) = 1
         dMn_1[i][1] = -1.e0; 
         dMn_2[i][1] = -1.e0; 
         dMn_3[i][1] = -1.e0;   
         // M_2(uui) = uui = u_i - int(u_i)
         Mn_1[i][0] = uX[i] - int(uX[i]);
         Mn_2[i][0] = uY[i] - int(uY[i]);
         Mn_3[i][0] = uZ[i] - int(uZ[i]);
         // M_2(uui+1) =  1 - uui
         Mn_1[i][1] = 1.e0 - uX[i] + int(uX[i]); 
         Mn_2[i][1] = 1.e0 - uY[i] + int(uY[i]);
         Mn_3[i][1] = 1.e0 - uZ[i] + int(uZ[i]);
    }
  
  

   
    /*------------------------------------------------------------
    * You can't use Mn[i][j] = M_{k}(uui+j) directly.
    *
    *
    * At the beginning of kth cycle : 
    * 
    * value stored in Mn[i][j] is M_{k-1}(uui+j) 
    *
    * value stored in Mn[i][j-1] is M_{k-1}(uui+j-1)
    *
    * At the End of the kth cycle
    *
    * value stored in Mn[i][j] is M_{k}(uui+j) 
    *
    * value stored in Mn[i][j-1] is M_{k}(uui+j-1)
    *
    * -----------------------------------------------------------
    */
    for (int k = 3; k <= n_order; ++k ){

        if ( k < n_order ) { 
            for (int i = 0; i < npar; ++i ){ 
                // M_{k-1}(uui+k) = 0
                Mn_1[i][k] =  0.e0;  
                Mn_2[i][k] =  0.e0;  
                Mn_3[i][k] =  0.e0;
            }
        }


        for (int j = k-1; j > 0; j=j-1){ 
            // not need M_{k}(uui+k)
  
            if (k == n_order){

                // dM_{n}(uui+n)=M_{n-1}(uui+n)-M_{n-1}(uui+n-1)
                for (int i = 0; i < npar; ++i ){ 
                    dMn_1[i][j] = Mn_1[i][j] - Mn_1[i][j-1];
                    dMn_2[i][j] = Mn_2[i][j] - Mn_2[i][j-1];
                    dMn_3[i][j] = Mn_3[i][j] - Mn_3[i][j-1];
                } 
            } 
     
  
            for (int i = 0; i < npar; ++i ){

                // uuX = uui + j
                double uuX = uX[i] + double(j) - int(uX[i]);
                double uuY = uY[i] + double(j) - int(uY[i]);
                double uuZ = uZ[i] + double(j) - int(uZ[i]);
                //                  uui+j                  k-(uui+j)
                // M_{k}(uui+j) = ------- M_{k-1}(uui+j) + --------- M_{k_1}(uui+j-1)
                //                   k-1                     k-1

                Mn_1[i][j] = (uuX*Mn_1[i][j] + (double(k)-uuX)*Mn_1[i][j-1])/double(k-1);
                Mn_2[i][j] = (uuY*Mn_2[i][j] + (double(k)-uuY)*Mn_2[i][j-1])/double(k-1);
                Mn_3[i][j] = (uuZ*Mn_3[i][j] + (double(k)-uuZ)*Mn_3[i][j-1])/double(k-1);

            } 
  
  
        } // loop for j
  
        /*------------------------------------------------------
        * 
        * The update of the j=0 item is different from the previous;
        *
        * bspX[i][-1] =  M_{k}(uui-1) = 0 
        *
        *-----------------------------------------------------------
        */
 
 
        if (k == n_order){

           for (int i = 0; i < npar; ++i ){ 
               // dM_{n}(uui) = M_{n-1}(uui)
               dMn_1[i][0] =  Mn_1[i][0];  
               dMn_2[i][0] =  Mn_2[i][0];
               dMn_3[i][0] =  Mn_3[i][0];
           }
        }

  
        for (int i = 0; i < npar; ++i ){ 
            //             uui
            // M_k(uui) = ---- M_{k-1}(uui)
            //             k-1
            Mn_1[i][0] = (uX[i]-int(uX[i]))*Mn_1[i][0]/double(k-1);
            Mn_2[i][0] = (uY[i]-int(uY[i]))*Mn_2[i][0]/double(k-1);
            Mn_3[i][0] = (uZ[i]-int(uZ[i]))*Mn_3[i][0]/double(k-1);
        }    
       

    } // loop for k
  
     
    return;
}



/////////////////////////////////////////////////////////
/*   calculate interaction mesh without truncation error */
/////////////////////////////////////////////////////////

void calc_interaction_mesh(double& lX,
		           double& lY,
		           double& lZ,
		           double& alpha,
		           int& Kmax1,
		           int& Kmax2,
		           int& Kmax3,
		           std::vector<int>& key1,
                           std::vector<int>& key2,
                           std::vector<int>& key3,
                           std::vector<std::complex<double>>& ww1,
                           std::vector<std::complex<double>>& ww2,
                           std::vector<std::complex<double>>& ww3,
                           std::vector<std::vector<std::vector<std::complex<double>>>>& interaction_mesh){

    /* unlike non-trunc PME3D, no truncation error exists on k3-axis, 
     * thus Kmax3 is not to be larger   */

    int pKmax1 = 128;
    int pKmax2 = 128;

    std::vector<std::vector<std::vector<std::complex<double>>>>
         interaction_mesh_dense(pKmax1, std::vector<std::vector<std::complex<double>>>
                               (pKmax2, std::vector<std::complex<double>>
                               (Kmax3) ) );

    std::vector<int> pkey1(pKmax1,0.e0);
    std::vector<int> pkey2(pKmax2,0.e0);
    std::vector<std::complex<double>> pww1(pKmax1,0.e0);
    std::vector<std::complex<double>> pww2(pKmax2,0.e0);

    calc_root_unity__bit_address(pKmax1, pKmax2, Kmax3,
                                 pkey1, pkey2, key3,
                                 pww1, pww2, ww3);

    int mm1, mm2, mm3;
    double kX, kY;
    double invlX = 1.e0/lX;
    double invlY = 1.e0/lY;
    double twopi = 2.e0*acos(-1.e0);
     
    double h, zgrid;
    double scalor = acosf(-1.e0)/(lX*lY);

    for (int m3 = 0; m3 < Kmax3; ++m3 ){
        mm3 = m3; if ( m3 > Kmax3/2 ){ mm3 = m3-Kmax3; };
        zgrid = lZ/double(Kmax3)*double(mm3);

        for (int m2 = 0; m2 < pKmax2; ++m2 ){
            mm2 = m2; if ( m2 > pKmax2/2 ){ mm2 = m2-pKmax2; };
            kY=twopi*double(mm2)*invlY;

            for (int m1 = 0; m1 < pKmax1; ++m1 ){
                mm1 = m1; if ( m1 > pKmax1/2 ){ mm1 = m1-pKmax1; };
                kX=twopi*double(mm1)*invlX;

                h = sqrtf(kY*kY + kX*kX);

                if ( h > 1.e-6 ){

                    interaction_mesh_dense[m1][m2][m3] = {( exp(-h*zgrid)*erfc(h/(2.e0*alpha)-alpha*zgrid) +
                                                   exp( h*zgrid)*erfc(h/(2.e0*alpha)+alpha*zgrid) )
                                                   /h*scalor, 0.e0  };
                }
                else{

                    interaction_mesh_dense[m1][m2][m3] = { 0.e0, 0.e0 };
                }
    }   }   } 
 
    bool is_negetive = false;

    for ( int m3 = 0; m3 < Kmax3; ++m3 ){

        std::vector<std::vector<std::complex<double>>>
             interaction_mesh2D_dense(pKmax1, std::vector<std::complex<double>>
             (pKmax2) );

        for ( int m2 = 0; m2 < pKmax2; ++m2 ){
            for ( int m1 = 0; m1 < pKmax1; ++m1 ){
            	interaction_mesh2D_dense[m1][m2] = interaction_mesh_dense[m1][m2][m3];
        }   }

        Fast_Fourier_Transform_2D(pKmax1, pKmax2,
                                  pkey1, pkey2,
                                  pww1, pww2,
                                  interaction_mesh2D_dense,
                                  is_negetive );

        for ( int m2 = 0; m2 < pKmax2; ++m2 ){
            for ( int m1 = 0; m1 < pKmax1; ++m1 ){
            	interaction_mesh_dense[m1][m2][m3] = interaction_mesh2D_dense[m1][m2];
        }   }

    }

    int multi1 = pKmax1/Kmax1;
    int multi2 = pKmax2/Kmax2;

    double coeffi = 1.e0/double(Kmax1*Kmax2*Kmax3);

    for ( int m3 = 0; m3 < Kmax3; ++m3 ){
	for ( int m2 = 0; m2 < Kmax2; ++m2 ){
 	    for ( int m1 = 0; m1 < Kmax1; ++m1 ){

                interaction_mesh[m1][m2][m3] =
		    interaction_mesh_dense[m1*multi1][m2*multi2][m3]*coeffi;
    }   }   }


    is_negetive = true;

    Fast_Fourier_Transform_3D(Kmax1, Kmax2, Kmax3,
                              key1,key2,key3,
                              ww1,ww2,ww3,
                              interaction_mesh,
                              is_negetive );

    return;
}



/////////////////////////////////////////////////////////
/*   calculate the energe and force in K-space */
/////////////////////////////////////////////////////////

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
                                std::vector<std::complex<double>>& Euler_b_1,
                                std::vector<std::complex<double>>& Euler_b_2,
                                std::vector<std::complex<double>>& Euler_b_3,
                                std::vector<std::vector<std::vector<std::complex<double>>>>& interaction_mesh,
                				double& E,
                                std::vector<double>& FX,
                                std::vector<double>& FY,
                                std::vector<double>& FZ){

 
    double invlX = 1.e0/lX;
    double invlY = 1.e0/lY;
    double invlZ = 1.e0/lZ;

    double twopi = 2.e0*acos(-1.e0); 

    for (int i = 0; i < npar; ++i ){
        FX[i] = 0.e0;  FY[i] = 0.e0;  FZ[i] = 0.e0; 
    }

  
    /*---------------------------------------------
    *  STEP1 : Obtain scaled fraction coordinates
    *  
    *      uX = rX/lX*Kamx1    
    *----------------------------------------------
    */
  
    std::vector<double> uX(npar);
    std::vector<double> uY(npar);
    std::vector<double> uZ(npar);
  
    for (int i = 0; i < npar; ++i ){
        double uu;
        uu = rX[i] - round(rX[i]*invlX)*lX; uX[i] = double(Kmax1)*(invlX*uu+0.5e0);
        uu = rY[i] - round(rY[i]*invlY)*lY; uY[i] = double(Kmax2)*(invlY*uu+0.5e0);
        uu = rZ[i]; uZ[i] = double(Kmax3)*(invlZ*uu);
    }

    
    /*-------------------------------------------------------------------
    * STEP 2: Obtain B-splines interpolation coefficient 
    *            M_n  <==>  Mn
    *         and its derivative d M_n  <==> dMn
    *
    *         Obtain charge_grid(k1,k2,k3) in real space see Eq.(4.6)
    *----------------------------------------------------------------
    */ 
   
    std::vector<std::vector<double>> Mn_1(npar, std::vector<double>(n_order));
    std::vector<std::vector<double>> Mn_2(npar, std::vector<double>(n_order));
    std::vector<std::vector<double>> Mn_3(npar, std::vector<double>(n_order));

    std::vector<std::vector<double>> dMn_1(npar, std::vector<double>(n_order));
    std::vector<std::vector<double>> dMn_2(npar, std::vector<double>(n_order));
    std::vector<std::vector<double>> dMn_3(npar, std::vector<double>(n_order));
  
    calc_B_splines_coeffs(npar, 
                          uX,uY,uZ, 
                          n_order,
                          Mn_1, Mn_2, Mn_3,
                          dMn_1, dMn_2, dMn_3); 
 

    std::vector<std::vector<std::vector<double>>>
    	charge_grid(Kmax1, std::vector<std::vector<double>>
                   (Kmax2, std::vector<double>
                   (Kmax3,0) ) );
   
    /*--------------------------------------------
    *
    * for j
    *
    *  bspX[i][j] = M_n( u_i-int(u_i) + j )
    *  
    *  k = int(u_i) - j
    *
    *  Q(int(u_i) - j) = qi*M_n( u_i-int(u_i) + j )
    *
    *  k = int(u_i) - j  <==> j = int(u_i) - k
    *
    *  Q(k) = qi*M_n( u_i-int(u_i) + int(u_i) - k )
    *       = qi*M_n( u_i - k ) 
    *
    *  Because of 0 <= u_i <= Kmax , 
    *            ( int(u_i) - j  ) >= Kmax  cannot happen;
    *
    *  When k < 0  M_n( u_i - k - n_1*Kmax )  
    *       n_1 = 1, k = k + Kmax      
    *---------------------------------------------
    */

  
    for (int i = 0; i < npar; ++i ){
        int k1, k2, k3;
        int iuXi=int(uX[i]), iuYi = int(uY[i]), iuZi = int(uZ[i]);    
  
        for (int j = 0; j < n_order; ++j ){
            k3 = iuZi-j; if( k3 >= Kmax3 ){ k3 = 0; } 
                         if( k3 < 0 ){ k3 = k3 + Kmax3;}
  
            for (int k = 0; k < n_order; ++k){
                k2 = iuYi-k; if( k2 >= Kmax2 ){ k2 = 0; }
                             if( k2 < 0 ){ k2 = k2 + Kmax2; }
  
                for (int l = 0; l < n_order; ++l ){
                    k1 = iuXi-l; if( k1 >= Kmax1 ){ k1 = 0; }
                                 if( k1 < 0 ){ k1 = k1 + Kmax1; }
          
                    charge_grid[k1][k2][k3] = charge_grid[k1][k2][k3]
                              + q[i]*Mn_1[i][l]*Mn_2[i][k]*Mn_3[i][j];
    }   }   }   }
/* 
    std::cout.precision(16);
    for (int k3 = 0; k3 < Kmax3; ++k3) {
    for (int k2 = 0; k2 < Kmax2; ++k2) {
    for (int k1 = 0; k1 < Kmax1; ++k1) {
          //  if (charge_grid_complex[k1][k2][k3].real() < 0.0) {
                std::cout << "[" << k1 << "][" << k2 << "][" << k3
                          << "] = " << charge_grid[k1][k2][k3] << std::endl;
   } } }    //        }
*/
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
        charge_grid_complex(Kmax1, std::vector<std::vector<std::complex<double>>>
                           (Kmax2, std::vector<std::complex<double>>
                           (Kmax3) ) );

  
    for (int k3 = 0; k3 < Kmax3; ++k3 ){
        for (int k2 = 0; k2 < Kmax2; ++k2 ){
            for (int k1 = 0; k1 < Kmax1; ++k1 ){
                charge_grid_complex[k1][k2][k3] = {charge_grid[k1][k2][k3],0.e0};
    }   }   }
  
/*
    std::cout.precision(16);
    for (int k3 = 0; k3 < Kmax3; ++k3) {
   for (int k2 = 0; k2 < Kmax2; ++k2) {
   for (int k1 = 0; k1 < Kmax1; ++k1) {
            if (charge_grid_complex[k1][k2][k3].real() < 0.0) {
                std::cout << "[" << k1 << "][" << k2 << "][" << k3
                          << "] = " << charge_grid_complex[k1][k2][k3].real() << std::endl;
   } } }            }
*/

    bool is_negetive = true;
    Fast_Fourier_Transform_3D(Kmax1, Kmax2, Kmax3,
                              key1,key2,key3, 
                              ww1,ww2,ww3,
                              charge_grid_complex,
             			      is_negetive );
 
 
    /*-------------------------------------------------------
    *  STEP 4: Obtain 
    *    
    *    interaction_mesh * B * F[Q](-m1, -m2, -m3)
    *
    *   now in the code charge_grid_complex <==>
    *    
    *    interaction_mesh(m1,m2,m3)* B(m1,m2, m3)*F[Q](-m1, -m2, -m3)
    *
    *------------------------------------------------------
    */ 
 
    double B3, B2, B;
  
    for (int m3 = 0; m3 < Kmax3; ++m3 ){
        B3 = (Euler_b_3[m3]*conj(Euler_b_3[m3])).real();

        for (int m2 = 0; m2 < Kmax2; ++m2 ){
            B2 = B3* ((Euler_b_2[m2]*conj(Euler_b_2[m2])).real());
  
            for (int m1 = 0; m1 < Kmax1; ++m1 ){
                B = B2*((Euler_b_1[m1]*conj(Euler_b_1[m1])).real());
 
                charge_grid_complex[m1][m2][m3] = B*interaction_mesh[m1][m2][m3]*charge_grid_complex[m1][m2][m3];                   
    }   }   }
 

    /*---------------------------------------------------------
    *
    *  STEP 5): Complete transformation  F[f]*g = f*F[g]  see Eq.(B3)
    *
    *  \sum_{m1,m2,m3} F[Q](m1,m2,m3)*
    *
    *    exp(-\pi^2 m^2/ beta^2)
    * {  ---------------------- * B(m1,m2, m3)*F[Q](-m1, -m2, -m3) }
    *            m^2
    *
    *  =
    *
    *  \sum_{m1,m2,m3} Q(k1,k2,k3)*
    *
    *      exp(-\pi^2 m^2/ beta^2)
    *  F[ ------------------------ * B(m1,m2, m3)*F[Q](-m1, -m2, -m3) ]
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
                              charge_grid_complex,
                              is_negetive );
 
    std::complex<double> E_complex = {0.e0,0.e0};
   
    for (int m1 = 0; m1 < Kmax1; ++m1 ){
        for (int m2 = 0; m2 < Kmax2; ++m2 ){
            for (int m3 = 0; m3 < Kmax3; ++m3 ){
                E_complex = E_complex + charge_grid[m1][m2][m3]*charge_grid_complex[m1][m2][m3];
    }   }   }

 
    E = E_complex.real()*0.5e0*r4pie;

      
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
  
    double F_coeffs = -r4pie;
    double force_field;
    double dMn1, dMn2, dMn3;
  
    for (int i = 0; i < npar; ++i ){
        int m1, m2, m3;
        int iuXi=int(uX[i]), iuYi = int(uY[i]),  iuZi = int(uZ[i]);
  
        for (int j = 0; j < n_order; ++j ){
            m3 = iuZi-j; if( m3 >= Kmax3 ){ m3 = 0; }
            if ( m3 < 0 ){ m3 = m3 + Kmax3; }
    
                for (int k = 0; k < n_order; ++k){
                    m2 = iuYi-k; if( m2 >= Kmax2 ){ m2 = 0; }
                    if ( m2 < 0 ){ m2 = m2 + Kmax2; }
  
                        for (int l = 0; l < n_order; ++l ){
                            m1 = iuXi-l; if( m1 >= Kmax1 ){ m1 = 0; }
                            if ( m1 < 0 ){ m1 = m1 + Kmax1; }
   
                                force_field = charge_grid_complex[m1][m2][m3].real();
                                dMn1 = force_field*dMn_1[i][l]* Mn_2[i][k]* Mn_3[i][j]*double(Kmax1);
                                dMn2 = force_field* Mn_1[i][l]*dMn_2[i][k]* Mn_3[i][j]*double(Kmax2);
                                dMn3 = force_field* Mn_1[i][l]* Mn_2[i][k]*dMn_3[i][j]*double(Kmax3);
                                FX[i] = FX[i] + F_coeffs*q[i]*dMn1*invlX;
                                FY[i] = FY[i] + F_coeffs*q[i]*dMn2*invlY;
                                FZ[i] = FZ[i] + F_coeffs*q[i]*dMn3*invlZ;
    }   }   }   }
  
    

    // Constrains the net residual force to zero; this operation is optional.

    /*    
    double FXtot = 0.e0, FYtot = 0.e0, FZtot = 0.e0;
  
    for (int i= 0; i < npar; ++i ){
  
        FXtot = FXtot + FX[i]; 
        FYtot = FYtot + FY[i];
        FZtot = FZtot + FZ[i];
    }
  
    std::cout << "residual force = " << FXtot << std::endl;
    std::cout << "residual force = " << FYtot << std::endl;
    std::cout << "residual force = " << FZtot << std::endl;

    
    FXtot = FXtot/double(npar);
    FYtot = FYtot/double(npar);
    FZtot = FZtot/double(npar);
  
    for (int i= 0; i < npar; ++i ){
  
        FX[i] = FX[i] - FXtot;
        FY[i] = FY[i] - FYtot;
        FZ[i] = FZ[i] - FZtot;
    }
  
    */
    return;
}


int main(){

    int npar = 2;
    double lX = 20.e0;
    double lY = 20.e0;
    double lZ = 20.e0;
    // The distance setting in the dimension 
    // where the interface is located must be at least twice 
    // the distance between the two farthest particles in that dimension.


    std::vector<std::string> element(npar);
    std::vector<double> rX(npar);
    std::vector<double> rY(npar);
    std::vector<double> rZ(npar);
    std::vector<double> q(npar);

    rX[0] = 1.e0; rY[0] = 2.e0; rZ[0] = 3.e0;
    rX[1] = 4.e0; rY[1] = 5.e0; rZ[1] = 6.e0;
    q[0] = 1.e0;
    q[1] = -1.e0;


    double r4pie;
    calc_unit_conversion_coeffs( r4pie ); // Conversion of reduced units

    double alpha = 0.5e0; // Ewald parameter



/////////////////////////////////////////////////
// energy and force in K-space 
//////////////////////////////////////////////////
     
    int n_order = 4; // interpolation order 
    int Kmax1 = 16;
    int Kmax2 = 16;
    int Kmax3 = 16; // PME mesh size 


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


    //Calculate Euler Exponential Spline Coefficients 
    std::vector<std::complex<double>> Euler_b_1(Kmax1);
    std::vector<std::complex<double>> Euler_b_2(Kmax2);
    std::vector<std::complex<double>> Euler_b_3(Kmax3);

    calc_Euler_splines_coeffs(n_order,  
                              Kmax1, Kmax2, Kmax3, 
                              ww1, ww2, ww3,
                              Euler_b_1, Euler_b_2, Euler_b_3);

    std::vector<std::vector<std::vector<std::complex<double>>>>
         interaction_mesh(Kmax1, std::vector<std::vector<std::complex<double>>>
                         (Kmax2, std::vector<std::complex<double>>
                         (Kmax3) ) );

    // calculacte the FT of interaction mesh 
    calc_interaction_mesh(lX, lY, lZ, 
        		          alpha,
	             		  Kmax1, Kmax2, Kmax3,
			              key1, key2, key3,
			              ww1, ww2, ww3,
                          interaction_mesh);
   
    double recE;
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
                               Euler_b_1, Euler_b_2, Euler_b_3,
			                   interaction_mesh,
                               recE, 
                               recFX, recFY, recFZ);


    std::cout.precision(16);
    std:: cout << recE << std::endl;

////////////////////////////////////////////////
// energy and force in R-space
///////////////////////////////////////////////
    double pi = acos(-1.e0);
    double dirE;
    double l;   
 
    dirE = 0.e0;
    for ( int n1 = -10; n1 <= 10 ; ++n1 ){
        for ( int n2 = -10; n2 <= 10 ; ++n2 ){

                for ( int i = 0; i < npar ; ++i ){
                for ( int j = 0; j < npar ; ++j ){

                    if ( n1 == 0 && n2 == 0 && i==j ){
                    
                        dirE = dirE - alpha/sqrt(pi)*q[i]*q[j];   
                    }   
                    else{

                        l = pow( (rX[i]- rX[j] + n1*lX),2.e0) 
                          + pow( (rY[i]- rY[j] + n2*lY),2.e0)
                          + pow( (rZ[i]- rZ[j]), 2.e0 ); 
                        l = sqrt(l);
                        dirE = dirE + 0.5e0*q[i]*q[j]*erfc(alpha*l)/l;
                    } 
                }   } 
    }  }  

    dirE = dirE*r4pie;
 
    double IBE = 0.e0;
    for ( int i = 0; i < npar ; ++i ){
        for ( int j = 0; j < npar ; ++j ){
 
            if ( j == i ){            
	        IBE = IBE - q[i]*q[j]/(alpha*sqrtf(pi));	   
	    }    
            else{
		double z = rZ[i]- rZ[j];
                IBE = IBE - q[i]*q[j]*(z*erf(alpha*z) + exp(-pow(alpha*z,2.e0))/(alpha*sqrtf(pi)));
	    }

    }   }

    IBE = IBE*pi/(lX*lY)*r4pie;



//////////////////////////////////////////////

    double unit_convert = 0.01e0;
           dirE = dirE*unit_convert;
           recE = recE*unit_convert;  
           IBE = IBE*unit_convert;


    std::ofstream E_O("Energy.dat");
    E_O.precision(16);

        E_O << "totE = " << dirE+recE+IBE << " KJ/mol" << std::endl;

        E_O <<" E^R= " << dirE <<"KJ/mol E^K= " 
	    << recE << " KJ/mol E^IB= "
	    << IBE << std::endl;

    E_O.close();


    std::ofstream dat_O("Force.dat");
    dat_O.precision(16);
    dat_O <<"unit : 10J/mol/angstrom "<< std::endl;

        for (int i= 0; i < npar; ++i ){

            dat_O << recFX[i] <<"  "
                  << recFY[i] <<"   "
                  << recFZ[i] << std::endl;
        }

    dat_O.close();

    return 0;
}

