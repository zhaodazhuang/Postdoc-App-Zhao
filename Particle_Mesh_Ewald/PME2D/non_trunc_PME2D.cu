#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<sstream>
#include<complex>
#include<algorithm> 
#include<cufftXt.h>
#include<cuComplex.h>
#include <cublas_v2.h>

__device__ double my_round(double x){ return (x>0.0)? floor(x+0.5):ceil(x-0.5); }


void calc_unit_conversion_coeffs( double& r4pie ){

    double pi = acos(-1.e0);
    double ele_to_clbp19 = 1.602176634e0;
    double epslonp12 = 8.8541878128e0;
    //double kbp23 = 1.38064852e0; // Kg*m^2/s^2/k
    double nan23=6.02214076e0;

    r4pie = pow(ele_to_clbp19,2.e0)*nan23/(4.e0*pi*epslonp12)*1.0e6;

    return;
}


///////////////////////////////////
///  calculate root of unit 
//////////////////////////////////////

__global__ void compute_unit_root_gpu(cufftDoubleComplex* ww, int Kmax) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int half = Kmax/2;

    if (tid <= half) {
        double twopi = 2.e0*acos(-1.e0);  // 2π

        if (tid == 0) {
            ww[0].x = 1.0; ww[0].y = 0.0;
        } else if (tid == half) {
            ww[half].x = -1.0; ww[half].y = 0.0;
        } else {
            double arg = twopi * tid/Kmax;
            double cos_val = cos(arg);
            double sin_val = sin(arg);

            ww[tid].x = cos_val;
            ww[tid].y = sin_val;

            ww[Kmax-tid].x = cos_val;
            ww[Kmax-tid].y = -sin_val;
        }
    }
}


///////////////////////////////////////////////
/// calculate  Euler Splines corfficients 
/////////////////////////////////////////////////

__global__ void compute_Euler_gpu(int Kmax,
       		                  int n_order,
			          double* d_M_k,
			          cufftDoubleComplex* d_ww,
			          cufftDoubleComplex* d_Euler_b){

    int mi = blockIdx.x * blockDim.x + threadIdx.x; 
    
    if (mi < Kmax){

   	double sum_real = 0.e0;
	double sum_imag = 0.e0;

	for ( int k = 0; k < n_order-1 ; ++k ){
	    int idx = mi*k % Kmax;
            double Mk = d_M_k[k+1];

	    sum_real += Mk*d_ww[idx].x;
	    sum_imag += Mk*d_ww[idx].y;
	}

        int idx = mi*(n_order-1) % Kmax;

        d_Euler_b[mi].x = d_ww[idx].x * sum_real + d_ww[idx].y * sum_imag;
	d_Euler_b[mi].y = d_ww[idx].y * sum_real - d_ww[idx].x * sum_imag;

	double denom = sum_real*sum_real + sum_imag*sum_imag;

	d_Euler_b[mi].x = d_Euler_b[mi].x /denom;
	d_Euler_b[mi].y = d_Euler_b[mi].y /denom;
    }

}




void calc_Euler_splines_coeffs(int n_order,
                               int Kmax1,
                               int Kmax2,
                               int Kmax3,
                               cufftDoubleComplex* d_Euler_b_1,
                               cufftDoubleComplex* d_Euler_b_2,
                               cufftDoubleComplex* d_Euler_b_3){            


    cufftDoubleComplex *d_ww1, *d_ww2, *d_ww3;
    cudaMalloc(&d_ww1, Kmax1*sizeof(cufftDoubleComplex));
    cudaMalloc(&d_ww2, Kmax2*sizeof(cufftDoubleComplex));
    cudaMalloc(&d_ww3, Kmax3*sizeof(cufftDoubleComplex));

    compute_unit_root_gpu<<<(Kmax1/2+1+255)/256, 256>>>(d_ww1, Kmax1);
    compute_unit_root_gpu<<<(Kmax2/2+1+255)/256, 256>>>(d_ww2, Kmax2);
    compute_unit_root_gpu<<<(Kmax3/2+1+255)/256, 256>>>(d_ww3, Kmax3);


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
  
 
    double *d_M_k;
    cudaMalloc(&d_M_k, (n_order+1)*sizeof(double));

    cudaMemcpy(d_M_k, M_k.data(), (n_order + 1) * sizeof(double), 
                                          cudaMemcpyHostToDevice);


    compute_Euler_gpu<<<(Kmax1+255)/256, 256 >>>(Kmax1, n_order,
						 d_M_k, d_ww1, d_Euler_b_1);

    compute_Euler_gpu<<<(Kmax2+255)/256, 256 >>>(Kmax2, n_order,
                                                 d_M_k, d_ww2, d_Euler_b_2);

    compute_Euler_gpu<<<(Kmax3+255)/256, 256 >>>(Kmax3, n_order,
                                                 d_M_k, d_ww3, d_Euler_b_3);

    return;     
}  


/////////////////////////////////////////////////////////
/*   calculate interaction mesh without truncation error */
/////////////////////////////////////////////////////////


__global__ void compute_interaction_mesh_dense_gpu(int pKmax1,
		                                   int pKmax2,
						   int Kmax3,
		                                   double invlX,
	                                           double invlY,
						   double lZ,
						   double alpha,
                                                   double twopi,
                                                   double scalor,
                                                   cufftDoubleComplex* d_interaction_mesh_dense){

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = pKmax1 * pKmax2 * Kmax3;
    
    if (tid < total) {
        int m1 = tid % pKmax1;
        int m2 = (tid / pKmax1) % pKmax2;
        int m3 = tid / (pKmax1 * pKmax2);

        if ( m3 > Kmax3/2 ){ m3 = m3-Kmax3; };
            double zgrid = lZ/double(Kmax3)*double(m3);
        
        if ( m2 > pKmax2/2 ){ m2 = m2-pKmax2; };
            double kY = twopi*double(m2)*invlY;

        if ( m1 > pKmax1/2 ){ m1 = m1-pKmax1; };
            double kX = twopi*double(m1)*invlX;

            double h = sqrt(kY*kY + kX*kX);

            if ( h > 1.e-6 ){

                d_interaction_mesh_dense[tid].x = ( exp(-h*zgrid)*erfc(h/(2.e0*alpha)-alpha*zgrid) +
                                                    exp( h*zgrid)*erfc(h/(2.e0*alpha)+alpha*zgrid) )/h*scalor;
		d_interaction_mesh_dense[tid].y = 0.e0;
            }
            else{

                d_interaction_mesh_dense[tid].x = 0.e0;
		d_interaction_mesh_dense[tid].y = 0.e0;
            }
    }

}


__global__ void  extract_interaction_mesh(int pKmax1,
	                                  int pKmax2,
                                          int Kmax1,
					  int Kmax2,
					  int Kmax3,
                                          int multi1,
					  int multi2,
                                          double coeffi,
                                          cufftDoubleComplex *d_interaction_mesh_dese,
                                          cufftDoubleComplex *d_interaction_mesh){

    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int total = Kmax1*Kmax2*Kmax3;

    if ( tid < total ){
	int m1 = tid % Kmax1;
	int m2 = (tid/ Kmax1) % Kmax2;
	int m3 = tid / (Kmax1*Kmax2);
	
        int id_dense = m3*(pKmax1*pKmax2) + m2*multi2*pKmax1 + m1*multi1;

	d_interaction_mesh[tid].x = d_interaction_mesh_dese[id_dense].x*coeffi;
	d_interaction_mesh[tid].y = d_interaction_mesh_dese[id_dense].y*coeffi;
    }

}



void calc_interaction_mesh(double lX,
		           double lY,
		           double lZ,
		           double alpha,
		           int Kmax1,
		           int Kmax2,
		           int Kmax3,
			   cufftHandle plan,
                           cufftDoubleComplex* d_interaction_mesh){

    /* unlike non-trunc PME3D, no truncation error exists on k3-axis, 
     * thus Kmax3 is not to be larger   */

    int pKmax1 = 128;
    int pKmax2 = 128;

    cufftDoubleComplex *d_interaction_mesh_dense;
    cudaMalloc(&d_interaction_mesh_dense, pKmax1*pKmax2*Kmax3*sizeof(cufftDoubleComplex));

    double invlX = 1.e0/lX;
    double invlY = 1.e0/lY;
    double twopi = 2.e0*acosf(-1.e0);
    double scalor = acosf(-1.e0)/(lX*lY);

    compute_interaction_mesh_dense_gpu<<<(pKmax1*pKmax2*Kmax3+255)/256,256>>>(pKmax1, pKmax2, Kmax3,
		                                                              invlX, invlY, lZ,
		                                                              alpha, twopi, scalor,
									      d_interaction_mesh_dense);

    cufftHandle plan_2d;
    cufftPlan2d(&plan_2d, pKmax2, pKmax1, CUFFT_Z2Z);

    for ( int m3 = 0; m3 < Kmax3; ++m3 ){

        cufftDoubleComplex* d_begin = d_interaction_mesh_dense + m3 * (pKmax1 * pKmax2);

	cufftExecZ2Z(plan_2d, d_begin, d_begin, CUFFT_INVERSE);
    }


    int multi1 = pKmax1/Kmax1;
    int multi2 = pKmax2/Kmax2;

    double coeffi = 1.e0/double(Kmax1*Kmax2*Kmax3);

    extract_interaction_mesh<<<(pKmax1*pKmax2*Kmax3+255)/256, 256>>>(pKmax1, pKmax2,
		                                                     Kmax1, Kmax2, Kmax3,
		                                                     multi1, multi2,
		                                                     coeffi,
						                     d_interaction_mesh_dense,
						                     d_interaction_mesh);

    cufftExecZ2Z(plan, d_interaction_mesh, d_interaction_mesh, CUFFT_FORWARD);

    return;
}



__global__ void  compute_frac_coor(int npar,
                                   int Kmax1,
				   int Kmax2,
				   int Kmax3,
                                   double lX,
				   double lY,
				   double lZ,
                                   double* d_rX, 
				   double* d_rY,
				   double* d_rZ,
                                   double* d_uX,
				   double* d_uY,
				   double* d_uZ){

    int i = blockIdx.x*blockDim.x + threadIdx.x;

    if ( i < npar ){
	double uu = d_rX[i] - my_round(d_rX[i]/lX)*lX;
        d_uX[i] = double(Kmax1)*(uu/lX+0.5e0);

        uu = d_rY[i] - my_round(d_rY[i]/lY)*lY;
        d_uY[i] = double(Kmax2)*(uu/lY+0.5e0);

        uu = d_rZ[i];
        d_uZ[i] = double(Kmax3)*(uu/lZ);
    	
    }
}


////////////////////////////////////////
/// calculate B-splines corfficent
//////////////////////////////////////////

__global__ void compute_B_splines_coeffs_gpu(int npar,
	                                     int n_order,
					     double* d_u,
                                             double* d_Mn,
                                             double* d_dMn){

    int i = blockIdx.x*blockDim.x + threadIdx.x;

    if ( i > npar ) return;

    int idx = i * n_order + 1; 
    d_dMn[idx] = -1.e0;

    idx = i * n_order + 0;
    d_Mn[idx] = d_u[i] - int(d_u[i]);
         
    idx = i * n_order + 1;
    d_Mn[idx] = 1.e0 - d_u[i] + int(d_u[i]);
    

    for (int k = 3; k <= n_order; ++k ){

         if ( k < n_order ){
	    d_Mn[i*n_order+k] = 0.e0;
         }

        for (int j = k-1; j > 0; j=j-1){

	    int idj = i * n_order + j;

            if ( k == n_order ){
 
 	        d_dMn[idj] = d_Mn[idj] - d_Mn[idj-1];
            }
        

            double uu = d_u[i] + double(j) - int(d_u[i]);
            
            d_Mn[idj] = (uu*d_Mn[idj] + (double(k)-uu)*d_Mn[idj-1])/double(k-1);
 
	}


	int id0 = i * n_order + 0;


        if ( k == n_order ){

           d_dMn[id0] =  d_Mn[id0];
        }
 

        d_Mn[id0] = (d_u[i]-int(d_u[i]))*d_Mn[id0]/double(k-1); 
 
	} 
}


__global__ void assign_charge_grid_gpu(int npar,
                                          int n_order,
                                          double* d_uX,
                                          double* d_uY,
                                          double* d_uZ,
                                          double* d_q,
                                          double* d_Mn_1,
                                          double* d_Mn_2,
                                          double* d_Mn_3,
                                          int Kmax1,
                                          int Kmax2,
                                          int Kmax3,
                                          double* d_charge_grid){
  
    int i = blockIdx.x;       
    int point_i = threadIdx.x;
    
    if (i >= npar) return;    
    if (point_i >= n_order*n_order*n_order) return;      
   
    int j = point_i / (n_order * n_order);                
    int k = (point_i % (n_order * n_order)) / n_order;    
    int l = point_i % n_order;                            

    int iuXi = int(d_uX[i]);
    int iuYi = int(d_uY[i]);
    int iuZi = int(d_uZ[i]);

    int k3 = iuZi - j;
    if (k3 >= Kmax3) k3 = 0;
    if (k3 < 0) k3 = k3 + Kmax3;

    int k2 = iuYi - k;
    if (k2 >= Kmax2) k2 = 0;
    if (k2 < 0) k2 = k2 + Kmax2;

    int k1 = iuXi - l;
    if (k1 >= Kmax1) k1 = 0;
    if (k1 < 0) k1 = k1 + Kmax1;

    int idx_grid = k3 * (Kmax1 * Kmax2) + k2 * Kmax1 + k1;

    int idx_base = i * n_order;
    double coeff = d_Mn_1[idx_base + l] * d_Mn_2[idx_base + k] * d_Mn_3[idx_base + j];

    double charge_contrib = d_q[i] * coeff;
    atomicAdd(&d_charge_grid[idx_grid], charge_contrib);
}


__global__ void convert_real_to_complex_kernel(int Kmax1,
		                               int Kmax2,
					       int Kmax3,
		                               double* d_charge_grid_real,
                                               cufftDoubleComplex* d_charge_grid_complex) {
   
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = Kmax1 * Kmax2 * Kmax3;
    
    if (tid < total) {
        d_charge_grid_complex[tid].x = d_charge_grid_real[tid];  
        d_charge_grid_complex[tid].y = 0.0;         
    }
}


__global__ void apply_euler_and_interaction_kernel(int Kmax1, 
		                                   int Kmax2,
						   int Kmax3,
                                                   cufftDoubleComplex* d_charge_grid_complex,
                                                   cufftDoubleComplex* d_interaction_mesh,
                                                   cufftDoubleComplex* d_Euler_b_1,
                                                   cufftDoubleComplex* d_Euler_b_2,
                                                   cufftDoubleComplex* d_Euler_b_3){
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = Kmax1 * Kmax2 * Kmax3;
    
    if (tid >= total) return;
    
    int m1 = tid % Kmax1;
    int m2 = (tid / Kmax1) % Kmax2;
    int m3 = tid / (Kmax1 * Kmax2);
    
    cufftDoubleComplex e3 = d_Euler_b_3[m3];
    double B3 = e3.x * e3.x + e3.y * e3.y;
    
    cufftDoubleComplex e2 = d_Euler_b_2[m2];
    double B2 = B3 * (e2.x * e2.x + e2.y * e2.y);
    
    cufftDoubleComplex e1 = d_Euler_b_1[m1];
    double B = B2 * (e1.x * e1.x + e1.y * e1.y);
    
    cufftDoubleComplex interaction = d_interaction_mesh[tid];
    cufftDoubleComplex charge_val = d_charge_grid_complex[tid];
    
    double tmp_real = interaction.x * charge_val.x - interaction.y * charge_val.y;
    double tmp_imag = interaction.x * charge_val.y + interaction.y * charge_val.x;
    
    d_charge_grid_complex[tid].x = B * tmp_real;
    d_charge_grid_complex[tid].y = B * tmp_imag;
}


__global__ void compute_forces_kernel(int npar,
    				      int n_order,
                                      double invlX,
				      double invlY,
				      double invlZ,
                                      int Kmax1, 
				      int Kmax2,
				      int Kmax3,
                                      double F_coeffs,   
                                      double* d_uX,
				      double* d_uY,
				      double* d_uZ,
                                      double* d_q,
                                      double* d_Mn_1, 
				      double* d_Mn_2,
				      double* d_Mn_3,
                                      double* d_dMn_1,
				      double* d_dMn_2,
				      double* d_dMn_3,
                                      cufftDoubleComplex* d_charge_grid_complex,
                                      double* d_FX,
				      double* d_FY, 
				      double* d_FZ){
    int i = blockIdx.x;         
    int point_i = threadIdx.x; 

    int p3 = n_order * n_order * n_order;
    if (i >= npar || point_i >= p3) return;

    int j = point_i / (n_order * n_order);
    int k = (point_i % (n_order * n_order)) / n_order;
    int l = point_i % n_order;

    int iuXi = int(d_uX[i]);
    int iuYi = int(d_uY[i]);
    int iuZi = int(d_uZ[i]);

    int m3 = iuZi - j;
    if (m3 >= Kmax3) m3 = 0;
    if (m3 < 0) m3 += Kmax3;

    int m2 = iuYi - k;
    if (m2 >= Kmax2) m2 = 0;
    if (m2 < 0) m2 += Kmax2;

    int m1 = iuXi - l;
    if (m1 >= Kmax1) m1 = 0;
    if (m1 < 0) m1 += Kmax1;

    int idx = m3 * (Kmax1 * Kmax2) + m2 * Kmax1 + m1;

    double force_field = d_charge_grid_complex[idx].x;

    int base = i * n_order;
    double Mn1 = d_Mn_1[base + l];
    double Mn2 = d_Mn_2[base + k];
    double Mn3 = d_Mn_3[base + j];
    double dMn1 = d_dMn_1[base + l];
    double dMn2 = d_dMn_2[base + k];
    double dMn3 = d_dMn_3[base + j];

    double qi = d_q[i];
    double dMn1_val = force_field * dMn1 * Mn2 * Mn3 * double(Kmax1);
    double dMn2_val = force_field * Mn1 * dMn2 * Mn3 * double(Kmax2);
    double dMn3_val = force_field * Mn1 * Mn2 * dMn3 * double(Kmax3);

    atomicAdd(&d_FX[i], F_coeffs * qi * dMn1_val * invlX);
    atomicAdd(&d_FY[i], F_coeffs * qi * dMn2_val * invlY);
    atomicAdd(&d_FZ[i], F_coeffs * qi * dMn3_val * invlZ);
}



/////////////////////////////////////////////////////////
/*   calculate the energe and force in K-space */
/////////////////////////////////////////////////////////

void calc_recSpace_Energy_Force(int npar,
                                double lX,
                                double lY,
                                double lZ,
                                double r4pie,
                                double* d_q,
                                double* d_rX,
                                double* d_rY,
                                double* d_rZ,
                                double alpha,
                                int n_order,
                                int Kmax1,
                                int Kmax2,
                                int Kmax3,
                                cufftDoubleComplex* d_Euler_b_1,
                                cufftDoubleComplex* d_Euler_b_2,
                                cufftDoubleComplex* d_Euler_b_3,
                                cufftDoubleComplex* d_interaction_mesh,
				cufftHandle plan,
				double& E,
                                double* d_FX,
                                double* d_FY,
                                double* d_FZ){

 
    double invlX = 1.e0/lX;
    double invlY = 1.e0/lY;
    double invlZ = 1.e0/lZ;

    double twopi = 2.e0*acos(-1.e0); 


    /*---------------------------------------------
    *  STEP1 : Obtain scaled fraction coordinates
    *  
    *      uX = rX/lX*Kamx1    
    *----------------------------------------------
    */
  
    double *d_uX, *d_uY, *d_uZ;

    cudaMalloc(&d_uX, npar*sizeof(double));
    cudaMalloc(&d_uY, npar*sizeof(double));
    cudaMalloc(&d_uZ, npar*sizeof(double));
		   
    compute_frac_coor<<<(npar+255)/256, 256>>>(npar,
		                           Kmax1, Kmax2, Kmax3,
		                           lX, lY, lZ,
					   d_rX, d_rY, d_rZ,
		                           d_uX, d_uY, d_uZ);
    

 
    /*-------------------------------------------------------------------
    * STEP 2: Obtain B-splines interpolation coefficient 
    *            M_n  <==>  Mn
    *         and its derivative d M_n  <==> dMn
    *
    *         Obtain charge_grid(k1,k2,k3) in real space see Eq.(4.6)
    *----------------------------------------------------------------
    */ 

    double *d_Mn_1, *d_Mn_2, *d_Mn_3;
    double *d_dMn_1, *d_dMn_2, *d_dMn_3;    

    cudaMalloc(&d_Mn_1, npar*n_order*sizeof(double));
    cudaMalloc(&d_Mn_2, npar*n_order*sizeof(double));
    cudaMalloc(&d_Mn_3, npar*n_order*sizeof(double));

    cudaMalloc(&d_dMn_1, npar*n_order*sizeof(double));
    cudaMalloc(&d_dMn_2, npar*n_order*sizeof(double));
    cudaMalloc(&d_dMn_3, npar*n_order*sizeof(double));


    compute_B_splines_coeffs_gpu<<<(npar+255)/256, 256>>>(npar, n_order,
                                                          d_uX,
                                                          d_Mn_1,
                                                          d_dMn_1);


    compute_B_splines_coeffs_gpu<<<(npar+255)/256, 256>>>(npar, n_order,
                                                          d_uY,
                                                          d_Mn_2,
                                                          d_dMn_2);

    compute_B_splines_coeffs_gpu<<<(npar+255)/256, 256>>>(npar, n_order,
                                                          d_uZ,
                                                          d_Mn_3,
                                                          d_dMn_3);

    double *d_charge_grid;
    cudaMalloc(&d_charge_grid, Kmax1*Kmax2*Kmax3*sizeof(double));
    cudaMemset(d_charge_grid, 0, Kmax1*Kmax2*Kmax3*sizeof(double));

    int p3 = n_order*n_order*n_order;
    assign_charge_grid_gpu<<<npar, p3>>>(npar, n_order,
						      d_uX, d_uY, d_uZ, d_q,
                                                      d_Mn_1, d_Mn_2, d_Mn_3,
                                                      Kmax1, Kmax2, Kmax3,
						      d_charge_grid);
    cudaDeviceSynchronize(); 

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
 
    cufftDoubleComplex *d_charge_grid_complex;
    cudaMalloc(&d_charge_grid_complex, Kmax1*Kmax2*Kmax3*sizeof(cufftDoubleComplex)); 

    convert_real_to_complex_kernel<<<(Kmax1*Kmax2*Kmax3+255)/256, 256>>>(Kmax1, Kmax2, Kmax3,
                                                                         d_charge_grid,     
                                                                         d_charge_grid_complex);
    cudaDeviceSynchronize(); 

    cufftExecZ2Z(plan, d_charge_grid_complex, d_charge_grid_complex, CUFFT_FORWARD); 

    cudaDeviceSynchronize();


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
 

    apply_euler_and_interaction_kernel<<<(Kmax1*Kmax2*Kmax3+255)/256, 256>>>(Kmax1, Kmax2, Kmax3,
		                                                             d_charge_grid_complex,
									     d_interaction_mesh,
									     d_Euler_b_1, d_Euler_b_2, d_Euler_b_3);
 
    cudaDeviceSynchronize();

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
  
    cufftExecZ2Z(plan, d_charge_grid_complex, d_charge_grid_complex, CUFFT_INVERSE);

    cudaDeviceSynchronize();



    cublasHandle_t handle;
    cublasCreate(&handle);

    double dot_result;
    cublasDdot(handle, Kmax1*Kmax2*Kmax3, d_charge_grid, 1,                          
           (double*)d_charge_grid_complex, 2,
           &dot_result);               

    cublasDestroy(handle);

    E = dot_result * 0.5 * r4pie;


      
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
 

    cudaMemset(d_FX, 0.e0, npar * sizeof(double));
    cudaMemset(d_FY, 0.e0, npar * sizeof(double));
    cudaMemset(d_FZ, 0.e0, npar * sizeof(double));
 

    p3 = n_order * n_order * n_order;

    compute_forces_kernel<<<npar, p3>>>(npar, n_order,
		                        invlX, invlY, invlZ,
					Kmax1, Kmax2, Kmax3,
					-r4pie,
                                        d_uX, d_uY, d_uZ, d_q,
                                        d_Mn_1, d_Mn_2, d_Mn_3,
                                        d_dMn_1, d_dMn_2, d_dMn_3,
                                        d_charge_grid_complex,
                                        d_FX, d_FY, d_FZ);
    

    // Constrains the net residual force to zero; this operation is optional.

    /*    
    double FXtot = 0.e0, FYtot = 0.e0, FZtot = 0.e0;
  
    for (int i= 0; i < npar; ++i ){
  
        FXtot = FXtot + FX[i]; 
        FYtot = FYtot + FY[i];
        FZtot = FZtot + FZ[i];
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

  
    cufftDoubleComplex *d_Euler_b_1, *d_Euler_b_2, *d_Euler_b_3;
    cudaMalloc(&d_Euler_b_1, Kmax1*sizeof(cufftDoubleComplex));
    cudaMalloc(&d_Euler_b_2, Kmax2*sizeof(cufftDoubleComplex));
    cudaMalloc(&d_Euler_b_3, Kmax3*sizeof(cufftDoubleComplex));

    calc_Euler_splines_coeffs(n_order,  
                              Kmax1, Kmax2, Kmax3, 
                              d_Euler_b_1, d_Euler_b_2, d_Euler_b_3);

    cufftDoubleComplex *d_interaction_mesh;
    cudaMalloc(&d_interaction_mesh, Kmax1*Kmax2*Kmax3*sizeof(cufftDoubleComplex));


        // Declare a handle for the FFT plan.
    cufftHandle plan; // Declare a handle for the FFT plan.
        // Create a 3D double-precision complex-to-complex FFT plan(Z2Z).
        // Dimensions order: fastest-changing to slowest (nx, ny, nz).
    cufftPlan3d(&plan, Kmax3, Kmax2, Kmax1, CUFFT_Z2Z);

 
    	 // calculacte the FT of interaction mesh 
    calc_interaction_mesh(lX, lY, lZ, 
   		          alpha,
			  Kmax1, Kmax2, Kmax3,
			  plan,
                          d_interaction_mesh);

    std::vector<std::complex<double>> 
         interaction_mesh(Kmax3*Kmax2*Kmax1);

    double *d_q;
    double *d_rX, *d_rY, *d_rZ;

    cudaMalloc(&d_q,  npar*sizeof(double));
    cudaMalloc(&d_rX, npar*sizeof(double));
    cudaMalloc(&d_rY, npar*sizeof(double));
    cudaMalloc(&d_rZ, npar*sizeof(double));

    cudaMemcpy(d_q, q.data(), npar * sizeof(double), cudaMemcpyHostToDevice); 
    cudaMemcpy(d_rX, rX.data(), npar * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rY, rY.data(), npar * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_rZ, rZ.data(), npar * sizeof(double), cudaMemcpyHostToDevice); 

    double recE = 0.e0;
    double *d_recFX, *d_recFY, *d_recFZ;
  
    cudaMalloc(&d_recFX, npar*sizeof(double));
    cudaMalloc(&d_recFY, npar*sizeof(double));
    cudaMalloc(&d_recFZ, npar*sizeof(double));
 
    calc_recSpace_Energy_Force(npar,
                               lX, lY, lZ,
                               r4pie, d_q, 
                               d_rX, d_rY, d_rZ, 
                               alpha,
                               n_order,
                               Kmax1, Kmax2, Kmax3,
                               d_Euler_b_1, d_Euler_b_2, d_Euler_b_3,
			       d_interaction_mesh,
			       plan,
                               recE, 
                               d_recFX, d_recFY, d_recFZ);
   
    cudaDeviceSynchronize();      

    std::vector<double> recFX(npar);
    std::vector<double> recFY(npar); 
    std::vector<double> recFZ(npar);

    cudaMemcpy(recFX.data(), d_recFX, npar * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(recFY.data(), d_recFY, npar * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(recFZ.data(), d_recFZ, npar * sizeof(double), cudaMemcpyDeviceToHost); 

    std::cout.precision(16);
    std::cout << recE << std::endl;


    std::ofstream dat_O("Force.dat");
    dat_O.precision(16);
    dat_O <<"unit : 10J/mol/angstrom "<< std::endl;

        for (int i= 0; i < npar; ++i ){

            dat_O << recFX[i] <<"  "
                  << recFY[i] <<"   "
                  << recFZ[i] << std::endl;
        }

    dat_O.close();



////////////////////////////////////////////////
// energy and force in R-space
///////////////////////////////////////////////
    double pi = acosf(-1.e0);
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

    return 0;
}

