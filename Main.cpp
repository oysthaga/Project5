# include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <armadillo>
#include <cassert>
#include <complex>

int ij2k(int i, int j, int N) // Transforms the matrix indices ij of matrix into k index of vector. 
{
    return j*N + i; 
}

arma::sp_cx_mat FillMat(int M, std::complex<double> r, arma::cx_vec v) // Fills one matrix with vector v along diagonal.
{
    arma::sp_cx_mat Mat = arma::sp_cx_mat(pow(M-2, 2), pow(M-2, 2));
    for (int l=0; l<pow(M-2, 2); l++)
    {
        Mat(l,l) = v(l);
    }
    for (int l=M-2; l<pow(M-2, 2); l++)
    {
        Mat(l,l-(M-2)) = r;
        Mat(l-(M-2),l) = r;
    }
    for (int l=1; l<pow(M-2, 2); l++)
    {
        Mat(l,l-1) = r;
        Mat(l-1,l) = r;
    }
    for (int l=1; l<M-2; l++)
    {
        Mat((M-2)*l,(M-2)*l-1) = 0+0j;
        Mat((M-2)*l-1,(M-2)*l) = 0+0j;
    }
    return Mat;
}

// arma::sp_cx_mat FillMat(int M, std::complex<double> r, arma::cx_vec v) // Fills one matrix with vector v along diagonal.
// {
//     arma::sp_cx_mat Mat = arma::sp_cx_mat(pow(M-2, 2), pow(M-2, 2));
//     for (int l=0; l<pow(M-2, 2); l++)
//     {
//         Mat(l,l) = v(l);
//     }    
//     for (int l=3; l<pow(M-2, 2); l++)
//     {
//         Mat(l,l-3) = r;
//         Mat(l-3,l) = r;
//     }    
//     for (int l=1; l<pow(M-2, 2); l++)
//     {
//         Mat(l,l-1) = r;
//         Mat(l-1,l) = r;
//     }
//     for (int l=1; l<pow(M-2, 2)/3; l++)
//     {
//         Mat(3*l,3*l-1) = 0+0j;
//         Mat(3*l-1,3*l) = 0+0j;
//     }
//     return Mat;
// }

void fill(arma::sp_cx_mat &A, arma::sp_cx_mat &B, arma::sp_mat V, int M, double h, double Dt, std::complex<double> r) // Fills two matrices A and B. 
{
    arma::cx_vec v = arma::cx_vec(pow(M-2,2), arma::fill::zeros);
    for (int i=0; i<M-2; i++)
    {
        for (int j=0; j<M-2; j++)
        {
            int k = ij2k(i, j, M-2);
            v(k) = V(i,j);
        }
    }

    arma::cx_vec a = 1. + 4.*r + (1j*Dt/2)*v;
    arma::cx_vec b = 1. - 4.*r - (1j*Dt/2)*v;
    A = FillMat(M, -r, a);
    B = FillMat(M, r, b);
}

arma::cx_mat u_initial(double sigx, double sigy, double px, double py, double xc, double yc, int M, double h)
{
    arma::vec x = arma::regspace(0., h, 1.);
    arma::vec y = arma::regspace(0., h, 1.);
    arma::cx_mat uMat = arma::cx_mat(M-2,M-2, arma::fill::zeros);
    std::complex<double> power;
    for (int i=0; i<M-2; i++)
    {
        for (int j=0; j<M-2; j++)
        {
            // x and y is defined at boundary, u is not. Therefore, x_i -> x_{i+1}.
            power = -pow(x(i+1)-xc,2)/(2*pow(sigx,2)) - pow(y(j+1)-yc,2)/(2*pow(sigy,2)) + 1j*px*(x(i+1)-xc) + 1j*py*(y(j+1)-yc);
            uMat(i,j) = exp(power);
        }
    }
    // Normalization: 
    double abs_u = sqrt(arma::accu(pow(arma::real(uMat), 2) + pow(arma::imag(uMat), 2))); 
    uMat = uMat/abs_u;

     
    // arma::cx_vec u = arma::cx_vec(M*M, arma::fill::zeros);
    // for (int i=0; i<M-2; i++)
    // {
    //     for (int j=0; j<M-2; j++)
    //     {
    //         int k = ij2k(i, j, M-2);
    //         u(k) = uMat(i,j);
    //     }
    // }
    return uMat;
}

arma::sp_mat Potential(double v0, int M, int slits=2) // Defines the walls, where the potential is large. 
{
    double thickness = 0.02;
    double centre = 0.5;
    double length = 0.05;

    arma::sp_mat V = arma::sp_mat(M-2,M-2);
    if (slits==2) // double-slit
    {
        for (int i = (centre - thickness/2)*M; i < (centre + thickness/2)*M; i++)
        {
            for (int j = 0; j < (centre - length/2 - length)*M-1; j++) // Wall under first slit
            {
                V(i,j) = v0;
            }
            for (int j = (centre - length/2)*M-1; j < (centre + length/2)*M-1; j++) // Wall between slits
            {
                V(i,j) = v0;
            }
            for (int j = (centre + length/2 + length)*M-1; j < M-2; j++) // Wall over second slit
            {
                V(i,j) = v0;
            }
        }
    }
    if (slits==1) // single-slit
    {
        for (int i = (centre - thickness/2)*M; i < (centre + thickness/2)*M; i++)
        {
            for (int j = 0; j < (centre - length/2)*M-1; j++) // Wall under first slit
            {
                V(i,j) = v0;
            }
            
            for (int j = (centre + length/2)*M-1; j < M-2; j++) // Wall over second slit
            {
                V(i,j) = v0;
            }
        }
    }
    if (slits==3) // Triple-slit
    {
        for (int i = (centre - thickness/2)*M; i < (centre + thickness/2)*M; i++)
        {
            for (int j = 0; j < (centre - length/2 - 2*length)*M-1; j++) // Wall under first slit
            {
                V(i,j) = v0;
            }
            for (int j = (centre - length/2 - length)*M-1; j < (centre - length/2)*M-1; j++) // Wall between slits
            {
                V(i,j) = v0;
            }
            for (int j = (centre + length/2)*M-1; j < (centre + length/2 + length)-1; j++) // Wall over second slit
            {
                V(i,j) = v0;
            }
            for (int j = (centre + length/2 + 2*length)*M-1; j < M-2; j++) // Wall over second slit
            {
                V(i,j) = v0;
            }
        }
    }
    return V;
}


int main(int argc, char* argv[])
{
    // Command line input: 
    double h = atof(argv[1]);
    double Dt = atof(argv[2]);
    double T = atof(argv[3]);
    double x_c = atof(argv[4]);
    double sig_x = atof(argv[5]);
    double p_x = atof(argv[6]);
    double y_c = atof(argv[7]);
    double sig_y = atof(argv[8]);
    double p_y = atof(argv[9]);
    double v0 = atof(argv[10]);
    int slits = atoi(argv[11]);

    
/*
    // ------------
    // Test FillMat
    // ------------
    int N=3;
    arma::cx_mat a = arma::cx_mat(N,N, arma::fill::randu);
    arma::cx_mat b= arma::cx_mat(N,N, arma::fill::randu);

    arma::cx_vec a_vec = arma::cx_vec(N*N, arma::fill::zeros);
    arma::cx_vec b_vec = arma::cx_vec(N*N, arma::fill::zeros);
    
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            int k = ij2k(i, j, N);
            a_vec(k) = a(i,j);
            b_vec(k) = b(i,j);
        }
    }

    std::cout << "\n";
    std::cout << a_vec;
    std::cout << "\n";
    std::cout << b_vec;
    std::cout << "\n";

    arma::sp_cx_mat A = FillMat(5, 1, a_vec);
    arma::sp_cx_mat B = FillMat(5, 1, b_vec);

    std::cout << "\n";
    std::cout << A;
    std::cout << "\n";
    std::cout << B;
    std::cout << "\n";

    // ------------
    // Test Fill
    // ------------
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    int M = 5;
    arma::sp_cx_mat V = arma::sp_cx_mat((M-2), (M-2));
    double h = 1;
    double Dt = 1;
    double r = 1;
    fill(A, B, V, M, h, Dt, r);
*/


    ////---------------
    //// TEST u_initial
    ////---------------
    //arma::cx_mat Utest = u_initial(sig_x, sig_y, p_x, p_y, x_c, y_c, N, h);
    //std::cout << arma::accu(abs(Utest%Utest)); // Outputs one, hence the state is normalized. (Ran with input [1 1 1 1 1 1 1 1 1] and [0.001 0.01 1 0.5 1 1 0.5 1 1]).


    // // ------------
    // // Save to test Potential in Python.
    // // ------------
    // arma::mat V = arma::mat(Potential(1e+10, 100));
    // V.save("Mat1.bin");
    
    int M = 1/h; // Number of x-steps and y-steps. 
    int Nt = T/Dt; // Number of time-steps. 
    std::complex<double> r = 1j*Dt/(2*h*h); 

    arma::sp_mat V = Potential(v0, M, slits); // Potential/walls
    arma::cx_mat U0 = u_initial(sig_x, sig_y, p_x, p_y, x_c, y_c, M, h); // Initial state
    arma::sp_cx_mat A;
    arma::sp_cx_mat B;
    fill(A, B, V, M, h, Dt, r);

    arma::cx_mat U = U0;
    arma::cx_vec u((M-2)*(M-2));
    for (int i=0; i<(M-2); i++)
    {
        for (int j=0; j<(M-2); j++)
        {
            int k = ij2k(i, j, M-2);
            u(k) = U(i,j);
        }
    }

    arma::cx_cube U_evol = arma::cx_cube(M-2, M-2, Nt);
    U_evol.slice(0) = U0; // Initial state
    for (int n=1; n < Nt; n++) // Each time-step
    {
        arma::cx_vec b = B*u; 
        u = arma::spsolve(A, b); // Solve matrix equation
        for (int i=0; i<M-2; i++)
        {
        for (int j=0; j<M-2; j++)
            {
                int k = ij2k(i, j, M-2);
                U_evol.slice(n)(i,j) = u(k);
            }
        }
        
    }
    U_evol.save("Us.bin");


    return 0;
}