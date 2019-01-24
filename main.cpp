// 2.5 GHz Intel Core i5, clang-800.0.42.1, macOS High Sierra version 10.13.6

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}
// #define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp"

using namespace std;

const double pi = 3.141592654;

/*
    Gauss-Jordan Elimination from Numerical Recipes 3rd ed.
 */
template<class T>
void gaussj(int n, arrayt<T> &a, int m, arrayt<T> &b) {
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv;
    arrayt<T> indxc(n),indxr(n),ipiv(n);
    
    for (j=0;j<n;j++) ipiv(j)=0;
    for (i=0;i<n;i++) {
        big=0.0;
        for (j=0;j<n;j++)
            if (ipiv(j) != 1)
                for (k=0;k<n;k++) {
                    if (ipiv(k) == 0) {
                        if (abs(a(j,k)) >= big) {
                            big=abs(a(j,k));
                            irow=j;
                            icol=k;
                        }
                    }
                }
        ++(ipiv(icol));
        if (irow != icol) {
            for (l=0;l<n;l++) SWAP(a(irow,l),a(icol,l));
            for (l=0;l<m;l++) SWAP(b(irow,l),b(icol,l));
        }
        indxr(i)=irow;
        indxc(i)=icol;
        if (a(icol,icol) == 0.0) throw("gaussj: Singular Matrix");
        pivinv=1.0/a(icol,icol);
        a(icol,icol)=1.0;
        for (l=0;l<n;l++) a(icol,l) *= pivinv;
        for (l=0;l<m;l++) b(icol,l) *= pivinv;
        for (ll=0;ll<n;ll++)
            if (ll != icol) {
                dum=a(ll,icol);
                a(ll,icol)=0.0;
                for (l=0;l<n;l++) a(ll,l) -= a(icol,l)*dum;
                for (l=0;l<m;l++) b(ll,l) -= b(icol,l)*dum;
            }
    }
    for (l=n-1;l>=0;l--) {
        if (indxr(l) != indxc(l))
            for (k=0;k<n;k++)
                SWAP(a(k,indxr(l)),a(k,indxc(l)));
    }
}

/*
 Returns the mth function of our 7 parameter model
 */
double f(int m, double t) {
    switch (m) {
        case (0): return 1;
        case (1): return t;
        case (2): return t*t;
        case (3): return cos((2*pi*t)/12.);
        case (4): return sin((2*pi*t)/12.);
        case (5): return cos((2*pi*t)/6.);
        case (6): return sin((2*pi*t)/6.);
        default: return 0;
    }
}

double sigma_sqr(double y) {
    return 0.002*0.002*y*y;
}

double f_sum(int m, arrayt<double> &x, double t) {
    double sum = 0.;
    for (int i=0;i<m;i++) {
        sum += x(i)*f(i,t);
    }
    return sum;
}

double chi(int N, int m, arrayt<double> &t, arrayt<double> &x, arrayt<double> &y) {
    double sum = 0.;
    for (int i=0;i<N;i++) {
        double fsum = 0.;
        for (int j=0;j<m;j++) {
            fsum += x(j)*f(j,t(i));
        }
        double sigma2 = sigma_sqr(y(i));
        double dx = (y(i)-fsum);
        sum += dx*dx/sigma2;
    }
    sum /= (N-m);
    return sum;
}

double F_lk(arrayt<double> &t, arrayt<double> &y, int N, int l, int k) {
    double sum = 0;
    for (int i=0; i<N; i++) {
        sum += f(l,t(i))*f(k,t(i))/sigma_sqr(y(i));
    }
    return sum;
}

double b_l(arrayt<double> &t, arrayt<double> &y, int N, int l) {
    double sum = 0;
    for (int i=0; i<N; i++) {
        sum += f(l,t(i))*y(i)/sigma_sqr(y(i));
    }
    return sum;
}

int main()
{
    /*
     read in CO2 data file from La Jolla Pier, California from
     http://cdiac.ornl.gov/ftp/trends/co2/ljo.dat
     */
    int i, j, npts, year, t, nval;
    double co2, ymin, ymax;
    // use dynamically sized container classes
    string cline;
    vector<double> x, y;
    ifstream fp; // input file stream
    fp.open( "ljo.dat" );
    if( fp.fail() ) { cout << "Can’t open file."<< endl; exit(0); }
    //--------- read data from file in complicted format -------
    // skip first 16 lines
    for( i=0; i<16; i++) getline( fp, cline ); // read a whole line
    t = 0; // time in months
    npts = 0; // number of data points
    ymin = 1000.0;
    ymax = -ymin;
    for( i=0; i<70; i++) {
        fp >> year;
        if( 0 == i ) nval = 11; else nval = 12; // line is short(?)
        for( j=0; j<nval; j++) {
            fp >> co2;
            if( co2 > 0.0 ) {
                x.push_back( t ); // use auto sizing because we don’t
                y.push_back(co2); // know how many elements there will be
                if( y[npts] > ymax ) ymax = y[npts]; // x and y index like an array
                if( y[npts] < ymin ) ymin = y[npts]; // could use co2 here also
                npts++;
            }
            t += 1;
        }
        if( year >= 2007 ) break; // end of file
        getline( fp, cline ); // read rest of line
    }
    
    int N = npts;
    int M = 7;
    arrayt<double> b(M,1);
    arrayt<double> A(M,M);
    arrayt<double> X(M);
    arrayt<double> time(N);
    arrayt<double> Y(N);
    arrayt<double> least_sqrs(N);
    
    // populate time and Y
    for (int i=0;i<N;i++) {
        time(i) = x[i];
        Y(i) = y[i];
    }
    
    // Compute initial values for A and b
    for (int l = 0; l < M; l++) {
        b(l,0) = b_l(time, Y, N, l);
        for (int k = 0; k < M; k++) {
            A(l,k) = F_lk(time, Y, N, l, k);
        }
    }

    // compute least squares coefficients
    gaussj<double>(M, A, 1, b);
    
    cout << "Coefficients: "<< endl;
    for (int l = 0; l < M; l++) {
        X(l) = b(l,0);
        cout << "x(" << l << ") = " << X(l) << endl;
    }
    
    /*
       Compute our model function
     */
    for (int i=0;i<N;i++) {
        double fsum = 0.;
        for (int j=0;j<M;j++) {
            fsum += X(j)*f(j,time(i));
        }
        least_sqrs(i) = fsum;
    }
    
    // A inverse gives us our squared error terms along the diagonal
    cout << "Error: "<< endl;
    for (int l = 0; l < M; l++) {
        for (int k = 0; k < M; k++) {
            if (l==k) {
                cout << "A(" << l << "," << k << ") = " << sqrt(A(l,k)) << endl;
            }
        }
    }
    
    double chi_sqrd = chi(N, M, time, X, Y);
    cout << "Chi Squared = " << chi_sqrd << endl;
    
    // OUTPUT
    ofstream least_squares;
    ofstream times;
    ofstream original;
    
    if (M == 7) {
        least_squares.open("least_squares_7.dat",ofstream::out);
        if (least_squares.fail() || least_squares.bad()) {
            cout << "Cannot open file" << endl;
            exit(1);
        }
    }
    else if (M == 5) {
        least_squares.open("least_squares_5.dat",ofstream::out);
        if (least_squares.fail() || least_squares.bad()) {
            cout << "Cannot open file" << endl;
            exit(1);
        }
    }
    else if (M == 3) {
        least_squares.open("least_squares_3.dat",ofstream::out);
        if (least_squares.fail() || least_squares.bad()) {
            cout << "Cannot open file" << endl;
            exit(1);
        }
    }
    else {
        least_squares.open("least_squares.dat",ofstream::out);
        if (least_squares.fail() || least_squares.bad()) {
            cout << "Cannot open file" << endl;
            exit(1);
        }
    }
    
    times.open("time.dat",ofstream::out);
    if (times.fail() || times.bad()) {
        cout << "Cannot open file" << endl;
        exit(1);
    }
    
    original.open("original.dat",ofstream::out);
    if (original.fail() || original.bad()) {
        cout << "Cannot open file" << endl;
        exit(1);
    }
    
    for (int i=0;i<N;i++) {
        times << time(i) << endl;
        original << y[i] << endl;
        least_squares << least_sqrs(i) << endl;
    }
    
    times.close();
    original.close();
    least_squares.close();
    return( EXIT_SUCCESS );
}

