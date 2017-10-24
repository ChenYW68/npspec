#include <Rcpp.h>
#include <math.h>
#include <iostream>


using namespace std;
using namespace Rcpp;


// this does the transformation L^{-1}y
// [[Rcpp::export]]
NumericVector LinvMult(NumericMatrix covarray, NumericVector y,
                       NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    //double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = y.length();
    int dim = locs.ncol();
    NumericVector z(n);
    int nvec[dim];
    nvec[0] = covarray.nrow();
    if( dim > 1 ) nvec[1] = covarray.ncol();
    // dim > 2 not supported yet.

    // number of neighbors + 1
    int m = NNarray.ncol();

    int h[dim];
    double d;
    double ysub[m];
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];


    // should really do i=m separately. Right now, I need to make
    // sure that the (m+1)th row of NNarray is
    // m+1,m,m-1,...,1
    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            if( dim == 3 ) locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covarray(0,0), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                for(el=0; el<dim; el++){
                    h[el] = locsub[k-1][el] - locsub[j-1][el];
                    if(h[el] < 0) h[el] += nvec[el];
                }
                sig[k-1] = covarray( h[0] , h[1] );
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( covarray(0,0) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        d = 0.0;
        if(i==m){
            // fill in first m entries of z
            for(k=1; k<m+1; k++){
                d = 0.0;
                for(el=1; el<k+1; el++){
                    d += Li[k-1][el-1]*ysub[el-1];
                }
                z(NNarray(m-1,m-k)-1) = d;
            }
        } else {
            // fill in ith entry of z
            d = 0.0;
            for(k=1;k<m+1;k++){
                d += Li[m-1][k-1]*ysub[k-1];
            }
            z(i-1) = d;
        }
    }

    return z;
}









// [[Rcpp::export]]
NumericVector LinvTransMult(NumericMatrix covarray, NumericVector z,
                       NumericMatrix locs, IntegerMatrix NNarray) {

    // here, z is the result of Linv * z
    // return x = Linv^T z
    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    //double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = z.length();
    int dim = locs.ncol();
    NumericVector x(n);
    for(j=0;j<n;j++){ x[j] = 0.0; }
    int nvec[dim];
    nvec[0] = covarray.nrow();
    if( dim > 1 ) nvec[1] = covarray.ncol();
    // dim > 2 not supported yet.

    // number of neighbors + 1
    int m = NNarray.ncol();

    int h[dim];
    double d;
    int isub[m];
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];


    // should really do i=m separately. Right now, I need to make
    // sure that the (m+1)th row of NNarray is
    // m+1,m,m-1,...,1
    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            if( dim == 3 ) locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
            isub[m-1-j] = NNarray(i-1,j)-1;
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covarray(0,0), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                for(el=0; el<dim; el++){
                    h[el] = locsub[k-1][el] - locsub[j-1][el];
                    if(h[el] < 0) h[el] += nvec[el];
                }
                sig[k-1] = covarray( h[0] , h[1] );
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( covarray(0,0) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        d = 0.0;
        if(i==m){
            // fill in first m entries of z
            for(k=1; k<m+1; k++){
                d = 0.0;
                for(el=1; el<k+1; el++){
                    x[isub[el-1]] += Li[k-1][el-1]*z[k-1];
                }
            }
        } else {
            // fill in ith entry of z
            for(k=1;k<m+1;k++){
                x[isub[k-1]] += Li[m-1][k-1]*z[i-1];
            }
        }
    }
    return x;
}

// [[Rcpp::export]]
NumericVector vecchiaPrecond(NumericMatrix covarray, NumericVector y,
                              NumericMatrix locs, IntegerMatrix NNarray){

    int n = y.length();
    NumericVector x(n);

    x = LinvMult(covarray, y, locs, NNarray);
    x = LinvTransMult(covarray, x, locs, NNarray);

    return x;
}



// this does the transformation L^{-1}y
// [[Rcpp::export]]
NumericMatrix getLinvEntries(NumericMatrix covarray,
                       NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    //double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = locs.nrow();
    int dim = locs.ncol();
    int nvec[dim];
    nvec[0] = covarray.nrow();
    if( dim > 1 ) nvec[1] = covarray.ncol();
    // dim > 2 not supported yet.

    // number of neighbors + 1
    int m = NNarray.ncol();

    int h[dim];
    double d;
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];

    NumericMatrix LinvEntries(n,m);


    // should really do i=m separately. Right now, I need to make
    // sure that the (m+1)th row of NNarray is
    // m+1,m,m-1,...,1
    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            if( dim == 3 ) locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covarray(0,0), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                for(el=0; el<dim; el++){
                    h[el] = locsub[k-1][el] - locsub[j-1][el];
                    if(h[el] < 0) h[el] += nvec[el];
                }
                sig[k-1] = covarray( h[0] , h[1] );
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( covarray(0,0) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        d = 0.0;
        if(i==m){
            // fill in first m rows of LinvEntries
            for(k=1; k<m+1; k++){
                for(el=1; el<k+1; el++){
                    LinvEntries(k-1,el-1) = Li[k-1][k-el];
                }
            }
        } else {
            // fill in ith row of LinvEntries
            d = 0.0;
            for(k=1;k<m+1;k++){
                LinvEntries(i-1,k-1) = Li[m-1][m-k];
            }
        }
    }

    return LinvEntries;
}




// [[Rcpp::export]]
NumericVector LinvTransMultFromEntries(NumericMatrix LinvEntries, NumericVector z,
                            IntegerMatrix NNarray) {

    // here, z is the result of Linv * z
    // return x = Linv^T z
    //const double PI = 3.141592653589793238463;
    int i;
    int j;

    int n = z.length();
    NumericVector x(n);
    for(j=0;j<n;j++){ x[j] = 0.0; }

    // number of neighbors + 1
    int m = NNarray.ncol();

    // is it this simple (first m rows)
    for(i=1; i<m; i++){
        for(j=1; j<i+1; j++){
            x( NNarray(i-1,j-1)-1 ) += z(i-1)*LinvEntries(i-1,j-1);
        }
    }

    // all rows after m
    for(i=m; i<n+1; i++){
        for(j=1; j<m+1; j++){
            x( NNarray(i-1,j-1)-1 ) += z(i-1)*LinvEntries(i-1,j-1);
        }
    }

    return x;
}

// [[Rcpp::export]]
NumericVector LinvMultFromEntries(NumericMatrix LinvEntries, NumericVector z,
                                       IntegerMatrix NNarray) {

    // here, z is the result of Linv * z
    // return x = Linv^T z
    //const double PI = 3.141592653589793238463;
    int i;
    int j;

    int n = z.length();
    NumericVector x(n);
    for(j=0;j<n;j++){ x[j] = 0.0; }

    // number of neighbors + 1
    int m = NNarray.ncol();

    // is it this simple (first m rows)
    for(i=1; i<m; i++){
        for(j=1; j<i+1; j++){
            x( i - 1 ) += z( NNarray(i-1,j-1) - 1 )*LinvEntries(i-1,j-1);
        }
    }

    // all rows after m
    for(i=m; i<n+1; i++){
        for(j=1; j<m+1; j++){
            x( i-1 ) += z( NNarray(i-1,j-1) - 1 )*LinvEntries(i-1,j-1);
        }
    }

    return x;
}




// [[Rcpp::export]]
NumericVector vecchiaPrecondFromEntries(NumericMatrix LinvEntries, NumericVector y,
                              IntegerMatrix NNarray){

    int n = y.length();
    NumericVector x(n);

    x = LinvMultFromEntries(LinvEntries, y, NNarray);
    x = LinvTransMultFromEntries(LinvEntries, x, NNarray);

    return x;
}




// [[Rcpp::export]]
NumericVector vecchiaLik(NumericMatrix covarray, NumericVector y,
                       NumericMatrix locs, IntegerMatrix NNarray) {

    //const double PI = 3.141592653589793238463;
    int i;
    int j;
    int k;
    int el;
    //double cparms[4] = {covparms[0], covparms[1], covparms[2], covparms[3]};

    int n = y.length();
    int dim = locs.ncol();
    NumericVector ll(1);
    int nvec[dim];
    nvec[0] = covarray.nrow();
    if( dim > 1 ) nvec[1] = covarray.ncol();
    // dim > 2 not supported yet.

    // number of neighbors + 1
    int m = NNarray.ncol();

    int h[dim];
    double d;
    double ysub[m];
    double locsub[m][3];

    double Li[m][m];

    double g[m];
    double sig[m];


    // should really do i=m separately. Right now, I need to make
    // sure that the (m+1)th row of NNarray is
    // m+1,m,m-1,...,1
    for(i=m; i<n+1; i++){

        // first, fill in ysub and locsub in reverse order
        for(j=m-1; j>=0; j--){
            locsub[m-1-j][0] = locs( NNarray(i-1,j)-1, 0 );
            locsub[m-1-j][1] = locs( NNarray(i-1,j)-1, 1 );
            if( dim == 3 ) locsub[m-1-j][2] = locs( NNarray(i-1,j)-1, 2 );
            ysub[m-1-j] = y[ NNarray(i-1,j)-1 ];
        }

        // initialize Li
        for(k=0;k<m;k++){ for(j=0;j<m;j++){ Li[k][j] = 0.0; }}

        // first row has just a single nonzero entry
        Li[1-1][1-1] = pow( covarray(0,0), -0.5 );

        for(j=2; j<m+1; j++){  // j = row of Li

            // initialize g
            for(k=1; k<m+1; k++){
                g[k-1] = 0.0;
            }

            // get first j-1 entries of jth row of L (not Linverse!)
            for(k=1; k<j; k++){

                for(el=0; el<dim; el++){
                    h[el] = locsub[k-1][el] - locsub[j-1][el];
                    if(h[el] < 0) h[el] += nvec[el];
                }
                sig[k-1] = covarray( h[0] , h[1] );
                g[k-1] = 0.0;
                for(el=1; el<k+1; el++){
                    g[k-1] += Li[k-1][el-1]*sig[el-1];
                }
            }

            // get diagonal entry
            d = 0.0;
            for( k=1;k<j;k++ ){ d += g[k-1]*g[k-1]; }
            Li[j-1][j-1] = pow( covarray(0,0) - d, -0.5 );

            // now get first j-1 entries jth row of Linverse
            for(k=1; k<j; k++){
                for(el=1; el<j+1; el++){
                    Li[j-1][k-1] += g[el-1]*Li[el-1][k-1];
                }
                Li[j-1][k-1] = -Li[j-1][k-1]*Li[j-1][j-1];
            }

        }

        d = 0.0;
        if(i==m){
            // fill in first m entries of z
            for(k=1; k<m+1; k++){
                d = 0.0;
                for(el=1; el<k+1; el++){
                    d += Li[k-1][el-1]*ysub[el-1];
                }
                ll(0) += -d*d/2 + log( Li[k-1][k-1] );
            }
        } else {
            // fill in ith entry of z
            d = 0.0;
            for(k=1;k<m+1;k++){
                d += Li[m-1][k-1]*ysub[k-1];
            }
            ll(0) += -d*d/2 + log( Li[m-1][m-1] );
        }
    }

    return ll;
}










