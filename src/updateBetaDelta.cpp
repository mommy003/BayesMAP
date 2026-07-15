#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]

List updateBetaDelta_cpp(
    NumericMatrix X,
    NumericVector ycorr,
    NumericVector beta,
    IntegerVector delta,
    NumericVector pi_j,
    NumericVector xpx,
    double vare,
    double sigmaBetaSq
){

    int n = X.nrow();
    int m = X.ncol();

    int nnz = 0;


    for(int j = 0; j < m; j++){

        double old_effect = beta[j] * delta[j];


        // remove previous SNP effect
        if(old_effect != 0){

            for(int i=0;i<n;i++){

                ycorr[i] += X(i,j) * old_effect;

            }
        }


        double rhs = 0.0;


        for(int i=0;i<n;i++){

            rhs += X(i,j) * ycorr[i];

        }


        rhs = rhs / vare;


        double invLhs =
            1.0 /
            (xpx[j]/vare + 1.0/sigmaBetaSq);


        double betaHat =
            invLhs * rhs;


        double logp0 =
            log(1.0-pi_j[j]);


        double logp1 =
            0.5 *
            (log(invLhs)
            - log(sigmaBetaSq)
            + betaHat*rhs)
            + log(pi_j[j]);


        double p_delta1 =
            1.0 /
            (1.0 + exp(logp0-logp1));


        if(!R_finite(p_delta1))
            p_delta1 = 0.5;


        if(p_delta1 < 1e-12)
            p_delta1 = 1e-12;

        if(p_delta1 > 1-1e-12)
            p_delta1 = 1-1e-12;


        delta[j] =
            R::rbinom(1,p_delta1);



        if(delta[j]==1){

            beta[j] =
                R::rnorm(betaHat,
                         sqrt(invLhs));


            for(int i=0;i<n;i++){

                ycorr[i] -= X(i,j)*beta[j];

            }


            nnz++;

        }
        else{

            beta[j]=0;

        }

    }


    return List::create(
        _["beta"]=beta,
        _["delta"]=delta,
        _["ycorr"]=ycorr,
        _["nnz"]=nnz
    );

}