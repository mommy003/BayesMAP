#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector updateGamma_cpp(
    IntegerVector gamma,
    IntegerVector delta,
    IntegerVector Delta,
    NumericMatrix A,
    NumericMatrix B,
    NumericMatrix L,
    double alpha1,
    double alpha2,
    double alpha3,
    double mu_pi,
    double mu_Pi,
    NumericVector baseline,
    NumericVector pathways,
    double rho
){

    int m = A.nrow();
    int C = A.ncol();
    int G = B.nrow();


    for(int c = 0; c < C; c++){

        double logpost0 = 0.0;
        double logpost1 = 0.0;


        for(int j = 0; j < m; j++){

            double eta0 =
                mu_pi + baseline[j];


            // L Delta contribution
            for(int g=0; g<G; g++){
                eta0 += L(j,g) * Delta[g] * alpha1;
            }


            // A gamma excluding current c
            for(int cc=0; cc<C; cc++){

                if(cc != c)
                    eta0 += A(j,cc) * gamma[cc] * alpha2;

            }


            double psnp0 =
                R::pnorm(eta0,0,1,1,0);


            double eta1 =
                eta0 + A(j,c)*alpha2;


            double psnp1 =
                R::pnorm(eta1,0,1,1,0);


            psnp0 = std::min(std::max(psnp0,1e-12),1-1e-12);
            psnp1 = std::min(std::max(psnp1,1e-12),1-1e-12);


            logpost0 += delta[j]*log(psnp0)
                      +(1-delta[j])*log(1-psnp0);


            logpost1 += delta[j]*log(psnp1)
                      +(1-delta[j])*log(1-psnp1);

        }



        // Gene component

        for(int g=0; g<G; g++){

            double xi0 =
                mu_Pi + pathways[g];


            for(int cc=0;cc<C;cc++){

                if(cc != c)
                    xi0 += B(g,cc)*gamma[cc]*alpha3;

            }


            double pg0 =
                R::pnorm(xi0,0,1,1,0);


            double xi1 =
                xi0 + B(g,c)*alpha3;


            double pg1 =
                R::pnorm(xi1,0,1,1,0);


            pg0 = std::min(std::max(pg0,1e-12),1-1e-12);
            pg1 = std::min(std::max(pg1,1e-12),1-1e-12);


            logpost0 += Delta[g]*log(pg0)
                      +(1-Delta[g])*log(1-pg0);


            logpost1 += Delta[g]*log(pg1)
                      +(1-Delta[g])*log(1-pg1);

        }


        logpost0 += log(1-rho);
        logpost1 += log(rho);



        double pr1 =
            1.0/(1.0+exp(logpost0-logpost1));


        if(!R_finite(pr1))
            pr1 = 0.5;


        pr1 = std::min(std::max(pr1,1e-12),1-1e-12);


        gamma[c] =
            R::rbinom(1,pr1);

    }


    return gamma;

}