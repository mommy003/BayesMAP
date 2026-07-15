#include <Rcpp.h>

using namespace Rcpp;


// [[Rcpp::export]]

IntegerVector updateDelta_cpp(
    IntegerVector Delta,
    IntegerVector delta,
    NumericMatrix L,
    IntegerVector gamma,
    NumericMatrix A,
    NumericMatrix B,
    double alpha1,
    double alpha2,
    double alpha3,
    double mu_pi,
    double mu_Pi,
    NumericVector baseline,
    NumericVector pathways,
    double rho
){

    int m = L.nrow();
    int G = L.ncol();
    int C = A.ncol();


    for(int g = 0; g < G; g++){

        double sum0 = 0.0;
        double sum1 = 0.0;


        for(int j = 0; j < m; j++){

            double eta =
                mu_pi +
                baseline[j];


            for(int gg=0; gg<G; gg++){

                if(gg != g)
                    eta += L(j,gg) * Delta[gg] * alpha1;

            }


            for(int c=0;c<C;c++){

                eta += A(j,c) * gamma[c] * alpha2;

            }


            double p0 =
                R::pnorm(eta,0,1,1,0);


            double eta1 =
                eta + L(j,g)*alpha1;


            double p1 =
                R::pnorm(eta1,0,1,1,0);


            p0 = std::min(std::max(p0,1e-12),1-1e-12);
            p1 = std::min(std::max(p1,1e-12),1-1e-12);


            sum0 += delta[j]*log(p0)
                    +(1-delta[j])*log(1-p0);

            sum1 += delta[j]*log(p1)
                    +(1-delta[j])*log(1-p1);

        }


        double xi =
            mu_Pi +
            pathways[g];


        for(int c=0;c<C;c++){

            xi += B(g,c)*gamma[c]*alpha3;

        }


        double Pi =
            R::pnorm(xi,0,1,1,0);


        Pi = std::min(std::max(Pi,1e-12),1-1e-12);


        sum0 += log(1-Pi);
        sum1 += log(Pi);



        double pr1 =
            1.0/(1.0+exp(sum0-sum1));


        if(!R_finite(pr1))
            pr1=0.5;


        pr1 =
            std::min(std::max(pr1,1e-12),1-1e-12);


        Delta[g] =
            R::rbinom(1,pr1);

    }


    return Delta;

}