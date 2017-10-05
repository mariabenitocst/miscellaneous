/*
 * ============================================================================
 *
 *       Filename:  sampling_biGaussian.cpp
 *
 *    Description:  sampling from a biGaussian distribution (Gaussian
 *                  distribution with a lower and upper standard deviations)
 *
 * ============================================================================
 */


/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <string>
#include <random>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

int main()
{

    double eta = 0.68;
    double sigma_eta_minus = 0.19;
    double sigma_eta_plus  = 0.05;
   

    int nrolls=20000;  // number of experiments

    // Sample uniformly from [0, 1)
    std::default_random_engine generator_uniform;
    std::uniform_real_distribution<double> distribution_uniform(0.0, 1.0);

    // Sample lower Gaussian
    std::default_random_engine generator_minus;
    std::normal_distribution<double> distribution_minus(eta, sigma_eta_minus);

    // Sample upper Gaussian
    std::default_random_engine generator_plus;
    std::normal_distribution<double> distribution_plus(eta, sigma_eta_plus);

    std::ofstream out("eta_bigaussian.dat");

    if(out.is_open())
    {
        for (int i=0; i<nrolls; ++i) 
        {
            double uniform    = distribution_uniform(generator_uniform);
            double sample_eta = 0;
            if(uniform - 0.5 < 0.)
            {
                 sample_eta = 0.70;
                 while(sample_eta >= 0.68)
                 {
                     sample_eta = distribution_minus(generator_minus);
                 }
            }
            else
            {
                 sample_eta = 0.5;
                 while(sample_eta <= 0.68)
                 {
                     sample_eta = distribution_plus(generator_plus);
                 }
            }
            // Write into outputfile
            out << sample_eta << endl;
        }
    }
    out.close();

    // Return
    return 0;
}


double integrand_opt_depth_bulge(double Dl, 
                                 void * parameters)
{
    params_opt p = *(params_opt *) parameters;
    double     z = Dl*sin(p.b);
    double     x = p.R_0-Dl*cos(p.b)*cos(p.l);
    double     y = -Dl*cos(p.b)*sin(p.l);

    // Return
    return BULGE(x, y, z, p.rhoBL, p.eta, p.zeta)*Dl*(1-Dl/p.Ds);
}


double integrate_opt_depth_bulge(double Ds,
                                 double l_0,
                                 double b_0,
                                 double R_0,
                                 double rhoBL,
                                 double eta,
                                 double zeta)
{   
    double sigma = 0.; // this two parameters not use within these bulge functions
    double Rd    = 0.;
    params_opt parameters = {l_0, b_0, Ds, R_0, sigma, Rd, rhoBL, eta, zeta};
    
    gsl_function F; 
    F.function = &integrand_opt_depth_bulge;
    F.params   = &parameters;
    
    double result, error;
    size_t neval;
    
    gsl_integration_cquad_workspace * w = gsl_integration_cquad_workspace_alloc (1000);
    
    gsl_integration_cquad(&F, 0., Ds, 1e-5, 0,w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free (w);
    
    // Return
    return result;
}

double integrand_Ds_bulge (double Ds,
                           void * parameters)
{
    p_integrand_Ds p       = *(p_integrand_Ds *)parameters;
    double         x       = p.R_0-Ds*cos(p.b_0)*cos(p.l_0);
    double         y       = -Ds*cos(p.b_0)*sin(p.l_0);
    double         z       = Ds*sin(p.b_0);
    gsl_interp_accel * acc = gsl_interp_accel_alloc();
    double f=gsl_spline_eval(&p.spline, Ds, acc)*BULGE(x, y, z, p.rhoBL, p.eta, p.zeta)*
             pow(Ds, 2+2*p.beta);
    gsl_interp_accel_free(acc);

    // Return
    return f;
}



double opt_depth_bulge(double l_0,
                       double b_0,
                       double beta,
                       double R_0,
                       double rhoBL,
                       double eta,
                       double zeta)
{
    int          j       = 1000;
    gsl_spline * spline  = gsl_spline_alloc(gsl_interp_cspline, j);
    double       min_Ds  = 0;
    double       max_Ds  = 20;
    double       step_Ds = (max_Ds-min_Ds)/(j-1);
    double integration_opt_depth_Ds0[j];
    double integration_Ds_values[j];

    for (int i=1; i<j; i++)
    {
        integration_Ds_values[i]     = min_Ds+i*step_Ds;
        integration_opt_depth_Ds0[i] = integrate_opt_depth_bulge(min_Ds+i*step_Ds,
                                            l_0, b_0, R_0, rhoBL, eta, zeta);
    }
    gsl_spline_init(spline, integration_Ds_values, integration_opt_depth_Ds0, j);

    p_integrand_Ds p = {*spline, l_0, b_0, beta, R_0};
    gsl_function F;
    F.function       = &integrand_Ds_bulge;
    F.params         = &p;

    double result, error;
    size_t neval;

    gsl_integration_cquad_workspace *w=gsl_integration_cquad_workspace_alloc(1000);
    gsl_integration_cquad(&F, 0., 20, 1e-5, 0, w, &result, &error, &neval);
    gsl_integration_cquad_workspace_free(w);

    // 4*Pi*G/c²=6.01335e-16 kpc/M_sun

    double denominator_=denominator(l_0, b_0, R_0, beta, rhoBL, eta, zeta);
    gsl_spline_free(spline);

    // Return
    return 6.01335e-16*result/denominator_;
}

