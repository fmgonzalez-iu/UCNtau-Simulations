#include "../inc/fields_nate.h"

void shift(double *x, double *y, double *z, double t, trace* tr) {
    int iLow = (int)(t/SAMPDT);
    double frac = (t - iLow*SAMPDT)/(SAMPDT);
    int iHi = iLow + 1;
    iLow = (iLow % tr->num + tr->num) % tr->num;
    iHi = (iHi % tr->num + tr->num) % tr->num;

    *x = *x + HEATMULT*(tr->x[iLow] + frac*(tr->x[iHi] - tr->x[iLow]));
    *y = *y + HEATMULT*(tr->y[iLow] + frac*(tr->y[iHi] - tr->y[iLow]));
    *z = *z + HEATMULT*(tr->z[iLow] + frac*(tr->z[iHi] - tr->z[iLow]));
}

void force(double *x_in, double *y_in, double *z_in, double *fx, double *fy, double *fz, double *totalU, double* t, trace* tr) //analytical form of halbach field force, mu*del(mod(B))
{
//    printf("%f %f %e\n", *t, *freq, AMPLITUDE * sin(2*M_PI * FREQ * (*t)));
//    printf("%e\n", *x_in);
//    printf("%p\n", totalU);
    double A = 4*B_REM/(M_PI*sqrt(2));

    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;

    double gx=0.0, gy=0.0, gz=0.0, R, r, Rprime, rprime;
    
    R = 0.5 + 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    r = 1.0 - 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    Rprime = 0.5*KAPPA*(1.0 - 1.0/(1 + exp(-KAPPA*x)))*1.0/(1 + exp(-KAPPA*x));
    rprime = -Rprime;

    double rho = sqrt(y*y+z*z);
    double r_zeta = sqrt((rho-R)*(rho-R)+x*x);

    if (z < -1.0 && r_zeta < r)
    {
        double eta = r*atan(x/(sqrt(y*y + z*z) - R));
        double zeta = r - sqrt(x*x + (sqrt(y*y + z*z) - R)*(sqrt(y*y + z*z) - R));
        double sum_cos=0.0, sum_sin=0.0, sum_k_cos=0.0, sum_k_sin=0.0;
        double cos_term=0.0, sin_term=0.0;
        
        double k_n;

        for (int n = 1; n <= N_TERMS; n += 1)
//        for (int n = N_TERMS; n > 0; n -= 1)
        {
            k_n = 2*M_PI*(4.0*n-3.0)/MAG_SPACE;
            
            cos_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*cos(k_n*eta);
            sin_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*sin(k_n*eta);
            
//            printf("%.15e %.15e %.15e\n", k_n, cos_term, sin_term);
            
            sum_cos += cos_term;
            sum_k_cos += k_n*cos_term;
            sum_sin += sin_term;
            sum_k_sin += k_n*sin_term;
        }
        
//        printf("%.15e %.15e %.15e %.15e\n", sum_cos, sum_sin, sum_k_cos, sum_k_sin);
        
        double b_zeta = A*sum_cos;
        double b_eta = A*sum_sin;
        double b_hold = B_HOLD*(r+R)
                            /
                     (sqrt(y*y + z*z));
        
        double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
    
        //zeta, eta
        double d_BZeta[2] = {-1*A*sum_k_cos, -1*A*sum_k_sin};
        double d_BEta[2] = {-1*A*sum_k_sin, A*sum_k_cos};

        //x,y,z
        double d_Bh[3] = {0.0,
                          -B_HOLD*(r+R)*y/(pow(y*y + z*z, 3.0/2.0)),
                          -B_HOLD*(r+R)*z/(pow(y*y + z*z, 3.0/2.0))};
        double d_Zeta[3] = {rprime
                                +
                            (Rprime*(sqrt(y*y + z*z) - R) - x)
                                        /
                            sqrt(x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))),

                            -y*(sqrt(y*y + z*z) - R)
                                /(
                                  sqrt(y*y + z*z)
                                  *sqrt(x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
                                 ),

                            -z*(sqrt(y*y + z*z) - R)
                                /(
                                  sqrt(y*y + z*z)
                                  *sqrt(x*x + ((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
                                 )};
        double d_Eta[3] = {-rprime*atan(
                                x/(R - sqrt(y*y + z*z))
                            )
                                -
                            (r*(R - sqrt(y*y + z*z) - x*Rprime))
                            /
                            (x*x + y*y + z*z + R*R - 2*R*sqrt(y*y + z*z)),

                           -r*x*y
                            /(
                               sqrt(y*y + z*z)
                               *((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))
                               *(1.0 + x*x/((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
                             ),

                           -r*x*z
                            /(
                               sqrt(y*y + z*z)
                               *((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z)))
                               *(1.0 + x*x/((R - sqrt(y*y + z*z))*(R - sqrt(y*y + z*z))))
                             )};
        
        gx = MU_N*(1.0/b_tot)*(
            b_zeta*(d_BZeta[0]*d_Zeta[0] + d_BZeta[1]*d_Eta[0])
            + b_eta*(d_BEta[0]*d_Zeta[0] + d_BEta[1]*d_Eta[0]));
        
        gy = MU_N*(1.0/b_tot)*(
            b_zeta*(d_BZeta[0]*d_Zeta[1] + d_BZeta[1]*d_Eta[1])
            + b_eta*(d_BEta[0]*d_Zeta[1] + d_BEta[1]*d_Eta[1])
            + b_hold*d_Bh[1]);
        
        gz = MU_N*(1.0/b_tot)*(
            b_zeta*(d_BZeta[0]*d_Zeta[2] + d_BZeta[1]*d_Eta[2])
            + b_eta*(d_BEta[0]*d_Zeta[2] + d_BEta[1]*d_Eta[2])
            + b_hold*d_Bh[2]);
        
        *totalU = -MU_N*b_tot + GRAV*MASS_N*z_grav;
        
        gz -= GRAV*MASS_N;
    }
    
    else {
        gx = NAN;
        gy = NAN;
        gz = NAN;
        *totalU = NAN;
    }

    *fx = gx;
    *fy = gy;
    *fz = gz;
}

void fieldstrength(double *x_in, double *y_in, double *z_in, double *totalB, double* t, trace* tr)
{
     double A = 4*B_REM/(M_PI*sqrt(2));

    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;

    double R, r, Rprime, rprime;
    
    R = 0.5 + 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    r = 1.0 - 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    Rprime = 0.5*KAPPA*(1.0 - 1.0/(1 + exp(-KAPPA*x)))*1.0/(1 + exp(-KAPPA*x));
    rprime = -Rprime;

    double rho = sqrt(y*y+z*z);
    double r_zeta = sqrt((rho-R)*(rho-R)+x*x);

    if (z < -1.0 && r_zeta < r)
    {
        double eta = r*atan(x/(sqrt(y*y + z*z) - R));
        double zeta = r - sqrt(x*x + (sqrt(y*y + z*z) - R)*(sqrt(y*y + z*z) - R));
        double sum_cos=0.0, sum_sin=0.0;
        double cos_term=0.0, sin_term=0.0;
        
        double k_n;

        for (int n = 1; n <= N_TERMS; n += 1)
//        for (int n = N_TERMS; n > 0; n -= 1)
        {
            k_n = 2*M_PI*(4.0*n-3.0)/MAG_SPACE;
            
            cos_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*cos(k_n*eta);
            sin_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*sin(k_n*eta);
            
//            printf("%.15e %.15e %.15e\n", k_n, cos_term, sin_term);
            
            sum_cos += cos_term;
            sum_sin += sin_term;
        }
        
//        printf("%.15e %.15e %.15e %.15e\n", sum_cos, sum_sin, sum_k_cos, sum_k_sin);
        
        double b_zeta = A*sum_cos;
        double b_eta = A*sum_sin;
        double b_hold = B_HOLD*(r+R)
                            /
                     (sqrt(y*y + z*z));
        
        double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
    
        *totalB = b_tot;
    }
    
    else {
        *totalB = NAN;
    }
}

void potential(double *x_in, double *y_in, double *z_in, double *totalU, double* t, trace* tr) //-mu*mod(B) + g*z. remember that mu is already negative.
{
    double A = 4*B_REM/(M_PI*sqrt(2));

    double x = *x_in;
    double y = *y_in;
    double z = *z_in;
    shift(&x, &y, &z, *t, tr);
    double z_grav = *z_in;

    double R, r, Rprime, rprime;
    
    R = 0.5 + 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    r = 1.0 - 0.5/(1 + exp(-KAPPA*x));// + AMPLITUDE * sin(2*M_PI * FREQ * (*t));
    Rprime = 0.5*KAPPA*(1.0 - 1.0/(1 + exp(-KAPPA*x)))*1.0/(1 + exp(-KAPPA*x));
    rprime = -Rprime;

    double rho = sqrt(y*y+z*z);
    double r_zeta = sqrt((rho-R)*(rho-R)+x*x);

    if (z < -1.0 && r_zeta < r)
    {
        double eta = r*atan(x/(sqrt(y*y + z*z) - R));
        double zeta = r - sqrt(x*x + (sqrt(y*y + z*z) - R)*(sqrt(y*y + z*z) - R));
        double sum_cos=0.0, sum_sin=0.0;
        double cos_term=0.0, sin_term=0.0;
        
        double k_n;

        for (int n = 1; n <= N_TERMS; n += 1)
//        for (int n = N_TERMS; n > 0; n -= 1)
        {
            k_n = 2*M_PI*(4.0*n-3.0)/MAG_SPACE;
            
            cos_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*cos(k_n*eta);
            sin_term = (n%2 == 0 ? 1 : -1)/(4.0*n-3.0)*(1-exp(-k_n*MAG_THICK))*exp(-k_n*zeta)*sin(k_n*eta);
            
//            printf("%.15e %.15e %.15e\n", k_n, cos_term, sin_term);
            
            sum_cos += cos_term;
            sum_sin += sin_term;
        }
        
//        printf("%.15e %.15e %.15e %.15e\n", sum_cos, sum_sin, sum_k_cos, sum_k_sin);
        
        double b_zeta = A*sum_cos;
        double b_eta = A*sum_sin;
        double b_hold = B_HOLD*(r+R)
                            /
                     (sqrt(y*y + z*z));
        
        double b_tot = sqrt(b_zeta*b_zeta + b_eta*b_eta + b_hold*b_hold);
        
        *totalU = -MU_N*b_tot + GRAV*MASS_N*z_grav;
    }
    
    else {
        *totalU = NAN;
    }
}
