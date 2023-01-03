#include "mex.h"
#include "matrix.h"
/**---------------------------------------------------------------------
 *
 *     ANIREC
 *
 * Matlab implemenation of the program anirec - synthetic seismogram software
 *
 * Author :
 *        Joseph Byrnes, University of Oregon
 *
 * Original Template-Code Author :
 *        Craig Rasmussen, University of Oregon
 *
 * History :
 *        2012-08-15, Template
 *        2012-10-10  Implementation
 *
 ************************************************************************/

/* C declaration for Fortran function implementation */
void anirec ( double * theta_in, double * phi_in, double * z_in, double * vp_in, double * A_in,
        double * B_in, double * vs_in, double * C_in, double * rho_in,                              // Model Inputs
        int * nlayers, double * cc_in, double * baz_in, double * dt_in, int * phase,       //Input parameters
             double * z_comp, double * r_comp, double * t_comp                                 );   // Output traces

void mexFunction( int nlhs,      mxArray * plhs[],
                 int nrhs, const mxArray * prhs[] )
{
    
    
    int nlayers, phase;                                                                    //input ints
    double cc_in, baz_in, dt_in;                                                    //input doubles
    double *theta_in, *phi_in, *z_in, *vp_in, *A_in, *B_in, *vs_in, *C_in, *rho_in; //input arrays
    double *z_comp, *r_comp, *t_comp;                                               //output 1D arrays
    
    double returnvalue;
    
    int i, mrows, ncols; //internal variables
    
    returnvalue = 5;
    
    /* Error checking
     */
    if (nrhs != 14) {
        mexErrMsgTxt("Anirec:  14 input arguments required");
    }
    
    if (nlhs != 3) {
        mexErrMsgTxt("Anirec:  3 output arguments required");
    }
    
    /* Some inputs must be noncomplex scalars.*/
    
    for(i = 9; i < 14; i++) {
        
        mrows = mxGetM(prhs[i]);
        ncols = mxGetN(prhs[i]);
        
        if(mxIsComplex(prhs[i]) ||
           !(mrows==1 && ncols==1) ) {
            mexErrMsgIdAndTxt( "MATLAB:anirec:inputNotRealScalarDouble",
                              "Input must be a noncomplex scalar double.");
        }
        
    }
        
    theta_in = mxGetPr(prhs[0]);
        
    phi_in = mxGetPr(prhs[1]);
    
    z_in = mxGetPr(prhs[2]);
        
    vp_in = mxGetPr(prhs[3]);
        
    A_in = mxGetPr(prhs[4]);
        
    B_in = mxGetPr(prhs[5]);
        
    vs_in = mxGetPr(prhs[6]);
    
    C_in = mxGetPr(prhs[7]);
    
    rho_in = mxGetPr(prhs[8]);
    
    //double variable assignments
    
    nlayers = (int) mxGetScalar(prhs[9]);
    
    cc_in = (double) mxGetScalar(prhs[10]);
    
    baz_in = (double) mxGetScalar(prhs[11]);
        
    dt_in = (double) mxGetScalar(prhs[12]);
    
    phase = (int) mxGetScalar(prhs[13]);
            
    //output array delcrations
    //sizes are static, taken from EQSTRESSv2.5.f
    
    plhs[0] = mxCreateDoubleMatrix(2001, 1, mxREAL);
    
    plhs[1] = mxCreateDoubleMatrix(2001, 1, mxREAL);
    
    plhs[2] = mxCreateDoubleMatrix(2001, 1, mxREAL);
        
    z_comp = mxGetData(plhs[0]);
    
    r_comp = mxGetData(plhs[1]);
    
    t_comp = mxGetData(plhs[2]);
    
    // call anirec
    
    anirec ( theta_in, phi_in, z_in, vp_in, A_in, B_in, vs_in, C_in, rho_in,                      // Model Inputs
                     (int*)&nlayers, (double*)&cc_in, (double*)&baz_in, (double*)&dt_in, (int*)&phase,   //Input parameters
                     z_comp, r_comp, t_comp                                        );                      // Output arrays
    
    
    //mexPrintf("returnvalue - :%f\n", returnvalue);
    
    //mexPrintf("z_comp values - :%f, %f\n", z_comp[0], z_comp[1000] );
        
    //mexPrintf("Some values from the scalers: %f\n", cc_in);
    
    //mexPrintf("MEX CALL COMPLETED\n");
    
}
