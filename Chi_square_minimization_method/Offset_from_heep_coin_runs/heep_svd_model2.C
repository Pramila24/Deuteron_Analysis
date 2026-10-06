#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

//This code uses equations W, Em and Pmx, Pmz to solve a linear system to determine the
//parameters dEf/Ef and dth_e that minimizes the observed and predicted quantities

//****IMPORTANT*** :: ONLY PARAMETERS a0, a1, a2, a3 --> dEb/Eb,  dEf/Ef,  dth_e, dth_p are considered
//Here,  the 12 equations are treated separately, NOT summed over all runs.

void heep_svd_model2()
{

  ifstream ifile;


  
  int runs[6] = {20841, 20846,20851,20858,20861,20869};

  //Define the variables
  Double_t dtr = TMath::Pi() / 180.;
  Double_t Mp = 0.938272;
  Double_t Eb = 10.542;
  Double_t Ef = 8.5207974;
  Double_t Pp[6] = {3.489213 ,3.137317 , 2.776992 , 2.41229,2.0456974 , 1.66214853};
  Double_t th_e[6] = {14.153, 12.94, 11.700, 10.435,9.125,7.704};
  Double_t th_p[6] = {33.344, 35.75, 38.549, 41.812,45.667,50.500};
  //assume out of plane angles in e- and p are zero. (variations are only due to in-plane quantities)
  Double_t ph_e = 0;  
  Double_t ph_p = 0;

  //Define the measured DATA, SIMC variables
  Double_t Wdata[6] = {9.428e-01, 9.433e-01, 9.458e-01, 9.467e-01, 9.471e-01,9.878e-02};
  Double_t Wdata_err[6] = {2.898e-02,  3.77469e-04, 6.01623e-05, 4.13486e-05,9.467e-01,9.878e-03};

  Double_t Wsimc[6] = {9.439e-01, 9.44460e-01, 9.467e-01, 9.48e-01,9.478e-01,9.48e-01};
  Double_t Wsimc_err[6] = {2.931e-02, 2.226e-02,  1.701e-02, 1.298e-02,1.069e-02,9.364e-03};

  Double_t Emdata[6] = {7.534e-03,8.301e-03, 8.486e-03, 9.119e-03, 1.099e-02, 1.228e-02};
  Double_t Emdata_err[6] = {8.3e-03,7.229e-03,  6.581e-04,  3.65940e-05, 4.14000e-05,7.229e-03,7.077e-03}; 

  Double_t Emsimc[6] = {2.938e-03,2.118e-03, 5.802e-03, 8.726e-03, 9.011e-03,8.143e-03};
  Double_t Emsimc_err[6] = {8.405e-03,6.006e-0, 8.347e-03,  7.673e-03, 6.875e-03,7.229e-03,7.312e-03};
  Double_t PmXdata[6] = {-1.328e-03, -1.178e-03, -7.809e-04, 3.741e-04,-7.7286e-04,-3.674e-04};
  Double_t PmXdata_err[6] = {3.30069e-05, 1.75464e-04, 2.87979e-05, 2.08787e-05};

  Double_t PmXsimc[6] = {1.239e-03, 7.288e-04, -2.885e-04, -9.342e-04,-9.057e-04,-2.969e-04};
  Double_t PmXsimc_err[6] = {1.248e-02, 1.058e-02, 9.125e-03, 7.939e-03,6.938e-03, 5.83e-03};

  Double_t PmZdata[6] = {8.512e-03, 8.358e-03, 7.842e-03, 8.513e-03,1.04e-2, 1.135e-02};
  Double_t PmZdata_err[6] = {1.008e-02, 7.456e-03, 6.776e-03, 6.714e-03,7.248e-03,7.891e-03};

  Double_t PmZsimc[6] = {4.288e-03, 3.195e-03, 5.784e-03, 7.84e-03, 8.166e-03, 7.543e-03};
  Double_t PmZsimc_err[6] = {8.718e-03, 7.291e-03, 8.138e-03, 7.08e-03,6.741e-03,7.807e-03};



  //Define derivatives ( and matrix coefficients)
  Double_t dW_dEb, dW_dEf, dW_dthe, dEm_dEb, dEm_dEf, dEm_dPp, dPmx_dEf, dPmx_dthe, dPmx_dPp, dPmx_dthp, dPmz_dEb, dPmz_dEf, dPmz_dthe, dPmz_dPp, dPmz_dthp;

  //Create Matrix and vectors to store the summed elements
  //General Form: [ m x n ] [n] =  [m],  where [mxn] -> has the coefficients Cij, [n] -> column vector (a1,a2)^Transpose 
  //with param to be minimized and [m] is a column vector of length m that carries the obsereved quantities

  //General Equation:  C * aVec = bVec + eVec

  //STEP1: SOLVE linear least squares via SVD
  
  TMatrixD C(24,4);     //4 quations / run, for 6 elastic runs --> total 24 equations
  TVectorD bVec(24);
  TVectorD bVec_err(24);
  TVectorD bVec_cpy(24);

  TVectorD aVec(4);  //vector to store the parameters (a0, a1, a2, a3)
  TVectorD eVec(24);   //vector to store residuals
  TVectorD eVec_err(24);

  int row1, row2, row3, row4; //keep track of the 24 equations

  TMatrixD N_inv(24,24);   //diagonal elements are 1/sig_meas**2, (inverse of the covariance matrix for measured errors, assumed measured errors are independent)
  
  //Define the parameters to be minimized
  Double_t a0, a1, a2, a3;   //a0 --> dEb/Eb,  a1 --> dEf/Ef,  a2 --> dth_e [rad]  a3 --> dth_p [rad]

  //Define the observed quantities
  Double_t dW_obs;
  Double_t dEm_obs;
  Double_t dPmx_obs;
  Double_t dPmz_obs;

  //Define the observed errors
  Double_t dW_obs_err2;
  Double_t dEm_obs_err2;
  Double_t dPmx_obs_err2;
  Double_t dPmz_obs_err2;

  //Systematic Error from Magnet CYcling PROCEDURE (Estimated to be 0.1 MeV or dEf/Ef = 1e-4 GeV)
  Double_t p_syst = 1e-4;

  //Loop over each run
  for (int i=0; i<3; i++)
    {

      cout << "Analyzing Run " << runs[i] << endl;
      
      //Convert to radians
      th_e[i] = th_e[i] * dtr;
      th_p[i] = th_p[i] * dtr;

      
      //Calculate the observed differences bewteen DATA and SIMC
      dW_obs = Wdata[i] - Wsimc[i];      
      dEm_obs = Emdata[i] - Emsimc[i];
      dPmx_obs = PmXdata[i] - PmXsimc[i];
      dPmz_obs = PmZdata[i] - PmZsimc[i];
      
      //Calculate errors of the observed differences bewteen DATA and SIMC
      dW_obs_err2   = pow(Wdata_err[i],2)   + pow(Wsimc_err[i], 2)   + pow(p_syst*Ef, 2);  
      dEm_obs_err2  = pow(Emdata_err[i],2)  + pow(Emsimc_err[i], 2)  + pow(p_syst*Ef, 2);
      dPmx_obs_err2 = pow(PmXdata_err[i],2) + pow(PmXsimc_err[i], 2) + pow(p_syst*Ef, 2);
      dPmz_obs_err2 = pow(PmZdata_err[i],2) + pow(PmZsimc_err[i], 2) + pow(p_syst*Ef, 2);


      //Define all 4 ith rows per run number
      row1 = 4*(i+1) - 4;
      row2 = 4*(i+1) - 3;
      row3 = 4*(i+1) - 2;
      row4 = 4*(i+1) - 1;
      
      //Calculate the derivatives for each run (and the coefficients)

      //---Equation 1---
      dW_dEb = Ef / Eb;
      C[row1][0] = dW_dEb * Eb;

      dW_dEf = - Eb / Ef;
      C[row1][1] = dW_dEf * Ef;   //the '+' is to increment coefficient (sum over all runs)

      dW_dthe = -2.*Eb*Ef/Mp * sin(th_e[i]/2.) * cos(th_e[i]/2.);
      C[row1][2] = dW_dthe;

      C[row1][3] = 0;   

      //---Equation 2---
      dEm_dEb = 1;
      C[row2][0] = dEm_dEb * Eb;

      dEm_dEf = -1;
      C[row2][1] = dEm_dEf * Ef;

      C[row2][2] = 0.;

      C[row2][3] = 0.;

      //---Equation 3---
      C[row3][0] = 0.;

      dPmx_dEf = -sin(th_e[i])*cos(th_e[i]);
      C[row3][1] = dPmx_dEf * Ef;

      dPmx_dthe = -Ef*cos(ph_e)*cos(th_e[i]);
      C[row3][2] = dPmx_dthe;

      dPmx_dthp = -Pp[i]*cos(th_p[i])*cos(ph_p);
      C[row3][3] = dPmx_dthp ;

      //---Equation 4---
      dPmz_dEb = 1.;
      C[row4][0] = dPmz_dEb *Eb;

      dPmz_dEf = -cos(th_e[i])*cos(ph_e);
      C[row4][1] = dPmz_dEf * Ef;

      dPmz_dthe = Ef*cos(ph_e)*sin(th_e[i]);
      C[row4][2] = dPmz_dthe;

      dPmz_dthp = Pp[i]*sin(th_p[i])*cos(ph_p);
      C[row4][3] = dPmz_dthp;


        //Fill the observed vector
      bVec[row1] = dW_obs;
      bVec[row2] = dEm_obs;
      bVec[row3] = dPmx_obs;
      bVec[row4] = dPmz_obs;
      //Fill observed vector error
      bVec_err[row1] = sqrt(dW_obs_err2);
      bVec_err[row2] = sqrt(dEm_obs_err2);
      bVec_err[row3] = sqrt(dPmx_obs_err2);
      bVec_err[row4] = sqrt(dPmz_obs_err2);

      //Fill copy of the observed vector (to be used to solve linear system)
      bVec_cpy[row1] = dW_obs;
      bVec_cpy[row2] = dEm_obs;
      bVec_cpy[row3] = dPmx_obs;
      bVec_cpy[row4] = dPmz_obs;

      N_inv[row1][row1] = 1./dW_obs_err2;
      N_inv[row2][row2] = 1./dEm_obs_err2;
      N_inv[row3][row3] = 1./dPmx_obs_err2;
      N_inv[row4][row4] = 1./dPmz_obs_err2;

      
      cout << "dW_obs = " << dW_obs << " +/-  " << sqrt(dW_obs_err2) << endl;
      cout << "dEm_obs = " << dEm_obs << " +/- " << sqrt(dEm_obs_err2) <<endl;
      cout << "dPmx_obs = " << dPmx_obs << " +/- " << sqrt(dPmx_obs_err2) <<endl;
      cout << "dPmz_obs = " << dPmz_obs << " +/- " << sqrt(dPmz_obs_err2) <<endl;

    }
  

  cout << "MODEL MATRIX" << endl;
  C.Print();
  cout << "Observed Vector" << endl;
  bVec.Print();
  cout << "Inverse Observed Covariance Matrix, N^{-1}" << endl;
  N_inv.Print();
  //Single - Value Decomposition
  TDecompSVD svd_decomp(C);
  bool ok = svd_decomp.Decompose();
  cout << "decomposition ok? " << ok << endl;
  ok = svd_decomp.Solve(bVec_cpy);  //the solution x = A^-1 *bVec is retured in the bVec as: x[0] = bVec[0] . . .in the first n elements of the original  (m x n)

  cout << "solving linear system ok? " << ok << endl;
  bVec_cpy.Print();
  //Get the Optimized Parameters
  aVec[0] = bVec_cpy[0];   //dEb / Eb
  aVec[1] = bVec_cpy[1];  //dEf / Ef
  aVec[2] = bVec_cpy[2];  //dth_e
  aVec[3] = bVec_cpy[3];  //dth_p



  TMatrixD R1(12,4);
  R1.Zero();

  //  V = (C^{T} * N_inv *C) = C^{T} * R1, where R1=N_inv *C
  //R1(24x4) = N_inv(24x24) * C(24x4)
  for(int i=0; i<24; i++){

    for(int j=0; j<4; j++){
      
      R1[i][j] = 0.;  //reset to sum over new elements
      
      for(int k=0; k<24; k++)
	{
	  R1[i][j] += N_inv[i][k] * C[k][j];
	}

    }

  }

  N_inv.Print();
  C.Print();
  R1.Print();
  
  TMatrixD CT(4,24);
  CT.Transpose(C);
  CT.Print();
    
  TMatrixD V(4,4);
  V.Zero();

  //Get Covariance Matrix: V(4x4) = CT(4x24) * R1(24x4) 
  for(int i=0; i<4; i++){

    for(int j=0; j<4; j++){
      
      V[i][j] = 0.;  //reset to sum over new elements
      
      for(int k=0; k<24; k++)
	{
	  V[i][j] += CT[i][k] * R1[k][j];
	}

    }

  }

  V.Print();

  TMatrixD Vinv(4,4);
  Vinv.Zero();
  Vinv = V.Invert();


  //Get Parameter Errors
  Double_t a0_err, a1_err, a2_err, a3_err, a4_err;
  //Off-Diagonal Covariance Errors
  Double_t a01_err, a02_err, a03_err, a04_err, a12_err, a13_err, a14_err, a23_err, a24_err, a34;
  
  a0_err = sqrt(Vinv[0][0]);    //uncertaintiy in dEb/Eb offset
  a1_err = sqrt(Vinv[1][1]);    //uncertaintiy in dEf/Ef offset
  a2_err = sqrt(Vinv[2][2]);    //uncertaintiy in dth_e offset
  a3_err = sqrt(Vinv[3][3]);    //uncertaintiy in dth_p offset

  //Write residuals to file
  ofstream ofs;
  ofs.open("residuals_model2.txt");
  ofs << "# The dObs_ means differential of observable dW, dEm, dPmx, dPmz " << endl;
  ofs << "#!x[i,1]/    residual[f,2]/    residual_err[f,3]/   dObs_meas[f,4]/   dObs_meas_err[f,5]/    dObs_model[f,6]/    dObs_model_err[f,7]/\n";
  
  Double_t dObs_model, dObs_model_err, residual, residual_err;
  Double_t chi2_init, chi2dof_init; 
  Double_t chi2, chi2dof; 
  chi2=0;
  chi2_init=0;

  for(int i=0; i<12; i++)
    {
      
      //Model (Using otimized parameters as input)
      dObs_model = C[i][0] * aVec[0] + C[i][1] * aVec[1] + C[i][2] * aVec[2] + C[i][3] * aVec[3];
      //Model Error (Include correlated errors)
      dObs_model_err = sqrt( pow(C[i][0] * a0_err, 2) + pow(C[i][1] * a1_err,2) + pow(C[i][2] * a2_err, 2) +  pow(C[i][3] * a3_err, 2)
			     + 2.*C[i][0]*C[i][1]*Vinv[0][1] + 2.*C[i][0]*C[i][2]*Vinv[0][2] + 2.*C[i][0]*C[i][3]*Vinv[0][3]
			     + 2.*C[i][1]*C[i][2]*Vinv[1][2] + 2.*C[i][1]*C[i][3]*Vinv[1][3] + 2.*C[i][2]*C[i][3]*Vinv[2][3]);

      //initial chi2 is obtained by setting the parameters to zero --> eVec = -bVEc
      chi2_init += bVec[i]*bVec[i]*N_inv[i][i];
        
      //Residuals: data - model
      eVec[i] = bVec[i] - dObs_model;

      //Residuals Errors: sqrt(data_rr**2 + model_err**2)
      eVec_err[i] = sqrt( pow(bVec_err[i], 2) + pow(dObs_model_err, 2) );

      chi2 +=  eVec[i]*eVec[i]*N_inv[i][i];   //sum over all chi2's

      ofs << "  " << i << "   " << eVec[i] << "   " << eVec_err[i] << "   " << bVec[i] << "   " << bVec_err[i] << "   " << dObs_model << "   " << dObs_model_err << endl;
      
    }
  ofs.close();

  chi2dof = chi2/20.;
  chi2dof_init = chi2_init/20.;

  //INvert V matrix
  cout << "****COVARIANCE MATRIX*****" << endl;

  Vinv.Print();
  

  cout << "---Optimized Parameters---" << endl;
  cout << "total equations x total runs = 6x4 = 24 observations, # parameters = 4, dof = 20 " << endl;
  cout << "initial chi2 = " << chi2_init << " initial chi2/dof = " << chi2dof_init << endl;
  cout << "chi2 = " << chi2 << "  chi2/dof = " << chi2dof << endl;
  cout << "dEb / Eb = " << aVec[0] << endl;
  cout << "dEf / Ef = " << aVec[1] << endl;
  cout << "dth_e = "    << aVec[2] << endl;
  cout << "dth_p = " << aVec[3] << endl;



  cout << "=== Uncertainty in Parameters (Diagonal Elements) ===" << endl;
  cout << "sqrt[cov(0,0)] = dEb / Eb = " << a0_err << endl;
  cout << "sqrt[cov(1,1)] = dEf / Ef = " << a1_err << endl;
  cout << "sqrt[cov(2,2)] = dth_e [rad] = " << a2_err << endl;
  cout << "sqrt[cov(3,3)] = dth_p [rad] = " << a3_err << endl;
  cout << "=== Uncertainty in Parameters (Off-Diagonal Elements) ===" << endl;
  cout << "cov(0,1) = dEb_Eb * dEf_Ef = " << Vinv[0][1] << endl;
  cout << "cov(0,2) = dEb_Eb * dth_e = " << Vinv[0][2] << endl;
  cout << "cov(0,3) = dEb_Eb * dth_p = " << Vinv[0][3] << endl;
  cout << "cov(1,2) = dEf_Ef * dth_e = " << Vinv[1][2] << endl;
  cout << "cov(1,3) = dEf_Ef * dth_p = " << Vinv[1][3] << endl;
  cout << "cov(2,3) = dth_e * dth_p = " << Vinv[2][3] << endl;

  //Get the Correlation Matrix
  //For covariance matrix, Sij, where the diagonals are Si **2, 

  //The correlation matrix is obtained by:  Cm_ij = S_ij / (Si * Sj), where Si, Sj are diagonal terms corresponding to (irow,icol), (jrow,jcol)
  //The diagonals are then:  Cm_ii = S_ii / (Si*Si) = 1 (100% correlation with itself)

  
  TMatrixD Cm(4,4);
  Cm.Zero();

  for(int i=0; i<4; i++){

    for(int j=0; j<4; j++){

      Cm[i][j] =  Vinv[i][j] / (sqrt(Vinv[i][i]) * sqrt(Vinv[j][j]) );

    }

  }


  cout << "****CORRELATION MATRIX*****" << endl;

  Cm.Print();
 

  
}




