//This subroutine details the functions required to compute the Euler 1D equations for a converging-diverging nozzle.

#include "euler_functions.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<Eigen/Dense>
#include<algorithm>

using namespace std;
using namespace Eigen;

//Finds the maximum of three values (here eigenvectors)
double maximum(double lambda1, double lambda2, double lambda3){
  double max = lambda1;

  if(lambda2 > max){
    max = lambda2;

  }
  if(lambda3 > max){
    max = lambda3;

  }
  return max;
}
//Computes flux
MatrixXd computeF(MatrixXd V, int N){
  int i(0); MatrixXd F(3,N); double rhoi,ui,Pi,ei;

  for (i=0; i < N; i++){
    rhoi = V(0,i); ui = V(1,i); Pi = V(2,i); ei = V(3,i);
    F(0,i) = rhoi*ui;
    F(1,i) = rhoi*(ui*ui)+Pi;
    F(2,i) = (ei+Pi)*ui;
  }
  return F;
}

//Computes the source vector
MatrixXd computeQ(VectorXd S_plus, MatrixXd V, int N){
  MatrixXd Q(3,N);
  int i; double Pi;

  for (i=0; i<N; i++){
    Pi = V(2,i);
    Q(1,i) = Pi*(S_plus(i+1)-S_plus(i));
    Q(0,i) = 0;
    Q(2,i) = 0;
  }
  return Q;
}

//Determines flux at control surfaces for Roe FLux Splitting
MatrixXd computeFplus(double eps_stable,VectorXd lambda, MatrixXd Wold, MatrixXd F, int N){
  int i;
  MatrixXd F_plus(3,N-1);
  for(i=0; i<(N-1);i++){
    F_plus(0,i) = (0.5*(F(0,i)+F(0,i+1))) - (0.5*eps_stable*lambda(i)*(Wold(0,i+1)-Wold(0,i)));
    F_plus(1,i) = (0.5*(F(1,i)+F(1,i+1))) - (0.5*eps_stable*lambda(i)*(Wold(1,i+1)-Wold(1,i)));
    F_plus(2,i) = (0.5*(F(2,i)+F(2,i+1))) - (0.5*eps_stable*lambda(i)*(Wold(2,i+1)-Wold(2,i)));
}
return F_plus;
}

//Determines eigenvectors
VectorXd computeLambda(MatrixXd V, int N){
  int i;
  double ui, ci, ui_1, ci_1, lambda1, lambda2, lambda3;
  VectorXd lambda(N-1);
  for(i=0; i<(N-1); i++){
    ui = V(1,i); ui_1 = V(1,i+1); ci = V(4,i); ci_1 = V(4,i+1);
    lambda1 = (ui+ui_1)/2;
    lambda2 = ((ui+ui_1)/2) + ((ci+ci_1)/2);
    lambda3 = ((ui+ui_1)/2) - ((ci+ci_1)/2);
    if (lambda1 > 0 ){
      lambda(i) = lambda2;
    } else {
      lambda(i) = maximum(lambda1, lambda2, lambda3);
    }
  }
  return lambda;
}

//Finds the residual
MatrixXd computeRes(MatrixXd Wold, double eps_stable,int N, double gam, VectorXd S_plus,double R_unv, MatrixXd F_plus, MatrixXd Q){

  double FR, SR, FL, SL;
  int i,j;
  MatrixXd Res(3,N);
  for(i=0;i<3;i++){
    Res(i,0) = 0;
    Res(i,N-1) = 0;
  }
  for(i=0; i<(N-2);i++){
    SR = S_plus(i+2);
    SL = S_plus(i+1);
    for(j=0;j<3;j++){
      FR = F_plus(j,i+1);
      FL = F_plus(j,i);
      Res(j,i+1) = (FR*SR)-(FL*SL)-Q(j,i+1);
    }
  }
  return Res;
}


//Computes the timestep for each control volume
VectorXd computedt(double CFL, MatrixXd V, double dx, int N){
  VectorXd dt(N); double ui, ci;
  int i;

  for (i=0; i<N; i++){
    ui = V(1,i); ci = V(4,i);
    dt(i) = (CFL*dx)/(ui+ci);
  }
  return dt;
}

//Creates a vector used to pass state variables between functions
MatrixXd computeState(MatrixXd Win, double gam, int N, double R_unv){
  double rhoi,ui,ei,epsiloni,Pi,Ti,ci; int i; MatrixXd V(6,N);
  for (i=0; i<N; i++){
    rhoi = Win(0,i);
    ui = Win(1,i)/Win(0,i);
    ei = Win(2,i);
    epsiloni = (ei/rhoi)-0.5*(ui*ui);
    Pi = (gam-1)*rhoi*epsiloni;
    Ti = Pi/(rhoi*R_unv);
    ci = sqrt(gam*R_unv*Ti);
    V(0,i) = rhoi; V(1,i) = ui; V(2,i) = Pi; V(3,i) = ei; V(4,i) = ci; V(5,i) = Ti;
  }
  return V;
}

//Explicit Solution to solve solution vector
MatrixXd updateW(MatrixXd Wold, VectorXd dt, VectorXd Vol, MatrixXd Res, int N){
  int i,j;
  MatrixXd Wnew(3,N);
  for(j=0;j<3;j++){
    Wnew(j,0) = Wold(j,0);
    Wnew(j,N-1) = Wold(j,N-1);
  }

  for(j=0;j<3;j++){
    for(i=1;i<(N-1);i++){
      Wnew(j,i) = Wold(j,i) - (dt(i)/Vol(i))*Res(j,i);
    }
  }
  return Wnew;
}

//Updates inlet solution vector
MatrixXd updateWin(MatrixXd Wold, VectorXd dt, double a_star, double gam, double gamgam, double P_t, int N, VectorXd c, VectorXd P, double dx, double T_t, double c_v, double R_unv){
  double ui, ci, rhoi, Pi, Pi_1, ci_1, ui_1,Ti,ei;
  double dpdu, del_u, lambdai; MatrixXd Wnew(3,N);
  Wnew = Wold;
  ui= Wold(1,0)/Wold(0,0); ui_1 = Wold(1,1)/Wold(0,1);
  rhoi = Wold(0,0);
  ci = c(0); ci_1 = c(1);
  Pi = P(0); Pi_1 = P(1);
  dpdu = P_t*gamgam*pow(1-((gam-1)/(gam+1))*(pow(ui,2)/a_star),1/(gam-1))*(-2*(gam-1)/(gam+1))*(ui/a_star);
  lambdai = (dt(0)/dx)*(((ui_1+ui)/2)-((ci_1+ci)/2));
  del_u = -lambdai*(Pi_1-Pi-(rhoi*ci)*(ui_1-ui))/(dpdu-rhoi*ci);
  ui = ui + del_u;
  Ti = T_t*(1-(((gam-1)*(ui*ui))/((gam+1)*(a_star))));
  Pi = P_t*pow((Ti/T_t),gamgam);
  rhoi = Pi/(R_unv*Ti);
  ei = rhoi*((c_v*Ti)+(0.5*ui*ui));
  Wnew(0,0) = rhoi;
  Wnew(1,0) = rhoi*ui;
  Wnew(2,0) = ei;

  return Wnew;
}

//Updates outlet solution vector
MatrixXd updateWout(MatrixXd Wold, VectorXd dt, double dx, int N, VectorXd P, VectorXd c, double c_v,double R_unv, double gam){
  double lambda1, lambda2, lambda3, Po, To, rhoo, rhoo_1, uo, co, R1, R2, R3, Mo, eo, co_1, Po_1, uo_1, del_p, del_rho, del_u;
  MatrixXd Wnew(3,N);
  Wnew = Wold;
  rhoo = Wold(0,N-1); rhoo_1 = Wold(0,N-2);
  uo = Wold(1,N-1)/Wold(0,N-1);
  eo = Wold(2,N-1);
  Po = P(N-1); Po_1 = P(N-2);
  co = c(N-1); co_1 = c(N-2);
  uo_1 = Wold(1,N-2)/Wold(0,N-2);

  lambda1 = (dt(N-1)/dx)*((uo+uo_1)/2);
  lambda2 = (dt(N-1)/dx)*(((uo_1+uo)/2)+((co+co_1)/2));
  lambda3 = (dt(N-1)/dx)*(((uo_1+uo)/2)-((co+co_1)/2));
  R1 = -lambda1*(rhoo-rhoo_1-((1/(co*co))*(Po-Po_1)));
  R2 = -lambda2*(Po-Po_1+(rhoo*co*(uo-uo_1)));
  R3 = -lambda3*(Po-Po_1-(rhoo*co*(uo-uo_1)));
  Mo = ((uo+uo_1)/2)/((co+co_1)/2);

  if(Mo > 1){
    del_p = 0.5*(R2+R3);
  } else {
    del_p = 0;
  }
  del_rho = R1 + ((del_p)/(co*co));
  del_u = (R2-del_p)/(rhoo*co);
  rhoo = rhoo + del_rho;
  uo = uo + del_u;
  Po = Po + del_p;
  To = Po/(R_unv*rhoo);
  eo = rhoo*(c_v*To+0.5*pow(uo,2));

  Wnew(0,N-1) = rhoo;
  Wnew(1,N-1) = rhoo*uo;
  Wnew(2,N-1) = eo;

  return Wnew;
}

//Implicit solution
MatrixXd computeWimplicit(VectorXd Vol, double CFL, double dx, VectorXd S_plus, int N, double eps_stable, MatrixXd Wold, double gam, double R_unv, double P_t, double T_t, VectorXd lambda, double max_newton_err, int max_newton_iter){
  MatrixXd Wnew_int(3,N-2),Wnew(3,N); int i,j; double c_v;
  c_v = R_unv/(gam-1);
  Wnew_int = computeWimplicit_int(Vol,CFL,dx,S_plus,N,eps_stable,Wold,gam,R_unv,c_v,P_t,T_t,lambda,max_newton_err,max_newton_iter);
  Wnew = Wold;
  for(i=0;i<3;i++){
    for(j=1;j<N-1;j++){
      Wnew(i,j) = Wnew_int(i,j-1);
    }
  }
  return Wnew;
}

//Function to update internal ("int") control volumes. "int" or "i" refer to internal to the domain, as opposed to control volumes at the boundaries
MatrixXd computeWimplicit_int(VectorXd Vol, double CFL, double dx, VectorXd S_plus, int N, double eps_stable, MatrixXd Wold,double gam, double R_unv,double c_v,double P_t,double T_t,VectorXd lambda, double max_newton_err, int max_newton_iter){
  MatrixXd Wnew_int(3,N-2), V(6,N), dR_dW_int(3,3);
  int i,j;
  MatrixXd Qi(3,N),Fi(3,N),F_plusi(3,N),Res(3,N);
  VectorXd lambdai(N), dti(N),Resi(3*(N-2)), dWi(3*(N-2)), V_dt_old(N-2);
  MatrixXd dW(3,N-2), del_abs(3,N-2), del(3,N-2), dW_old(3,N-2), V_dt(3,N-2), dR_dW(3*(N-2),3*(N-2)), V_dt_mat(3*(N-2),3*(N-2)),A(3*(N-2),3*(N-2));
  double err;
  int iter;
  for(i=0; i<3; i++){
    for(j=0; j<(N-2); j++){
        Wnew_int(i,j) = Wold(i,j+1);
      }
  }
  iter = 1;
  err = 100;
  V = computeState(Wold,gam,N,R_unv);
  dti = computedt(CFL,V,dx,N);

  while (err > max_newton_err && iter < max_newton_iter){
    for (i=0; i<3; i++){
      for (j=1; j<N-1; j++){
        Wold(i,j) = Wnew_int(i,j-1);
      }
    }
    V=computeState(Wold,gam,N,R_unv);
    Qi = computeQ(S_plus,V,N);
    Fi = computeF(V,N);
    F_plusi = computeFplus(eps_stable,lambda,Wold,Fi,N);
    Res = computeRes(Wold, eps_stable, N, gam, S_plus, R_unv, F_plusi, Qi);
    for(i=1;i<(N-1);i++){
      j = 3*i-3;
      Resi(j) = Res(0,i);
      Resi(j+1) = Res(1,i);
      Resi(j+2) = Res(2,i);
    }
    dR_dW = computeJacobian(eps_stable,lambda,gam,V,R_unv,S_plus,N,c_v);
    for(i=1; i<(N-1); i++){
      V_dt_old(i-1) = Vol(i)/dti(i);
    }
    V_dt_mat.fill(0.0);
    for(i=0;i<(N-2);i++){
      j = 3*(i-1)+3;
      V_dt_mat(j,j) = V_dt_old(i);
      V_dt_mat(j+1,j+1) = V_dt_old(i);
    V_dt_mat(j+2,j+2) = V_dt_old(i);
    }
    A = V_dt_mat + dR_dW;
    dWi = A.colPivHouseholderQr().solve(-Resi); //Solution Step.
    for(i=0;i<(N-2);i++){
      j = 3*(i-1)+3;
        dW(0,i) = dWi(j);
      dW(1,i) = dWi(j+1);
      dW(2,i) = dWi(j+2);
    }

  if (iter==1);
  else {
    for(i=0;i<3;i++){
      for(j=0;j<N-2;j++){
        del(i,j) = (dW(i,j)-dW_old(i,j))/dW(i,j);
        if (del(i,j)>=0);
        else{
          del_abs(i,j) = -1*del(i,j);
        }
      }
    }
    err = del.maxCoeff();
  }

  dW_old = dW;

  iter = iter + 1;
  Wnew_int = Wnew_int + dW;
  }


    return Wnew_int;
}


MatrixXd computeJacobian(double eps_stable, VectorXd lambdai, double gam, MatrixXd V, double R_unv, VectorXd S_plus, int N, double c_v){
  MatrixXd dR_dW(3*(N-2),3*(N-2)), dQ_dW(3*(N-2),3*(N-2)),dF_dW(3*(N-1),3*(N-1)), dW_dW(3*(N-1),3*(N-1)), dF_plus_dW(3*(N-1),3*(N-1)), dFR(3*(N-2),3*(N-2)), dFL(3*(N-2),3*(N-2));

  dQ_dW = computeSourceJacobian(R_unv,S_plus,gam,V,N);
  dF_dW = computeFluxJacobian(R_unv,gam,V,N);
  dW_dW = computeStateJacobian(c_v,lambdai,V,N);
  dF_plus_dW = 0.5*(dF_dW-eps_stable*(dW_dW));
  dFR = computeRoeFluxPosJacobian(dF_plus_dW,S_plus,N);
  dFL = computeRoeFluxNegJacobian(dF_plus_dW,S_plus,N);
  dR_dW = dFR - dFL - dQ_dW;
  return dR_dW;
}

MatrixXd computeSourceJacobian(double R_unv, VectorXd S_plus, double gam, MatrixXd V, int N){
  MatrixXd dQ_dW(3*(N-2),3*(N-2)); int i,k; double Ti, ui;

  dQ_dW.fill(0.0);

  for(i =1; i<(N-1); i++){
    k = 3*(i-2)+3;
    ui = V(1,i); Ti = V(5,i);
    dQ_dW(k+1,k) = R_unv*Ti*(S_plus(i+1)-S_plus(i));
    dQ_dW(k+1,k+1) = -(gam-1)*(ui/2)*(S_plus(i+1)-S_plus(i));
    dQ_dW(k+1,k+2) = (gam-1)*(S_plus(i+1)-S_plus(i));
  }

  return dQ_dW;
}

MatrixXd computeFluxJacobian(double R_unv, double gam, MatrixXd V, int N){
  double c_v, ui,Ti,ui_1,Ti_1; MatrixXd dF_dW(3*(N-1),3*(N-1)); int i,k;
  c_v = R_unv/(gam-1);
  dF_dW.fill(0.0);
  for(i=0;i<(N-1);i++){
    k = 3*(i-1)+3;
    ui = V(1,i); Ti = V(5,i); ui_1 = V(1,i+1); Ti_1 = V(5,i+1);
    dF_dW(k,k) = ui + ui_1;
    dF_dW(k+1,k) = (pow(ui,2)+R_unv*Ti)+(pow(ui_1,2)+R_unv*Ti);
    dF_dW(k+2,k) = ui*(c_v*Ti+0.5*pow(ui,2)+R_unv*Ti)+ui_1*(c_v*Ti_1+0.5*pow(ui_1,2)+R_unv*Ti_1);
    dF_dW(k,k+1) = 1 + 1;
    dF_dW(k+1,k+1) = ui + ui_1;
    dF_dW(k+2,k+1) = (0.5*pow(ui,2)+c_v*gam*Ti) + (0.5*pow(ui_1,2)+c_v*gam*Ti_1);
    dF_dW(k,k+2) = (2/ui)+(2/ui_1);
    dF_dW(k+1,k+2) = (gam+1) + (gam+1);
    dF_dW(k+2,k+2) = (ui*gam) + (ui_1*gam);
  }
  return dF_dW;
}

MatrixXd computeStateJacobian(double c_v, VectorXd lambda, MatrixXd V, int N){
  MatrixXd dW_dW(3*(N-1),3*(N-1)); int i,k; double ui, Ti, ui_1, Ti_1;
  dW_dW.fill(0.0);
  for(i=0;i<(N-1);i++){
    k=3*(i-1)+3;
    ui = V(1,i); Ti = V(5,i); ui_1 = V(1,i+1); Ti_1 = V(5,i+1);
    dW_dW(k,k) = (1*lambda(i))-(1*lambda(i));
    dW_dW(k+1,k) = (ui_1-ui)*lambda(i);
    dW_dW(k+2,k) = ((c_v*Ti_1+0.5*pow(ui_1,2))-(c_v*Ti+0.5*pow(ui,2)))*lambda(i);
    dW_dW(k,k+1) = 0;
    dW_dW(k+1,k+1) = 0;
    dW_dW(k+2,k+1) = ((ui_1/2)-(ui/2))*lambda(i);
    dW_dW(k,k+2) = ((1/(c_v*Ti_1+0.5*pow(ui_1,2)))-(1/(c_v*Ti+0.5*pow(ui,2))))*lambda(i);
    dW_dW(k+1,k+2) = ((2/ui_1)-(2/ui))*lambda(i);
    dW_dW(k+2,k+2) = 0;
  }
  return dW_dW;
}

MatrixXd computeRoeFluxPosJacobian(MatrixXd dF_plus_dW, VectorXd S_plus, int N){
  MatrixXd dFR(3*(N-2),3*(N-2)); int i,k,k1;
  dFR.fill(0.0);
  for(i=0;i<(N-2);i++){
    k = 3*(i-1)+3;
    k1 = 3*i+3;
    dFR(k,k) = dF_plus_dW(k1,k1)*S_plus(i+2);
    dFR(k,k+1) = dF_plus_dW(k1,k1+1)*S_plus(i+2);
    dFR(k,k+2) = dF_plus_dW(k1,k1+2)*S_plus(i+2);
    dFR(k+1,k) = dF_plus_dW(k1+1,k1)*S_plus(i+2);
    dFR(k+1,k+1) = dF_plus_dW(k1+1,k1+1)*S_plus(i+2);
    dFR(k+1,k+2) = dF_plus_dW(k1+1,k1+2)*S_plus(i+2);
    dFR(k+2,k) = dF_plus_dW(k1+2,k1)*S_plus(i+2);
    dFR(k+2,k+1) = dF_plus_dW(k1+2,k1+1)*S_plus(i+2);
    dFR(k+2,k+2) = dF_plus_dW(k1+2,k1+2)*S_plus(i+2);
  }

  return dFR;
}

MatrixXd computeRoeFluxNegJacobian(MatrixXd dF_plus_dW, VectorXd S_plus, int N){
  MatrixXd dFL(3*(N-2),3*(N-2)); int i,k;
  dFL.fill(0.0);
  for(i=0; i<(N-2); i++){
    k=3*(i-1)+3;
    dFL(k,k) = dF_plus_dW(k,k)*S_plus(i+1);
    dFL(k,k+1) = dF_plus_dW(k,k+1)*S_plus(i+1);
    dFL(k,k+2) = dF_plus_dW(k,k+2)*S_plus(i+1);
    dFL(k+1,k) = dF_plus_dW(k+1,k)*S_plus(i+1);
    dFL(k+1,k+1) = dF_plus_dW(k+1,k+1)*S_plus(i+1);
    dFL(k+1,k+2) = dF_plus_dW(k+1,k+2)*S_plus(i+1);
    dFL(k+2,k) = dF_plus_dW(k+2,k)*S_plus(i+1);
    dFL(k+2,k+1) = dF_plus_dW(k+2,k+1)*S_plus(i+1);
    dFL(k+2,k+2) = dF_plus_dW(k+2,k+2)*S_plus(i+1);
  }
  return dFL;
}



