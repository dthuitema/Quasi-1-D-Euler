#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<Eigen/Dense>
#include<algorithm>
#include "euler_functions.h" //Header with most function definitions.

using namespace std;
using namespace Eigen;

/*Variable Initialization*/
int i,j,q,N,N_plus;
double t, dx,a,b,N_point;
double gam =1.4, T_t=531.2, P_t=2117.0, R_unv=1716.0;
double t_1=0.8,t_2=3,h=0.15;
double eps_stable, p_exit, Res_max, temp, del_rho, del_p, del_u, c_v, c_p, a_star, gamgam,p_ratio,dpdu, lambda1,lambda2,lambda3,R1,R2,R3,CFL;
vector<double> temp1(5),temp2(5),temp3(5),temp4(5),temp5(6),temp6(5),temp7(6),temp8(5),iter_plot(2),Res_plot(2);
int input_variable(0),max_iter, iter_num;
bool conv_check;
double max_newton_err;
int max_newton_iter, user_input,implicit_state;

/*Function defintions for Data Output*/
void writeFields_Geom(int N_plus){
  int k;
  ofstream write;
  write.open("geom_imp.csv");
  for(k = 0; k < N_plus; k++){
    write << temp1[k] << ", " << temp2[k] << endl;
  }
  write.close();
}
void writeFields_Error(){
  int k;
  ofstream write;
  write.open("error_imp.csv");
  for(k = 0; k < N; k++){
    write << temp1[k] << ", " << temp2[k] << ", " << temp3[k] << ", " << temp4[k] << endl;
  }
  write.close();
}
void writeFields_Flow(){
  int k;
  ofstream write;
  write.open("flow_imp.csv");
  for(k = 0; k < N; k++){
    write << temp1[k] << ", " << temp2[k] << ", " << temp3[k] << ", " << temp4[k] << ", " << temp5[k] << ", " << temp6[k] << endl;
  }
  write.close();
}

int main()
{
  //User input and program initialization
  cout << "Euler 1-D Code" << endl;
  cout << "Please input 1 for an explicit solve, and 2 for implicit: " << endl;
  cin >> user_input;
  if (user_input==1){
    implicit_state = 0;
  } else if (user_input==2){
    implicit_state = 1;
  } else{
    implicit_state=0;
  }
  cout << "Please indicate pressure ratio, exit to inlet: " << endl;
  cin >> p_ratio;
  cout << "Please input the number of grid points: " << endl;
  cin >> N;
  cout << "Please input epsilon for scalar dissipation: " << endl;
  cin >> eps_stable;
  N_plus = N+1;

  //Variable Initialization
  VectorXd x(N), x_plus(N_plus), Vol(N), M(N), P(N), T(N), u(N), rho(N), S(N), S_plus(N_plus), lambda(N), e(N), dt(N), c(N), epsilon(N), Pt(N), Res_Dens(N);
  MatrixXd Wold(3,N), Wnew(3,N), F(3,N), F_plus(3,N), Res(3,N), Q(3,N), Res_max_mat(3,N), V(6,N);

  //Temporary variables for data output and plotting
  temp1.resize(N+1);
  temp2.resize(N+1);
  temp3.resize(N);
  temp4.resize(N);
  temp5.resize(N);
  temp6.resize(N);
  temp7.resize(N);
  temp8.resize(N);

  //Max Iteration Number
  max_iter = 500000; //Maximum number of iterations.
  max_newton_iter = 20; //Maximum number of inner iterations for implicit solution (Newton's Method).
  max_newton_err = pow(10,-15); //Minimum error for Newton's Method.
  CFL = 0.1; //CFL

  double min_Res;
  min_Res = pow(10,-14); //Minimum error for outer (overall) convergence is 1*10^-14.

  //Grid Initialization
  N_point = static_cast<double>(N);
  a = 0.0; b = 1.0; //Nozzle length should be normalized to 0-1 in the x-domain.
  dx = (b-a)/N_point;
  x(0) = a + dx/2;
  x_plus(0) = a; //x_plus and S_plus refer to x-coordinate and surface area at the control surfaces.
  S(0) = 1 - h*pow((sin(M_PI*pow(x(0),t_1))),t_2);
  S_plus(0) = 1 - h*pow((sin(M_PI*pow(x_plus(0),t_1))),t_2);
  for(i = 1; i < N; i++){
    x(i) = x(i-1) + dx;
    x_plus(i) = x_plus(i-1) + dx;
    S(i) = 1 - h*pow((sin(M_PI*pow(x(i),t_1))),t_2);
    S_plus(i) = 1 - h*pow((sin(M_PI*pow(x_plus(i),t_1))),t_2);
  }
  x_plus(N) = b;
  S_plus(N) = 1 - h*pow((sin(M_PI*pow(x_plus(N),t_1))),t_2);
  for(i=0; i<N; i++){
    Vol(i) = dx*(S_plus(i+1)+S_plus(i))/2;
  }

  //Constants
  c_v = R_unv/(gam-1);
  c_p = c_v*gam;
  a_star = 2*gam*((gam-1)/(gam+1))*c_v*T_t;
  gamgam = gam/(gam-1); //Short-hand.

  //State Variable Initialization
  p_exit = p_ratio*P_t; //Program expects an input of a pressure ratio (exit to inlet total pressure.
  M.fill(0.01); //Initial guess for Mach number is 0.01. Can increase for quicker convergence.

  //All other variables are initialized from the Mach number guess and pressure ratio, for a given total inlet pressure and temperature, P_t and T_t.
  for(i = 0; i < N; i++){
    T(i) = pow(((1/T_t)+(gam*R_unv*((gam-1)/(gam+1))*pow(M(i),2)/a_star)),-1);
    c(i) = sqrt(gam*R_unv*T(i));
    u(i) = M(i) * c(i);
    P(i) = P_t*pow(1+((gam-1)/2)*pow(M(i),2),-gamgam);
    rho(i) = P(i)/(R_unv*T(i));
    e(i) = rho(i)*(c_v*T(i) + 0.5*pow(u(i),2));
  }
  P(N-1) = p_exit;
  rho(N-1) = P(N-1)/(R_unv*T(N-1));
  e(N-1) = rho(N-1)*(c_v*T(N-1) + 0.5*pow(u(N-1),2));

  //Solution Vector Initialization
  for(j = 0; j < N; j++){
    Wnew(0,j) = rho(j);
    Wnew(1,j) = rho(j)*u(j);
    Wnew(2,j) = e(j);
  }

  Res.fill(0.0);
  iter_num = 1; //keeps track of the number of iterations
  iter_plot[0] = 1;
  conv_check = 0; //conv_check is used as a Boolean variable for convergence. If the residual is lower than the minimum error, then conv_check will change to 1, and the solution saved/plotted.

  while(conv_check ==0){
  Wold = Wnew;
  V=computeState(Wold,gam,N,R_unv); //Creates a vector to pass state variables to each function.
  F = computeF(V,N); //Computes flux vector
  Q = computeQ(S_plus, V, N); //Computes source vector
  dt = computedt(CFL, V, dx, N); //Computes timestep for each control volume
  lambda = computeLambda(V,N); //Eigen values
  F_plus = computeFplus(eps_stable, lambda, Wold, F, N); //Flux at control surfaces
  Res = computeRes(Wold, eps_stable, N, gam, S_plus, R_unv, F_plus,Q); //Computes the residual for each control volume
  if (implicit_state==0){
    Wnew = updateW(Wold, dt, Vol, Res, N); //Explicit Solution
  } else{
    Wnew = computeWimplicit(Vol, CFL,dx,S_plus,N,eps_stable,Wold,gam,R_unv,P_t,T_t,lambda,max_newton_err,max_newton_iter); //Implicit Solution
  }

  //Update state variables for inner control volumes
  for(i=1;i<(N-1);i++){
  rho(i) = Wnew(0,i);
  u(i) = Wnew(1,i)/Wnew(0,i);
  e(i) = Wnew(2,i);
  epsilon(i) = (e(i)/rho(i)) - 0.5*pow(u(i),2);
  P(i) = (gam-1)*rho(i)*epsilon(i);
  T(i) = P(i)/(rho(i)*R_unv);
  c(i) = sqrt(gam*R_unv*T(i));
  M(i) = u(i)/c(i);
  }

  Wnew = updateWin(Wnew,dt,a_star,gam,gamgam,P_t,N,c,P,dx,T_t,c_v,R_unv); //Update solution vector at nozzle inlet.
  //Update state variables at inlet.
  rho(0) = Wnew(0,0);
  u(0) = Wnew(1,0)/Wnew(0,0);
  e(0) = Wnew(2,0);
  T(0) = T_t*(1-(((gam-1)*(u(0)*u(0)))/((gam+1)*a_star)));
  P(0) = P_t*pow((T(0)/T_t),gamgam);
  rho(0) = P(0)/(R_unv*T(0));
  c(0) = sqrt(gam*R_unv*T(0));
  M(0) = u(0)/c(0);

  Wnew = updateWout(Wnew,dt,dx,N,P,c,c_v,R_unv,gam); //Update solution vector at nozzle outlet.
  //Update state variables at outlet.
  rho(N-1) = Wnew(0,N-1); u(N-1) = Wnew(1,N-1)/Wnew(0,N-1); e(N-1) = Wnew(2,N-1);
  T(N-1) = ((e(N-1)/rho(N-1))-0.5*pow(u(N-1),2))/c_v;
  P(N-1) = T(N-1)*R_unv*rho(N-1);
  c(N-1) = sqrt(gam*R_unv*T(N-1));
  M(N-1) = u(N-1)/c(N-1);

  //Check for convergence
  conv_check = 1;

  for(j = 0; j<3; j++){
  for(i=0; i<N; i++){
  if(Res(j,i) < 0){
    Res_max_mat(j,i) = -1*Res(j,i);
  } else {
    Res_max_mat(j,i) = Res(j,i);
  }
  if (Res_max_mat(0,i) > min_Res){
  conv_check = 0;
  }
  }
  }
  Res_max = Res_max_mat.maxCoeff();
  Res_plot[iter_num-1] = Res_max_mat(0,0);
  for (i=1;i<N;i++){
    if(Res_max_mat(0,i)>Res_plot[iter_num-1]){
      Res_plot[iter_num-1] = Res_max_mat(0,i);
    }
  }
  iter_plot[iter_num-1] = iter_num;
  if(conv_check == 0){
  iter_num = iter_num + 1;
  cout << "iteration: " << iter_num << endl;
  iter_plot.resize(iter_num);
  Res_plot.resize(iter_num);
  }
  else if(conv_check == 1){
  cout << "Convergence achieved in " << iter_num << " iterations." << endl;
  cout << "Maximum residual is: " << Res_max << endl;
  break;
  }
  else{
  cout << "Error: Residual Error." << endl;
  }

  if(iter_num == max_iter){
  cout << "Error: Exceeded Max Iteration Number" << endl;
  cout << "Maximum Residual: " << Res_max << endl;
  break;
  }
  //conv_check = 1;
  }

  Pt = P/P_t; //Variable for plotting pressure distribution normalized by inlet total pressure.


  //Data Output to File: Geometry, Residuals, Flow Variables
  for (j = 0; j < N_plus; j++){
    temp1[j] = x_plus(j);
    temp2[j] = S_plus(j);
  }
  writeFields_Geom(N_plus); //Geometry of nozzle
  temp1.resize(N);
  temp2.resize(N);

  for(j = 0; j < iter_num; j++){
    temp1[j] = iter_plot[j];
    temp2[j] = Res_plot[j];
  }

  writeFields_Error(); //Density residual

  for (j = 0; j < N; j++){
    temp1[j] = x(j);
    temp2[j] = Pt(j);
    temp3[j] = T(j);
    temp4[j] = rho(j);
    temp5[j] = M(j);
    temp6[j] = e(j);
  }
  writeFields_Flow(); //State variables


}
