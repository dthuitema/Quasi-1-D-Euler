//This file lists the major function calls required to compute the quasi Euler 1-D in a converging-divering nozzle
#include<Eigen/Dense>


using namespace Eigen;

double maximum(double lambda1, double lambda2, double lambda3);
MatrixXd computeFplus(double eps_stable, VectorXd lambda, MatrixXd Wold, MatrixXd F, int N);
MatrixXd computeRes(MatrixXd Wold, double eps_stable, int N, double gam, VectorXd S_plus, double R_unv, MatrixXd F_plus, MatrixXd Q);
MatrixXd updateW(MatrixXd Wold, VectorXd dt, VectorXd Vol, MatrixXd Res, int N);
MatrixXd updateWin(MatrixXd Wold, VectorXd dt, double a_star, double gam, double gamgam, double P_t, int N, VectorXd c, VectorXd P, double dx, double T_t, double c_v, double R_unv);
MatrixXd updateWout(MatrixXd Wold, VectorXd dt, double dx, int N, VectorXd P, VectorXd c, double c_v, double R_unv, double gam);
MatrixXd computeF(MatrixXd V, int N);
MatrixXd computeQ(VectorXd S_plus, MatrixXd V, int N);
VectorXd computedt(double CFL, MatrixXd V, double dx, int N);
MatrixXd computeState(MatrixXd Win, double gam, int N, double R_unv);
MatrixXd computeWimplicit(VectorXd Vol, double CFL, double dx, VectorXd S_plus, int N, double eps_stable, MatrixXd Wold, double gam, double R_unv, double P_t, double T_t, VectorXd lambda, double max_newtown_err, int max_newton_iter);
MatrixXd computeWimplicit_int(VectorXd Vol, double CFL, double dx, VectorXd S_plus, int N, double eps_stable, MatrixXd Wold, double gam, double R_un, double c_v, double P_t, double T_t, VectorXd lambda, double max_newton_err, int max_newton_iter);
MatrixXd computeJacobian(double eps_stable, VectorXd lambdai, double gam, MatrixXd V, double R_unv, VectorXd S_plus, int N, double c_v);
MatrixXd computeSourceJacobian(double R_unv, VectorXd S_plus, double gam, MatrixXd V, int N);
MatrixXd computeFluxJacobian(double R_un, double gam, MatrixXd V, int N);
MatrixXd computeStateJacobian(double c_v, VectorXd lambda, MatrixXd V, int N);
MatrixXd computeRoeFluxPosJacobian(MatrixXd dF_plus_dW, VectorXd S_plus, int N);
MatrixXd computeRoeFluxNegJacobian(MatrixXd dF_plus_dW, VectorXd S_plus, int N);
VectorXd computeLambda(MatrixXd V, int N);
