#include <ctime>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/src/Core/products/Parallelizer.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <physics.hpp>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <string>
using Eigen::MatrixXcd;
using Eigen::VectorXcd;
static const std::complex<double> I(0,1);
static int D=2;
static int L=15;
static double alpha=2;
static 