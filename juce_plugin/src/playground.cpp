#include <iostream>
//#include <Eigen/Dense>
#include <highfive/highfive.hpp>
//#include <highfive/H5Easy.hpp>
#include <highfive/eigen.hpp>



// Create a dataset name "dset" of double 4x6
//
int main() {
    HighFive::File file("eigen_matrix.h5", HighFive::File::Truncate);

    // Create a matrix.
    Eigen::MatrixXd A(4, 3);
    // clang-format off
    A <<  1,  2,  3,
            4,  5,  6,
            7,  8,  9,
            10, 11, 12;
    // clang-format on
    //
    std::cout << "A = \n" << A << "\n\n";

    // Write it to the file:
    file.createDataSet("mat", A);

    // ... and read it back as fixed-size and row-major:
    using Matrix43d = Eigen::Matrix<double, 4, 3, Eigen::RowMajor>;
    auto B = file.getDataSet("mat").read<Matrix43d>();

    std::cout << "B = \n" << B << "\n";

    return 0;
}