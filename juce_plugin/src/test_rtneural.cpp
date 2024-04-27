#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <RTNeural/RTNeural.h>
#include <highfive/H5File.hpp>
#include <Eigen/Dense>
#include "../BatchNorm/BatchNorm1d.h"


template <typename T, typename DenseType>
void loadDense(const HighFive::File& hdf5File, const std::string& layerName, DenseType& dense)
{
    auto weightsDataset = hdf5File.getDataSet(layerName + "/weight");
    auto biasDataset = hdf5File.getDataSet(layerName + "/bias");

    auto weights = weightsDataset.read<std::vector<std::vector<T>>>();
    auto bias = biasDataset.read<std::vector<T>>();

    weightsDataset.read(weights);
    biasDataset.read(bias);

    dense.setWeights(weights);
    dense.setBias(bias.data());
}


template <typename T, typename BatchNormType>
void loadBatchNorm1d(const HighFive::File& hdf5File, const std::string& layerName, BatchNormType& batchNorm1d)
{
    auto runningMeanDataset = hdf5File.getDataSet(layerName + "/running_mean");
    auto runningVarianceDataset = hdf5File.getDataSet(layerName + "/running_var");
    auto epsilonDataset = hdf5File.getDataSet(layerName + "/eps");
    auto weightsDataset = hdf5File.getDataSet(layerName + "/weight");
    auto biasDataset = hdf5File.getDataSet(layerName + "/bias");

    auto runningMean = runningMeanDataset.read<std::vector<T>>();
    auto runningVariance = runningVarianceDataset.read<std::vector<T>>();
    auto epsilon = epsilonDataset.read<T>();
    auto weights = weightsDataset.read<std::vector<T>>();
    auto bias = biasDataset.read<std::vector<T>>();

    runningMeanDataset.read(runningMean);
    runningVarianceDataset.read(runningVariance);
    epsilonDataset.read(epsilon);
    weightsDataset.read(weights);
    biasDataset.read(bias);

    batchNorm1d.setRunningMean(runningMean);
    batchNorm1d.setRunningVariance(runningVariance);
    batchNorm1d.setEpsilon(epsilon);
    batchNorm1d.setGamma(weights);
    batchNorm1d.setBeta(bias);
}


int main(int argc, char* argv[])
{
    //============================================================================
    auto gn1     = RTNeural::BatchNorm1DT<float, 5>();
    auto linear1 = RTNeural::DenseT<float, 5, 16>();
    auto celu1   = RTNeural::ELuActivationT<float, 16>();
    auto linear2 = RTNeural::DenseT<float, 16, 16>();
    auto celu2   = RTNeural::ELuActivationT<float, 16>();
    auto linear3 = RTNeural::DenseT<float, 16, 2>();

    //============================================================================
    std::string hdf5Path = "/Users/yangshijie/Desktop/fuzz/models/2x16_celu.h5";
    HighFive::File hdf5File = HighFive::File(hdf5Path, HighFive::File::ReadOnly);

    loadBatchNorm1d<float>(hdf5File, "gn1", gn1);
    loadDense<float>(hdf5File, "linear1", linear1);
    loadDense<float>(hdf5File, "linear2", linear2);
    loadDense<float>(hdf5File, "linear3", linear3);


    //============================================================================
    // Test the model
    std::string testDataHDF5Path = "/Users/yangshijie/Desktop/fuzz/models/2x16_test.h5";
    HighFive::File testDataHDF5File(testDataHDF5Path, HighFive::File::ReadOnly);
    auto xDataset = testDataHDF5File.getDataSet("x");
    auto yDataset = testDataHDF5File.getDataSet("y");

    auto x = xDataset.read<std::vector<std::vector<float>>>();
    auto y = yDataset.read<std::vector<std::vector<float>>>();

    for (size_t i = 0; i < x.size(); ++i) {
        gn1.forward(Eigen::Map<Eigen::Matrix<float, 5, 1>> { x[i].data() });
        linear1.forward(gn1.outs);
        celu1.forward(linear1.outs);
        linear2.forward(celu1.outs);
        celu2.forward(linear2.outs);
        linear3.forward(celu2.outs + celu1.outs);

        // print the output and expected output
        std::cout << std::fixed << std::setprecision(10) << linear3.outs(0, 0) << " " << linear3.outs(1, 0) << "\n";
        std::cout << std::fixed << std::setprecision(10) << y[i][0] << " " << y[i][1] << std::endl;
    }

    return 0;
}