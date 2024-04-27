#include <iostream>
#include <Eigen/Dense>

template <typename DType, size_t num_features>
class BatchNorm1d {
public:
    explicit BatchNorm1d()
    : eps(static_cast<DType>(1e-5)) {
        running_mean.setZero();
        running_var.setOnes();
        weight.setOnes();
        bias.setZero();
    }

    ~BatchNorm1d() = default;

    void set_running_mean(const std::vector<DType>& mean) {
        if (mean.size() != num_features) {
            throw std::invalid_argument("Mean vector size must be equal to num_features.");
        }
        for (size_t i = 0; i < num_features; ++i) {
            running_mean(i, 0) = mean[i];
        }
    }

    void set_running_var(const std::vector<DType>& var) {
        if (var.size() != num_features) {
            throw std::invalid_argument("Variance vector size must be equal to num_features.");
        }
        for (size_t i = 0; i < num_features; ++i) {
            running_var(i, 0) = var[i];
        }
    }

    void set_eps(DType epsilon) {
        eps = epsilon;
    }

    void set_weight(const std::vector<DType>& w) {
        if (w.size() != num_features) {
            throw std::invalid_argument("Weight vector size must be equal to num_features.");
        }
        for (size_t i = 0; i < num_features; ++i) {
            weight(i, 0) = w[i];
        }
    }

    void set_bias(const std::vector<DType>& b) {
        if (b.size() != num_features) {
            throw std::invalid_argument("Bias vector size must be equal to num_features.");
        }
        for (size_t i = 0; i < num_features; ++i) {
            bias(i, 0) = b[i];
        }
    }

    Eigen::Matrix<DType, num_features, 1>& forward(const Eigen::Matrix<DType, num_features, 1>& x) {
        // 使用 Eigen::Matrix 对象的 cwiseSqrt() 函数进行开方，并将 eps 作为标量添加到每个元素
        outs = (x - running_mean).cwiseQuotient((running_var.array() + eps).matrix().cwiseSqrt());
        outs = outs.cwiseProduct(weight) + bias;
        return outs;
    }



private:
    Eigen::Matrix<DType, num_features, 1> running_mean;
    Eigen::Matrix<DType, num_features, 1> running_var;
    DType eps;
    Eigen::Matrix<DType, num_features, 1> weight;
    Eigen::Matrix<DType, num_features, 1> bias;

    Eigen::Matrix<DType, num_features, 1> outs;
};
