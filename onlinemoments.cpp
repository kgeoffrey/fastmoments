#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;


class online_moments {
public:
    double x;
    double mean;
    double variance;
    double n;
    double delta;
    double dN;
    double dN2;
    double M2;
    double M3;
    double M4;
    double skewness;
    double kurtosis;

    double term;

    online_moments() {
        mean = 0.0;
        variance = 0.0;
        n = 0.0;
        delta = 0.0;
        dN = 0.0;
        dN2 = 0.0;
        M2 = 0.0;
        M3 = 0.0;
        M4 = 0.0;
        skewness = 0.0;
        kurtosis = 0.0;

        term = 0.0;
    }

    public:
    void push(double dataPoint) {
        x = dataPoint;
        updateMoments();
    }

    void updateMoments() {
        n++;
        delta = x - mean;
        dN = delta / n;
        dN2 = pow(dN, 2);
        mean += dN;

        term = dN2 * delta * (n - 1.0);

        M4 += dN * term * (n*n - 3.0*n + 3.0) + 6.0 * dN2 * M2 - 4.0 * dN * M3;
        M3 += term * (n - 2.0) - 3.0 * dN * M2;
        M2 += delta * (x - mean);

        variance =  M2 / (n - 1);
        skewness = sqrt(n) * M3 / pow(M2, 1.5);
        kurtosis = (n * M4) / (M2 * M2); //3
    }
};

class moving_moments {
public:
    double window;
    double n = 0;

    double d = 0;
    double dk = 0;
    double xd = 0;

    double pre_mean = 0;
    double pre_var = 0;
    double pre_skew = 0;

    double x = 0;
    double x_out = 0;
    double mean = 0;
    double var = 0;
    double skew = 0;
    double kur = 0;
    double variance = 0;
    double kurtosis = 0;
    double skewness = 0;

    std::deque<double> series;

public:
    moving_moments(double Window) {
        window = Window;
    }

    public:
    void push(double dataPoint) {
        x = dataPoint;
        series.push_front(x);

        if (series.size() <= window) {
            n++;
            updateMoments();
        } else {
            x_out = series.back();
            series.pop_back();
            updateMovingMoments();
        }
    }

    void updateMovingMoments() {
        d = x - mean;
        dk = x_out - mean;
        xd = x - x_out;

        x_out = d;
        x = dk;

        pre_mean = xd / window;
        d *= x_out;
        dk *= x;
        pre_var = var + d - dk;
        d *= x_out;
        dk *= x;
        pre_skew = skew + d - dk;
        d *= x_out;
        dk *= x;

        mean += pre_mean;
        xd *= pre_mean;
        var = pre_var - xd;
        xd *= pre_mean;
        skew = pre_skew - 3*pre_mean*(pre_var) + 2*xd;
        xd *= pre_mean;
        kur += d - dk - 4*pre_mean*(pre_skew) + 6*pow(pre_mean,2)*(pre_var) - 3*xd;

        variance = var / (window - 1);
        skewness = sqrt(window) * skew / pow(var, 1.5);
        kurtosis = (window * kur) / (var * var); //3
    }

    void updateMoments() {
        d = x - mean;
        pre_mean = d / n;
        dk = pow(pre_mean, 2);
        mean += pre_mean;
        xd = dk * d * (n - 1.0);

        kur += pre_mean * xd * (n*n - 3.0*n + 3.0) + 6.0 * dk * var - 4.0 * pre_mean * skew;
        skew += xd * (n - 2.0) - 3.0 * pre_mean * var;
        var += d * (x - mean);

        // if window = series.length this should be enabled
        variance = var / (window - 1);
        skewness = sqrt(window) * skew / pow(var, 1.5);
        kurtosis = (window * kur) / (var * var); //3
    }
};


RCPP_EXPOSED_WRAP(online_moments)
RCPP_EXPOSED_AS(online_moments)
RCPP_EXPOSED_WRAP(moving_moments)
RCPP_EXPOSED_AS(moving_moments)

RCPP_MODULE(ParamModule){
using namespace Rcpp;
  class_<online_moments>("online_moments")
  .constructor()
  .method("push", &online_moments::push)
  .field("mean", &online_moments::mean)
  .field("variance", &online_moments::variance)
  .field("skewness", &online_moments::skewness)
  .field("kurtosis", &online_moments::kurtosis)
  ;
  class_<moving_moments>("moving_moments")
  .constructor<double>("Speeding up R")
  .method("push", &moving_moments::push)
  .field("mean", &moving_moments::mean)
  .field("variance", &moving_moments::variance)
  .field("skewness", &moving_moments::skewness)
  .field("kurtosis", &moving_moments::kurtosis)
  .field("series", &moving_moments::series)
  ;
}

// [[Rcpp::export]]
NumericVector get_online_moments(NumericVector data){
    online_moments obj;
    NumericVector moments(4);

    for (auto & element : data) {
        obj.push(element);
    }
    moments[0] = obj.mean;
    moments[1] = obj.variance;
    moments[2] = obj.skewness;
    moments[3] = obj.kurtosis;
    return moments;
}

// [[Rcpp::export]]
NumericVector get_moving_moments(NumericVector data, double window){
    moving_moments obj(window);
    NumericVector moments(4);

    for (auto & element : data) {
        obj.push(element);
    }
    moments[0] = obj.mean;
    moments[1] = obj.variance;
    moments[2] = obj.skewness;
    moments[3] = obj.kurtosis;
    return moments;
}
