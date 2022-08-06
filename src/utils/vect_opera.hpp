
#ifndef VECT_H
#define VECT_H

const double CalcPearsonCorr(const std::vector<float> &x, const std::vector<float> &y);
const double CalcPearsonCorr_old(const std::vector<float> &x, const std::vector<float> &y);
const double CalcPearsonDist(const std::vector<float> &x, const std::vector<float> &y);

const double CalcSpearmanCorr(const std::vector<float> &x, const std::vector<float> &y);
const double CalcSpearmanCorr_old(const std::vector<float> &x, const std::vector<float> &y);
const double CalcSpearmanDist(const std::vector<float> &x, const std::vector<float> &y);

const double CalcMACDist(const std::vector<float> &x, const std::vector<float> &y);


// --- utils for unittests ---

const void getOrder(const std::vector<float> &vec, std::vector<uint> &order);
const void orderToRank(const std::vector<float> &vec, const std::vector<uint> & order, std::vector<float> &rank);

// --- Basic stats ---

const void dual_mean_stddev(const std::vector<float> &vect, const std::vector<size_t> &categ, double &mean1, double &mean2, double &stddev1, double &stddev2, size_t &nb_catg1, size_t &nb_catg2);
const void mean_stddev_min(const std::vector<float> &vect, double &mean, double &stddev, double &min);
const void mean_stddev(const std::vector<float> &vect, double &mean, double &stddev);

#endif