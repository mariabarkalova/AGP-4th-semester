
#pragma once
#include "PointConstrained.h"
#include <cmath>
#include <limits>
#include <stdexcept>

class IntervalConstrained
{
public:
    PointConstrained left;
    PointConstrained right;
    double length;
    double R;
    int index;

    IntervalConstrained()
        : left(), right(), length(0.0), R(0.0), index(0) {}

    IntervalConstrained(const PointConstrained& a, const PointConstrained& b, int idx = 0)
        : left(a), right(b), length(b.x - a.x), R(0.0), index(idx) {}

    void updateInterv() { 
        length = right.x - left.x; 
    }

    double getDelta(int N) const {
        if (length <= 0.0) return 0.0;
        if (N == 1) return length;
        if (N == 2) return std::sqrt(length);
        if (N == 3) return std::cbrt(length);
        if (N == 4) return std::sqrt(std::sqrt(length));
        return std::pow(length, 1.0 / N);
    }

    bool isValid() const {
        return std::isfinite(left.x) && std::isfinite(right.x) && std::nextafter(left.x, right.x) < right.x;
    }

    double getNewPoint(double mu, double r, int N) const {
        if (!isValid()) 
            throw std::runtime_error("Interval cannot be subdivided");
        if (N < 1 || !std::isfinite(mu) || mu <= 0.0 || !std::isfinite(r) || r <= 1.0)
            throw std::invalid_argument("Invalid search parameters");
        //середина интервала
        double x = left.x + 0.5 * length;
        //края одного типа
        if (left.nu == right.nu) {
            const double diff = right.z - left.z;
            x -= std::copysign(0.5 * std::pow(std::abs(diff) / mu, N) / r, diff);
        }
        //защита от выхода за границы
        if (!std::isfinite(x) || x <= left.x || x >= right.x)
            x = left.x + 0.5 * length;
        if (x <= left.x) 
            x = std::nextafter(left.x, right.x);
        if (x >= right.x) 
            x = std::nextafter(right.x, left.x);
        return x;
    }
};
