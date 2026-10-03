/*
#pragma once
#include "PointConstrained.h"
#include <cmath>

class IntervalConstrained
{
public:
    PointConstrained left;
    PointConstrained right;
    double length;
    double R;
    int index;

    IntervalConstrained()
        : left(), right(), length(1.0), R(0.0), index(0) {}

    IntervalConstrained(const PointConstrained& left_gran, const PointConstrained& right_gran, int idx = 0)
        : left(left_gran), right(right_gran), index(idx)
    {
        updateInterv();
    }

    void updateInterv() {
        length = right.x - left.x;
    }

    double getDelta(int N) const {
        if (length <= 1e-10) return 0;
        return pow(length, 1.0 / N);
    }

    bool isValid() const {
        return length > 1e-12;
    }

    double getNewPoint(double mu, double r, int N, int total_intervals) const
    {
        // Для граничных интервалов - точка в середине
        if (index == 0 || index == total_intervals - 1) {
            return (left.x + right.x) / 2.0;
        }

        // Для внутренних интервалов
        double diff = right.z - left.z;
        double delta_z = fabs(diff);

        if (mu < 1e-10) mu = 1.0;
        if (delta_z < 1e-10) delta_z = 1e-10;

        double sign = (diff > 0) ? 1.0 : -1.0;

        double pow_val = pow(delta_z / mu, N);
        double new_x = (left.x + right.x) / 2.0 - sign * (1.0 / (2.0 * r)) * pow_val;

        if (new_x < left.x) new_x = left.x;
        if (new_x > right.x) new_x = right.x;

        return new_x;
    }
};*/



#pragma once
#include "PointConstrained.h"
#include <cmath>

class IntervalConstrained
{
public:
    PointConstrained left;
    PointConstrained right;
    double length;
    double R;
    int index;

    IntervalConstrained()
        : left(), right(), length(1.0), R(0.0), index(0) {}

    IntervalConstrained(const PointConstrained& left_gran, const PointConstrained& right_gran, int idx = 0)
        : left(left_gran), right(right_gran), index(idx)
    {
        updateInterv();
    }

    void updateInterv() {
        length = right.x - left.x;
    }

    double getDelta(int N) const {
        if (length <= 1e-10) return 0;
        return pow(length, 1.0 / N);
    }

    bool isValid() const {
        return length > 1e-12;
    }

    /*
    double getNewPoint(double mu, double r, int N, int total_intervals) const
    {
        // Для граничных интервалов - точка в середине
        if (index == 0 || index == total_intervals - 1) {
            double mid = (left.x + right.x) / 2.0;
            return mid;
        }

        // Для внутренних интервалов
        double diff = right.z - left.z;
        double delta_z = fabs(diff);

        if (mu < 1e-10) mu = 1.0;
        if (delta_z < 1e-10) delta_z = 1e-10;

        double sign = (diff > 0) ? 1.0 : -1.0;

        double pow_val = pow(delta_z / mu, N);
        double shift = (1.0 / (2.0 * r)) * pow_val;

        //double new_x = (left.x + right.x) / 2.0 - sign * shift;  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! было 
        double new_x = (left.x + right.x) / 2.0 - (right.z - left.z) / (2.0 * r * mu);

        // Не позволяем точке совпасть с границами
        if (new_x <= left.x) new_x = left.x + (right.x - left.x) * 0.001;
        if (new_x >= right.x) new_x = right.x - (right.x - left.x) * 0.001;

        return new_x;
    }*/
    double getNewPoint(double mu, double r, int N, int total_intervals) const
    {
        // Случай 1: индексы разные — середина интервала
        if (left.nu != right.nu) {
            return (left.x + right.x) / 2.0;
        }

        // Случай 2: индексы одинаковые — каноническая формула
        if (mu < 1e-10) mu = 1.0;
        double mid = (left.x + right.x) / 2.0;
        double new_x = mid - (right.z - left.z) / (2.0 * r * mu);

        // Точка строго внутри интервала
        double eps = (right.x - left.x) * 1e-9;
        if (new_x <= left.x)  new_x = left.x + eps;
        if (new_x >= right.x) new_x = right.x - eps;

        return new_x;
    }
};