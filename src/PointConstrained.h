
#pragma once
class PointConstrained
{
public:
    double x; // координата на кривой Пеано (от 0 до 1)
    double z;  // значение функции в этой точке
    int nu;// индекс: 0 или 1 - ограничения, 2 - критерий

    PointConstrained() : x(0), z(0), nu(2) {}

    PointConstrained(double x_val, double z_val, int nu_val = 2)
        : x(x_val), z(z_val), nu(nu_val) {}

    bool operator==(const PointConstrained& other) const {
        return x == other.x;
    }

    bool operator<(const PointConstrained& other) const {
        return x < other.x;
    }
};