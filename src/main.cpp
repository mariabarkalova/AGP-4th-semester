#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include "AGPmnog_constrained.h"

using namespace std;

struct ConstrainedTestParams {
    ProblemType type;
    ConstraintProblemType constrType;
    int dim;
    int func_count;
    double epsilon;
    double r;
    int tightness;
    double fraction;
    string name;
    int tau;
};

void runConstrainedTest(const ConstrainedTestParams& params, const string& filename) {
    cout << "\n========================================" << endl;
    cout << "Запуск: " << params.name << endl;
    cout << "Размерность: " << params.dim
        << ", tightness: " << params.tightness
        << ", epsilon: " << params.epsilon
        << ", r: " << params.r << endl;
    cout << "========================================\n" << endl;

    ofstream file(filename);
    if (!file) {
        cerr << "Не удалось открыть файл: " << filename << endl;
        return;
    }

    if (params.dim == 2) {
        file << "func_num;dim;class;constr_type;iterations;"
            << "x1_found;x2_found;f_found;"
            << "x1_true;x2_true;f_true;"
            << "error_f;constr0;constr1" << endl;
    }
    else if (params.dim == 3) {
        file << "func_num;dim;class;constr_type;iterations;"
            << "x1_found;x2_found;x3_found;f_found;"
            << "x1_true;x2_true;x3_true;f_true;"
            << "error_f;constr0;constr1" << endl;
    }
    else {
        file << "func_num;dim;class;constr_type;iterations;"
            << "x1_found;x2_found;x3_found;x4_found;f_found;"
            << "x1_true;x2_true;x3_true;x4_true;f_true;"
            << "error_f;constr0;constr1" << endl;
    }

    double total_iterations = 0;      // сюда копим точки
    double total_point_count = 0;     // сюда копим итерации
    double total_error_f = 0;
    vector<double> errors;
    int feasible_count = 0;
    int max_func = params.func_count;
    const double nan = numeric_limits<double>::quiet_NaN();

    for (int func_num = 1; func_num <= max_func; func_num++) {
        try {
            GlobalSearchConstrained search(params.type, params.constrType, func_num, params.dim,
                params.epsilon, params.r, params.tightness, params.fraction, params.tau);

            search.algorithm();

            vector<double> point, true_point, constraints;
            double value, true_value;
            size_t points_count;

            search.getResults(point, value, points_count, true_point, true_value);
            search.getConstraintsInfo(constraints);

            size_t iter_count = search.getIterationCount();

            bool has_point = (point.size() == (size_t)params.dim);
            bool has_true = (true_point.size() == (size_t)params.dim);

            bool feasible = true;
            for (size_t i = 0; i < constraints.size(); ++i) {
                if (constraints[i] > 1e-6) {
                    feasible = false;
                    break;
                }
            }

            double error_f = (feasible && has_true && std::isfinite(true_value) && has_point)
                ? std::fabs(value - true_value)
                : nan;

            total_iterations += points_count;      // копим точки
            total_point_count += iter_count;        // копим итерации
            if (std::isfinite(error_f)) {
                total_error_f += error_f;
                errors.push_back(error_f);
            }
            if (feasible) feasible_count++;

            string constrTypeStr;
            switch (params.constrType) {
            case CPT_Normal:        
                constrTypeStr = "normal";
                break;
            case CPT_InFeasibleDomain:
                constrTypeStr = "inside"; 
                break;
            case CPT_OutFeasibleDomain:
                constrTypeStr = "outside";
                break;
            case CPT_OnFeasibleBorder:
                constrTypeStr = "border";
                break;
            }

            string className;
            if (params.type == PT_Grishagin) className = "grishagin";
            else if (params.type == PT_GKLS_Simple) className = "simple";
            else className = "hard";

            // Общий префикс строки
            file << func_num << ";" << params.dim << ";" << className << ";" << constrTypeStr << ";"
                << iter_count << ";";

            // Найденная точка
            if (params.dim == 2) {
                if (has_point) file << point[0] << ";" << point[1];
                file << ";";
                if (feasible && std::isfinite(value)) file << value;
                file << ";";
                if (has_true) file << true_point[0] << ";" << true_point[1];
                file << ";";
                if (std::isfinite(true_value)) file << true_value;
                file << ";";
                if (std::isfinite(error_f)) file << error_f;
                file << ";";
                if (constraints.size() > 0) file << constraints[0];
                file << ";";
                if (constraints.size() > 1) file << constraints[1];
                file << endl;
            }
            else if (params.dim == 3) {
                if (has_point) file << point[0] << ";" << point[1] << ";" << point[2];
                file << ";";
                if (feasible && std::isfinite(value)) file << value;
                file << ";";
                if (has_true) file << true_point[0] << ";" << true_point[1] << ";" << true_point[2];
                file << ";";
                if (std::isfinite(true_value)) file << true_value;
                file << ";";
                if (std::isfinite(error_f)) file << error_f;
                file << ";";
                if (constraints.size() > 0) file << constraints[0];
                file << ";";
                if (constraints.size() > 1) file << constraints[1];
                file << endl;
            }
            else {
                if (has_point) {
                    for (int k = 0; k < params.dim; k++) file << point[k] << ";";
                }
                else {
                    for (int k = 0; k < params.dim; k++) file << ";";
                }
                if (feasible && std::isfinite(value)) file << value;
                file << ";";
                if (has_true) {
                    for (int k = 0; k < params.dim; k++) file << true_point[k] << ";";
                }
                else {
                    for (int k = 0; k < params.dim; k++) file << ";";
                }
                if (std::isfinite(true_value)) file << true_value;
                file << ";";
                if (std::isfinite(error_f)) file << error_f;
                file << ";";
                if (constraints.size() > 0) file << constraints[0];
                file << ";";
                if (constraints.size() > 1) file << constraints[1];
                file << endl;
            }

            cout << "  Выполнено " << func_num << "/" << max_func << " функций"
                << " (iters=" << iter_count
                << ", pts=" << points_count
                << ", f=" << value
                << ", feasible=" << feasible << ")" << endl;
        }
        catch (const string& s) {
            cout << "Ошибка для функции " << func_num << ": " << s << endl;
        }
        catch (const std::exception& e) {
            cout << "Ошибка для функции " << func_num << ": " << e.what() << endl;
        }
    }

    file.close();

    if (!errors.empty()) {
        sort(errors.begin(), errors.end());
        cout << "\n========== РЕЗУЛЬТАТЫ ДЛЯ " << params.name << " ==========\n";
        cout << "Допустимых решений: " << feasible_count << " из " << max_func << endl;
        cout << "Среднее количество точек: " << (total_iterations / max_func) << "\n";
        cout << "Среднее количество итераций: " << (total_point_count / max_func) << "\n";
        cout << "Средняя погрешность: "
            << (feasible_count > 0 ? total_error_f / feasible_count : 0.0) << "\n";
        cout << "Результаты сохранены в: " << filename << endl;
    }
}

int main() {
    setlocale(LC_ALL, "Russian");
    vector<ConstrainedTestParams> tests = {
        // Последовательная версия (tau=1)
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_seq", 1},
        //{PT_GKLS_Simple, CPT_Normal, 2, 100, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal_seq", 1},
        //{PT_GKLS_Hard, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_GKLS_Hard_Normal_seq", 1},
        //{PT_GKLS_Simple, CPT_Normal, 3, 20, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal_seq", 1},
        //{PT_GKLS_Simple, CPT_Normal, 4, 20, 0.01, 4.7, 10, 1.0, "4D_GKLS_Simple_Normal_seq", 1},

        // Параллельная версия (tau=4)
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_par2", 2},
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_par4", 4},
        //{PT_GKLS_Simple, CPT_Normal, 2, 100, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal_par4", 4},
        //{PT_GKLS_Simple, CPT_Normal, 2, 100, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal_par8", 8},
        //{PT_GKLS_Hard, CPT_Normal, 2, 100, 0.01, 4.5, 10, 1.0, "2D_GKLS_Hard_Normal_par4", 4},
        //{PT_GKLS_Simple, CPT_Normal, 3, 20, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal_par4", 4},
        {PT_GKLS_Simple, CPT_Normal, 3, 40, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal_par8", 8},
        //{PT_GKLS_Simple, CPT_Normal, 4, 5, 0.01, 4.7, 10, 1.0, "4D_GKLS_Simple_Normal_par8", 8}
    };

    for (const auto& test : tests) {
        string filename = "results_" + test.name + ".csv";
        runConstrainedTest(test, filename);
    }
    return 0;
}