
#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
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

    for (int func_num = 1; func_num <= max_func; func_num++) {
        try {
            GlobalSearchConstrained search(params.type, params.constrType, func_num, params.dim, params.epsilon, 
                params.r, params.tightness, params.fraction, params.tau);

            search.algorithm();

            vector<double> point;
            double value;
            size_t points_count;
            vector<double> true_point;
            double true_value;

            search.getResults(point, value, points_count, true_point, true_value);

            vector<double> constraints;
            search.getConstraintsInfo(constraints);

            size_t iter_count = search.getIterationCount();

            double error_f = fabs(value - true_value);
            total_iterations += points_count;      // копим точки
            total_point_count += iter_count;       // копим итерации
            total_error_f += error_f;
            errors.push_back(error_f);

            bool feasible = true;
            for (double c : constraints) {
                if (c > 1e-6) {
                    feasible = false;
                    break;
                }
            }
            if (feasible) feasible_count++;

            string constrTypeStr;
            switch (params.constrType) {
            case CPT_Normal:           constrTypeStr = "normal"; break;
            case CPT_InFeasibleDomain: constrTypeStr = "inside"; break;
            case CPT_OutFeasibleDomain: constrTypeStr = "outside"; break;
            case CPT_OnFeasibleBorder: constrTypeStr = "border"; break;
            }

            string className;
            if (params.type == PT_Grishagin) className = "grishagin";
            else if (params.type == PT_GKLS_Simple) className = "simple";
            else className = "hard";

            if (params.dim == 2) {
                file << func_num << ";" << params.dim << ";" << className << ";" << constrTypeStr << ";"
                    << points_count << ";" << point[0] << ";" << point[1] << ";" << value << ";"
                    << true_point[0] << ";" << true_point[1] << ";" << true_value << ";"
                    << error_f << ";" << constraints[0] << ";" << constraints[1] << endl;
            }
            else if (params.dim == 3) {
                file << func_num << ";" << params.dim << ";" << className << ";" << constrTypeStr << ";"
                    << points_count << ";" << point[0] << ";" << point[1] << ";" << point[2] << ";" << value << ";"
                    << true_point[0] << ";" << true_point[1] << ";" << true_point[2] << ";" << true_value << ";"
                    << error_f << ";" << constraints[0] << ";" << constraints[1] << endl;
            }
            else {
                file << func_num << ";" << params.dim << ";" << className << ";" << constrTypeStr << ";"
                    << points_count << ";" << point[0] << ";" << point[1] << ";" << point[2] << ";" << point[3] << ";" << value << ";"
                    << true_point[0] << ";" << true_point[1] << ";" << true_point[2] << ";" << true_point[3] << ";" << true_value << ";"
                    << error_f << ";" << constraints[0] << ";" << constraints[1] << endl;
            }

            //if (func_num % 10 == 0) {
              //  cout << "  Выполнено " << func_num << "/" << max_func << " функций" << endl;
            //}
            cout << "  Выполнено " << func_num << "/" << max_func << " функций" << endl;
        }
        catch (string s) {
            cout << "Ошибка для функции " << func_num << ": " << s << endl;
        }
    }

    file.close();

    if (!errors.empty()) {
        sort(errors.begin(), errors.end());
        cout << "\n========== РЕЗУЛЬТАТЫ ДЛЯ " << params.name << " ==========\n";
        //cout << "Допустимых решений: " << feasible_count << " из " << max_func << " (" << (100.0 * feasible_count / max_func) << "%)" << endl;
        cout << "Среднее количество точек: " << (total_iterations / max_func) << "\n";
        cout << "Среднее количество итераций: " << (total_point_count / max_func) << "\n";
        cout << "Средняя погрешность: " << (total_error_f / max_func) << "\n";
        cout << "Результаты сохранены в: " << filename << endl;
    }
}

int main() {
    setlocale(LC_ALL, "Russian");
    vector<ConstrainedTestParams> tests = {

        //{PT_Grishagin, CPT_InFeasibleDomain, 2, 5, 0.001, 6.0, 10, 0.3, "2D_Grishagin_Inside"},
        //{PT_Grishagin, CPT_OutFeasibleDomain, 2, 5, 0.001, 6.0, 10, 0.3, "2D_Grishagin_Outside" },
        // 2D GKLS Simple
        /*
        {PT_GKLS_Simple, CPT_Normal, 2, 5, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal"},
        {PT_GKLS_Simple, CPT_InFeasibleDomain, 2, 5, 0.001, 4.5, 10, 0.5, "2D_GKLS_Simple_Inside"},
        {PT_GKLS_Simple, CPT_OutFeasibleDomain, 2, 5, 0.001, 4.5, 10, 0.5, "2D_GKLS_Simple_Outside"},

        // 2D GKLS Hard
        {PT_GKLS_Hard, CPT_Normal, 2, 5, 0.01, 3.0, 10, 1.0, "2D_GKLS_Hard_Normal"},
        {PT_GKLS_Hard, CPT_InFeasibleDomain, 2, 5, 0.001, 6.0, 10, 0.5, "2D_GKLS_Hard_Inside"},
        {PT_GKLS_Hard, CPT_OutFeasibleDomain, 2, 5, 0.001, 6.0, 10, 0.5, "2D_GKLS_Hard_Outside"},

        // 3D GKLS Simple
        {PT_GKLS_Simple, CPT_InFeasibleDomain, 3, 5, 0.1, 8.0, 5, 0.3, "3D_GKLS_Simple_Inside"},
        {PT_GKLS_Simple, CPT_OutFeasibleDomain, 3, 5, 0.01, 4.0, 8, 0.3, "3D_GKLS_Simple_Outside"},

        // 3D GKLS Hard
        {PT_GKLS_Hard, CPT_InFeasibleDomain, 3, 10, 0.01, 10.0, 10, 0.3, "3D_GKLS_Hard_Inside"},
        {PT_GKLS_Hard, CPT_OutFeasibleDomain, 3, 10, 0.01, 10.0, 10, 0.3, "3D_GKLS_Hard_Outside"},
        */


        //{PT_Grishagin, CPT_Normal, 2, 5, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal"}, // да
        //{PT_GKLS_Simple, CPT_Normal, 2, 5, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal"}, //да
        //{PT_GKLS_Hard, CPT_Normal, 2, 3, 0.01, 4.5, 10, 1.0, "2D_GKLS_Hard_Normal"},   //большая погрешность
        //{PT_GKLS_Simple, CPT_Normal, 3, 100, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal"},  // по крайней мере лучщий рез
        //{PT_GKLS_Hard, CPT_Normal, 3, 100, 0.002, 4.3, 11, 1.0, "3D_GKLS_Hard_Normal"}  
        //{PT_GKLS_Simple, CPT_Normal, 4, 2, 0.1, 4.7, 10, 1.0, "4D_GKLS_Simple_Normal"}



        // Последовательная версия (tau=1)
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_seq", 1},
        {PT_GKLS_Simple, CPT_Normal, 2, 5, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal_seq", 1},
        //{PT_GKLS_Hard, CPT_Normal, 2, 3, 0.01, 4.5, 10, 1.0, "2D_GKLS_Hard_Normal_seq", 1},
        //{PT_GKLS_Simple, CPT_Normal, 3, 3, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal_seq", 1},
        //{PT_GKLS_Simple, CPT_Normal, 4, 3, 0.01, 4.7, 10, 1.0, "4D_GKLS_Simple_Normal_seq", 1},

        // Параллельная версия (tau=4)
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_par2", 2},
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_par4", 4},
        //{PT_Grishagin, CPT_Normal, 2, 20, 0.01, 4.5, 10, 1.0, "2D_Grishagin_Normal_par8", 8},
        {PT_GKLS_Simple, CPT_Normal, 2, 5, 0.01, 4.5, 10, 1.0, "2D_GKLS_Simple_Normal_par4", 4},
        //{PT_GKLS_Hard, CPT_Normal, 2, 3, 0.01, 4.5, 10, 1.0, "2D_GKLS_Hard_Normal_par4", 4},
        //{PT_GKLS_Simple, CPT_Normal, 3, 3, 0.01, 4.7, 10, 1.0, "3D_GKLS_Simple_Normal_par4", 4},
        //{PT_GKLS_Simple, CPT_Normal, 4, 3, 0.1, 4.7, 10, 1.0, "4D_GKLS_Simple_Normal_par8", 4}

    };

    for (const auto& test : tests) {
        string filename = "results_" + test.name + ".csv";
        runConstrainedTest(test, filename);
    }
    return 0;
}