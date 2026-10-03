
#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include "../evolvent/evolvent.hpp"
#include "GKLS/GKLSProblem.hpp"
#include "GKLS/GKLSConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "PointConstrained.h"
#include "IntervalConstrained.h"
#include <omp.h>

using namespace std;

enum ProblemType {
    PT_Grishagin = 0,
    PT_GKLS_Simple = 1,
    PT_GKLS_Hard = 2
};

enum ConstraintProblemType {
    CPT_Normal = 0,
    CPT_InFeasibleDomain = 1, //вот это
    CPT_OutFeasibleDomain = 2,
    CPT_OnFeasibleBorder = 3
};

class GlobalSearchConstrained
{
private:
    double r;
    double epsilon;
    int N;
    int m;
    int tau;   // количество испытаний за итерацию (для паралл)

    vector<PointConstrained> points;
    vector<IntervalConstrained> intervals;

    IConstrainedOptProblem* problem;
    ags::Evolvent* evolvent;

    PointConstrained min_point;
    vector<double> best_point_multidim;
    //size_t point_count;

    vector<double> true_opt_point;
    double true_opt_value;

    map<int, double> mu;
    map<int, double> best_z;

    size_t last_points_count_for_mu;
    int iteration_counter;
    size_t point_count;
    size_t iteration_count;
    int mu_update_frequency = 20;

    static constexpr int CONSTR0_IDX = 0;
    static constexpr int CONSTR1_IDX = 1;
    static constexpr int CRITERION_IDX = 2;





    int constr0_count;   // сколько раз считали первое ограничение
    int constr1_count;   //второе ограничение
    int function_count;  //целевую функцию

    vector<vector<double>> best_history;
    vector<double> best_values;
    vector<PointConstrained> all_points_history; //  все точек для Excel

public:

    GlobalSearchConstrained(ProblemType type, ConstraintProblemType constrType, int function_number, int dim = 2, double epsilon_ = 0.01,
        double r_ = 4.5, int tightness_ = 12, double fraction = 0.3, int tau_ = 4)
        : epsilon(epsilon_), r(r_), m(tightness_), tau(tau_), last_points_count_for_mu(0), iteration_counter(0), constr0_count(0), 
        constr1_count(0), function_count(0), iteration_count(0)
    {
        EConstrainedProblemType cpt;
        switch (constrType) {
        case CPT_Normal:
            cpt = cptNormal;
            break;
        case CPT_InFeasibleDomain:
            cpt = cptInFeasibleDomain;
            break;
        case CPT_OutFeasibleDomain:
            cpt = cptOutFeasibleDomain;
            break;
        case CPT_OnFeasibleBorder:
            cpt = cptOnFeasibleBorder;
            break;
        default:
            cpt = cptNormal;
        }

        if (type == PT_Grishagin) {
            problem = new GrishaginConstrainedProblem(cpt, fraction, 0, function_number);
        }
        else if (type == PT_GKLS_Simple) {
            problem = new TGKLSConstrainedProblem(cpt, fraction, 0, function_number, dim, Simple, TD);
        }
        else {
            problem = new TGKLSConstrainedProblem(cpt, fraction, 0, function_number, dim, Hard, TD);
        }

        initializeProblem();
    }

    ~GlobalSearchConstrained()
    {
        delete problem;
        delete evolvent;
    }

private:
    // Функция возвращает значение и тип точки за один вызов
    pair<double, int> computeFunctionAndType(const vector<double>& point)
    {
        double c0 = problem->ComputeConstraint(0, point);
#pragma omp atomic
        constr0_count++;
        if (c0 > 0) return { c0, CONSTR0_IDX };

        double c1 = problem->ComputeConstraint(1, point);
#pragma omp atomic
        constr1_count++;
        if (c1 > 0) return { c1, CONSTR1_IDX };

#pragma omp atomic   //(не знаю надо ли)???????????7
        function_count++;
        return { problem->ComputeFunction(point), CRITERION_IDX };
    }

    void initializeProblem() {
        N = problem->GetDimension();
        vector<double> lb, ub;
        problem->GetBounds(lb, ub);

        evolvent = new ags::Evolvent(N, m, lb.data(), ub.data(), ags::Simple);
        best_point_multidim.resize(N);
        true_opt_point = problem->GetOptimumPoint();
        true_opt_value = problem->GetOptimumValue();
        min_point = PointConstrained(0, 1e100, CRITERION_IDX);

        vector<double> point(N);

        // Точка t = 0
        evolvent->GetImage(0.0, point.data());
        auto [z0, nu0] = computeFunctionAndType(point);
        points.push_back(PointConstrained(0.0, z0, nu0));
        if (nu0 == CRITERION_IDX && z0 < min_point.z) {
            min_point = PointConstrained(0.0, z0, nu0);
            best_point_multidim = point;
        }

        // Точка t = 1
        evolvent->GetImage(1.0, point.data());
        auto [z1, nu1] = computeFunctionAndType(point);
        points.push_back(PointConstrained(1.0, z1, nu1));
        if (nu1 == CRITERION_IDX && z1 < min_point.z) {
            min_point = PointConstrained(1.0, z1, nu1);
            best_point_multidim = point;
        }

        sortPoints();
        updateIntervals();
        //лучшие значения для каждого из огр (ограничение0, ограничение1, критерий)
        for (size_t i = 0; i < points.size(); ++i)
        {
            const PointConstrained& p = points[i];
            if (best_z.find(p.nu) == best_z.end() || p.z < best_z[p.nu]) {
                best_z[p.nu] = p.z;
            }
        }

        ocenkaConstL();
        last_points_count_for_mu = points.size();
        point_count = points.size();
    }


    void sortPoints() {
        sort(points.begin(), points.end());
    }

    void updateIntervals() {
        intervals.clear();
        for (size_t i = 0; i < points.size() - 1; i++) {
            IntervalConstrained interval(points[i], points[i + 1], i);
            intervals.push_back(interval);
        }
    }
 
    void ocenkaConstL() {
        vector<int> indices = { CONSTR0_IDX, CONSTR1_IDX, CRITERION_IDX };
        for (int nu : indices) {
            double max_m = 0.0;
            int count = 0;
            for (size_t i = 0; i < intervals.size(); i++)
            {
                const IntervalConstrained& interval = intervals[i];
                if (interval.left.nu == nu && interval.right.nu == nu)
                {
                    /*double delta = interval.getDelta(N);
                    if (delta > 1e-10)
                    {
                        double delta_z = fabs(interval.right.z - interval.left.z);
                        double m_val = delta_z / delta;
                        if (m_val > max_m) {
                            max_m = m_val;
                        }
                        count++;
                    }*/
                    double delta = interval.length;
                    if (delta > 1e-10) {
                        double delta_z = fabs(interval.right.z - interval.left.z);
                        double m_val = delta_z / delta;
                        if (m_val > max_m) {
                            max_m = m_val;
                        }
                        count++;
                    }
                }
            }

            if (count == 0 || max_m < 1e-10) {
                mu[nu] = 1.0;
            }
            else {
                mu[nu] = max_m;
            }
        }
    }
    /*
    vector<double> kharakteristika_R() {
        iteration_counter++;
        if (points.size() != last_points_count_for_mu || iteration_counter % mu_update_frequency == 0) {
            ocenkaConstL();
            last_points_count_for_mu = points.size();
            iteration_counter = 0;
        }

        vector<double> R_values;
        int k = intervals.size(); //кол-во интервалов
        for (int i = 0; i < k; i++) {
            const auto& interval = intervals[i];
            double delta = interval.getDelta(N);
            double z_left = interval.left.z;
            double z_right = interval.right.z;
            int nu_left = interval.left.nu;
            int nu_right = interval.right.nu;
            double R = 0.0;

            if (nu_left == nu_right) {
                int nu = nu_left;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                double diff = z_right - z_left;
                R = delta + (diff * diff) / (r * r * mu_nu * mu_nu * delta) - 2.0 * (z_right + z_left - 2.0 * z_star) / (r * mu_nu);
            }
            else if (nu_left < nu_right) {
                int nu = nu_right;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                R = 2.0 * delta - 2.0 * (z_right - z_star) / (r * mu_nu);
            }
            else {
                int nu = nu_left;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                R = 2.0 * delta - 2.0 * (z_left - z_star) / (r * mu_nu);
            }

            intervals[i].R = R;
            R_values.push_back(R);

            if (!interval.isValid()) {
                intervals[i].R = -1e100;  // отрицательное, а не положительное
                R_values.push_back(-1e100);
                continue;
            }
        }
        return R_values;
    }*/
    vector<double> kharakteristika_R() {
        iteration_counter++;
        if (points.size() != last_points_count_for_mu || iteration_counter % mu_update_frequency == 0) {
            ocenkaConstL();
            last_points_count_for_mu = points.size();
            iteration_counter = 0;
        }

        vector<double> R_values;
        int k = intervals.size();

        for (int i = 0; i < k; i++) {
            const auto& interval = intervals[i];

            // ВАЖНО: используем length, а не getDelta(N)
            double delta = interval.length;
            if (delta < 1e-12) {
                intervals[i].R = -1e100;
                R_values.push_back(-1e100);
                continue;
            }

            int nu_left = interval.left.nu;
            int nu_right = interval.right.nu;
            double R = 0.0;

            // --- ГРАНИЧНЫЕ ИНТЕРВАЛЫ ---
            if (i == 0) {
                int nu = nu_right;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                R = 2.0 * delta - 4.0 * (interval.right.z - z_star) / (r * mu_nu);
            }
            else if (i == k - 1) {
                int nu = nu_left;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                R = 2.0 * delta - 4.0 * (interval.left.z - z_star) / (r * mu_nu);
            }
            // --- ВНУТРЕННИЕ ИНТЕРВАЛЫ ---
            else if (nu_left == nu_right) {
                int nu = nu_left;
                double mu_nu = mu[nu];
                double z_star = best_z[nu];
                double diff = interval.right.z - interval.left.z;
                R = delta
                    + (diff * diff) / (r * r * mu_nu * mu_nu * delta)
                    - 2.0 * (interval.right.z + interval.left.z - 2.0 * z_star) / (r * mu_nu);
            }
            else if (nu_left < nu_right) {
                int nu = nu_right;
                R = 2.0 * delta - 4.0 * (interval.right.z - best_z[nu]) / (r * mu[nu]);
            }
            else {
                int nu = nu_left;
                R = 2.0 * delta - 4.0 * (interval.left.z - best_z[nu]) / (r * mu[nu]);
            }

            intervals[i].R = R;
            R_values.push_back(R);
        }

        return R_values;
    }




/////////////////////////////////////////////////////////////////////////////////////////////////
    // Находит ? интервалов с наибольшими характеристиками R
    /*vector<size_t> findTopTIntervals(int t) {
        int k = intervals.size();
        if (k == 0) return {};

        vector<size_t> indices(k);
        std::iota(indices.begin(), indices.end(), 0);

        int count = std::min(t, k);
        std::partial_sort(indices.begin(), indices.begin() + count, indices.end(),
            [this](size_t a, size_t b) {
                return intervals[a].R > intervals[b].R;
            });

        indices.resize(count);
        return indices;
    }*/
    // это вроде проще записано
    vector<size_t> findTopTIntervals(int t) {
        vector<size_t> result;
        vector<bool> used(intervals.size(), false);

        for (int s = 0; s < t; s++) {
            size_t best = 0;
            double best_R = -1e100;
            for (size_t i = 0; i < intervals.size(); i++) {
                if (!used[i] && intervals[i].R > best_R) {
                    best_R = intervals[i].R;
                    best = i;
                }
            }
            used[best] = true;
            result.push_back(best);
        }
        return result;
    }

    // Результат одного испытания
    struct TrialResult {
        double new_x;
        double new_z;
        int new_nu;
        vector<double> multidim_point;
    };

    // Провести одно испытание (вызывается из параллельного потока)
    TrialResult performTrial(size_t interval_index) {
        const IntervalConstrained& interval = intervals[interval_index];

        int nu;
        if (interval.left.nu == interval.right.nu) {
            nu = interval.left.nu;
        }
        else {
            nu = max(interval.left.nu, interval.right.nu);
        }

        double new_x = interval.getNewPoint(mu[nu], r, N, intervals.size());

        vector<double> multidim_point(N);
        evolvent->GetImage(new_x, multidim_point.data());

        auto [new_z, new_nu] = computeFunctionAndType(multidim_point);

        return { new_x, new_z, new_nu, multidim_point };
    }
    ///////////////////////////////////////////////////////////////////////////////////////






    size_t findMaxRInterval() {
        if (intervals.empty()) return intervals.size(); // сигнал вызывающему
        size_t max_index = 0;
        double max_R = intervals[0].R;
        for (size_t i = 1; i < intervals.size(); i++) {
            if (intervals[i].R > max_R) {
                max_R = intervals[i].R;
                max_index = i;
            }
        }
        return max_index;
    }
   
    // Вставка точки с сохранением сортировки
    void insertPointSorted(const PointConstrained& new_point, const vector<double>& multidim_point) {
        auto pos = lower_bound(points.begin(), points.end(), new_point);
        points.insert(pos, new_point);
        updateIntervals();

        if (best_z.find(new_point.nu) == best_z.end() || new_point.z < best_z[new_point.nu]) {
            best_z[new_point.nu] = new_point.z;
        }
        if (new_point.nu == CRITERION_IDX && new_point.z < min_point.z) {
            min_point = new_point;
            best_point_multidim = multidim_point;
        }
    }

    // Функция для вычисления расстояния до истинного минимума
    double distanceToTrueOptimum(const vector<double>& point) const {
        double dist = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = point[i] - true_opt_point[i];
            dist += diff * diff;
        }
        return sqrt(dist);
    }

public:
    /*
    void algorithm() {
        int iteration = 0;
        const int MAX_ITERATIONS = 20000;
        const int MAX_POINTS = 20000;

        while (iteration < MAX_ITERATIONS && points.size() < MAX_POINTS) {
            iteration++;

            kharakteristika_R();
            size_t t = findMaxRInterval();
            if (t >= intervals.size()) break;

            IntervalConstrained& max_interval = intervals[t];

            int nu;
            if (max_interval.left.nu == max_interval.right.nu) {
                nu = max_interval.left.nu;
            }
            else {
                nu = max(max_interval.left.nu, max_interval.right.nu);
            }

            double new_x = max_interval.getNewPoint(mu[nu], r, N, intervals.size());
            double delta = max_interval.getDelta(N);

            vector<double> multidim_point(N);
            evolvent->GetImage(new_x, multidim_point.data());

            //double new_z = evaluateFunction(multidim_point);
            //int new_nu = getFunctionIndex(multidim_point);
            auto [new_z, new_nu] = computeFunctionAndType(multidim_point);

            PointConstrained new_point(new_x, new_z, new_nu);

            double old_best = min_point.z;
            insertPointSorted(new_point, multidim_point);

            //остановка
            // попадание в окрестность
            double dist_to_true = distanceToTrueOptimum(best_point_multidim);
            if (dist_to_true <= epsilon) {
                break;
            }
        }

        point_count = points.size();
    }*/



    void algorithm() {
        int iteration = 0;
        iteration_count = 0;
        const int MAX_ITERATIONS = 200000;
        const int MAX_POINTS = 200000;

#ifdef _OPENMP
        cout << "OpenMP: " << omp_get_max_threads() << " потоков, tau=" << tau << endl;
#else
        cout << "OpenMP НЕ включён" << endl;
#endif

        while (iteration < MAX_ITERATIONS && points.size() < MAX_POINTS) {
            iteration++;
            iteration_count++;
            // Шаги 1-3: характеристики (последовательно)
            kharakteristika_R();

            // Шаг 4: ? лучших интервалов
            int t_count = std::min(tau, (int)intervals.size());
            if (t_count == 0) break;
            vector<size_t> top_indices = findTopTIntervals(t_count);

            // Шаг 5: параллельные испытания
            vector<TrialResult> results(t_count);

#pragma omp parallel for
            for (int s = 0; s < t_count; s++) {
                results[s] = performTrial(top_indices[s]);
            }

            // Последовательная вставка
            for (const auto& res : results) {
                PointConstrained new_point(res.new_x, res.new_z, res.new_nu);
                insertPointSorted(new_point, res.multidim_point);
            }

            // Шаг 6: остановка
            double dist_to_true = distanceToTrueOptimum(best_point_multidim);
            if (dist_to_true <= epsilon) {
                break;
            }
        }

        point_count = points.size();
    }

    void getResults(vector<double>& x_min, double& f_min, size_t& count,
        vector<double>& true_x, double& true_f) {
        x_min = best_point_multidim;
        f_min = min_point.z;
        count = point_count;
        true_x = true_opt_point;
        true_f = true_opt_value;
    }
    size_t getIterationCount() const { return iteration_count; }
    size_t getPointCount() const { return point_count; }

    void getConstraintsInfo(vector<double>& constraints) {
        constraints.clear();
        for (int i = 0; i < problem->GetConstraintsNumber(); i++) {
            constraints.push_back(problem->ComputeConstraint(i, best_point_multidim));
        }
    }
};