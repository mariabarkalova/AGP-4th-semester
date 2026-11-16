
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "../evolvent/evolvent.hpp"
#include "GKLS/GKLSConstrainedProblem.hpp"
#include "Grishagin/GrishaginConstrainedProblem.hpp"
#include "PointConstrained.h"
#include "IntervalConstrained.h"

#ifdef _OPENMP
#include <omp.h>
#endif

enum ProblemType {
    PT_Grishagin = 0,
    PT_GKLS_Simple = 1,
    PT_GKLS_Hard = 2 
};

enum ConstraintProblemType {
    CPT_Normal = 0,          // обычные огр
    CPT_InFeasibleDomain = 1,// оптимум внутри недопустимой области
    CPT_OutFeasibleDomain = 2,// снаружи
    CPT_OnFeasibleBorder = 3 // на границе
};

//причина остановки
enum StopReason {
    STOP_NOT_STARTED = 0, // ещё не запускали
    STOP_ACCURACY = 1,  // достигнута точность
    STOP_POINT_LIMIT = 2,   // исчерпан лимит испытаний
    STOP_ITERATION_LIMIT = 3,// исчерпан лимит итер
    STOP_NO_INTERVALS = 4   // не осталось интервалов для дробления
};

class GlobalSearchConstrained
{
private:
    // Для каждого индекса ню (0, 1, 2) храним свой std::map точек, свою оценку мю_ню и своё z*_ню.
    struct MuData {
        std::map<double, double> samples;  // отсортированные по x точки
        double slope;  // текущая оценка mu_nu
        double best_z; // лучшее знач z*_ню
        double worst_z; //худшее

        MuData() {
            slope = 0.0;
            best_z = std::numeric_limits<double>::infinity();
            worst_z = -std::numeric_limits<double>::infinity();
        }

        // Возвращаем mu, но не меньше 1.0 (защита от нулевого наклона)
        double getMu() const { 
            return slope > 0.0 ? slope : 1.0; 
        }
    };


    struct Rank {
        double R;
        double left_x;

        bool operator<(const Rank& other) const {
            if (R != other.R)
                return R > other.R; // по убыванию R
            return left_x < other.left_x;
        }
    };

    //результат одного испытания.
    struct TrialResult {
        PointConstrained point; // (x, z, nu)
        std::vector<double> multidim_point;
        std::vector<size_t> evaluations; // сколько раз вызвали каждую функцию
    };




    double r;
    double epsilon;
    int N;
    int tau; // число стартовых интервалов = размер пачки = число потоков
    int criterion_index; // индекс критерия

    IConstrainedOptProblem* problem;
    ags::Evolvent* evolvent;

    std::vector<MuData> mu_data; // по одному MuData на каждый nu
    std::map<double, IntervalConstrained> intervals; // все интервалы, ключ — left.x
    std::set<Rank> ranking;  // ранжирование интервалов по R

    PointConstrained min_point;   // лучшая допустимая точка (x, z, nu)
    std::vector<double> best_point_multidim; // её N-мерный образ y
    std::vector<double> true_opt_point;      // эталонный оптимум
    double true_opt_value;

    // счётчики и лимиты
    std::vector<size_t> evaluation_counts; // сколько раз вызвана каждая функция
    size_t point_count;   // всего испытаний
    size_t iteration_count; // всего итераций
    size_t max_points;
    size_t max_iterations;
    StopReason stop_reason; // причина остановки
    int actual_threads;// сколько потоков реально было


    // по enum'ам создаёт нужную задачу.
    static IConstrainedOptProblem* createProblem(ProblemType type, ConstraintProblemType constrType,
        int function_number, int dim, double fraction) {
        EConstrainedProblemType cpt = cptNormal;
        if (constrType == CPT_InFeasibleDomain) 
            cpt = cptInFeasibleDomain;
        else if (constrType == CPT_OutFeasibleDomain) 
            cpt = cptOutFeasibleDomain;
        else if (constrType == CPT_OnFeasibleBorder)
            cpt = cptOnFeasibleBorder;

        if (type == PT_Grishagin) {
            return new GrishaginConstrainedProblem(cpt, fraction, 0, function_number);
        }
        return new TGKLSConstrainedProblem(cpt, fraction, 0, function_number, dim, (type == PT_GKLS_Simple) ? Simple : Hard, TD);
    }

    // одно испытание в точке x
    // Функция возвращает значение и тип точки
    TrialResult computeFunctionAndType(double x) const {
        TrialResult result;
        result.multidim_point.resize(N);
        result.evaluations.assign(3, 0);

        // x -> y (N-мерная точка)
        evolvent->GetImage(x, result.multidim_point.data());
        // ограничение 0
        double c0 = problem->ComputeConstraint(0, result.multidim_point);
        result.evaluations[0]++;
        if (c0 > 0) {
            result.point = PointConstrained(x, c0, 0);
            return result;
        }
        // ограничение 1
        double c1 = problem->ComputeConstraint(1, result.multidim_point);
        result.evaluations[1]++;
        if (c1 > 0) {
            result.point = PointConstrained(x, c1, 1);
            return result;
        }
        // все ограничения выполнены — считаем критерий
        result.evaluations[2]++;
        double z = problem->ComputeFunction(result.multidim_point);
        result.point = PointConstrained(x, z, 2);

        return result;
    }


    bool updateMu(double p_x, double p_z, double sample_x, double sample_z, double max_diff, MuData& data) const {
        double lo = std::min(p_x, sample_x);
        double hi = std::max(p_x, sample_x);
        double delta = IntervalConstrained(PointConstrained(lo, 0.0), PointConstrained(hi, 0.0)).getDelta(N);

        // Ранний выход: оставшиеся точки ещё дальше, а max_diff уже не может дать наклон больше текущего slope
        if (data.slope > 0.0 && max_diff / delta <= data.slope) 
            return false;

        double candidate = std::abs(p_z - sample_z) / delta;
        if (candidate > data.slope) data.slope = candidate;
        return true;
    }


    bool insertPointSorted(const TrialResult& result) {
        const PointConstrained& p = result.point;
        MuData& data = mu_data[p.nu];
        double old_mu = data.getMu();
        double old_best = data.best_z;
        double max_diff = std::max(std::abs(p.z - data.best_z), std::abs(p.z - data.worst_z));

        // Позиция для вставки в отсортированный map
        std::map<double, double>::iterator right = data.samples.lower_bound(p.x); //lower_bound(p.x) — первый элемент с ключом ? p.x

        // Обход вправо: пока updateMu возвращает true
        for (std::map<double, double>::iterator it = right; it != data.samples.end(); ++it) {
            if (!updateMu(p.x, p.z, it->first, it->second, max_diff, data)) break;
        }
        // Обход влево
        std::map<double, double>::iterator left = right;
        while (left != data.samples.begin()) {
            --left;
            if (!updateMu(p.x, p.z, left->first, left->second, max_diff, data)) break;
        }

        // Вставляем саму точку
        data.samples.insert(right, std::make_pair(p.x, p.z));
        data.best_z = std::min(data.best_z, p.z);
        data.worst_z = std::max(data.worst_z, p.z);

        // Обнов-м счётчики вызовов
        for (size_t nu = 0; nu < evaluation_counts.size(); ++nu)
            evaluation_counts[nu] += result.evaluations[nu];
        ++point_count;

        // Обнов-м глоб лучший рез (только допустимая т.)
        if (p.nu == criterion_index && p.z < min_point.z) {
            min_point = p;
            best_point_multidim = result.multidim_point;
        }

        // Сообщаем поменялись ли mu или best_z
        return old_mu != data.getMu() || old_best != data.best_z;
    }


    double kharakteristika_R(const IntervalConstrained& interval) const {
        // Невалидный интервал не участвует в ранжировании
        if (!interval.isValid()) return -std::numeric_limits<double>::infinity();

        double delta = interval.getDelta(N);
        const PointConstrained& a = interval.left;
        const PointConstrained& b = interval.right;
        int nu = std::max(a.nu, b.nu);       // «тип» интервала
        double M = r * mu_data[nu].getMu();  // r * mu_nu
        double z_star = mu_data[nu].best_z;  // лучшее z для этого nu

        if (a.nu == b.nu) {
            // Оба края одного типа 
            double diff = (b.z - a.z) / M;
            return delta + diff * diff / delta - 2.0 * ((a.z - z_star) / M + (b.z - z_star) / M);
        }
        // Края разных типов
        double z = (a.nu > b.nu) ? a.z : b.z;
        return 2.0 * delta - 4.0 * (z - z_star) / M;
    }


    void ocenkaConstL() {
        ranking.clear();
        std::map<double, IntervalConstrained>::iterator it = intervals.begin();
        for (; it != intervals.end(); ++it) {
            IntervalConstrained& interval = it->second;
            interval.R = kharakteristika_R(interval);
            if (std::isfinite(interval.R)) {
                Rank rk;
                rk.R = interval.R;
                rk.left_x = interval.left.x;
                ranking.insert(rk);
            }
        }
    }


    void initializeProblem(int tightness) {
        N = problem->GetDimension();
        criterion_index = problem->GetConstraintsNumber();

        std::vector<double> lb, ub;
        problem->GetBounds(lb, ub);
        evolvent = new ags::Evolvent(N, tightness, lb.data(), ub.data(), ags::Simple);

        mu_data.resize(criterion_index + 1);
        evaluation_counts.assign(criterion_index + 1, 0);
        min_point = PointConstrained(0.0, std::numeric_limits<double>::infinity(), criterion_index);

        true_opt_point = problem->GetOptimumPoint();
        true_opt_value = problem->GetOptimumValue();

        // Стартовая сетка: tau+1 равномерных точек, включая 0 и 1
        std::vector<TrialResult> initial(tau + 1);
#pragma omp parallel for if(tau > 1) num_threads(tau) schedule(static)
        for (int s = 0; s <= tau; ++s) {
            initial[s] = computeFunctionAndType(static_cast<double>(s) / tau);
        }

        // Вставляем все стартовые точки в историю
        for (size_t i = 0; i < initial.size(); ++i) {
            insertPointSorted(initial[i]);
        }

        // Строим tau начальных интервалов между соседними точками
        for (int s = 0; s < tau; ++s) {
            intervals.insert(std::make_pair(initial[s].point.x, IntervalConstrained(initial[s].point, initial[s + 1].point)));
        }

        // Первичный пересчёт R и построение ranking
        ocenkaConstL();
    }

public:

    GlobalSearchConstrained(ProblemType type, ConstraintProblemType constrType, int function_number, int dim = 2, double epsilon_ = 0.01,
        double r_ = 4.5, int tightness = 12, double fraction = 0.3, int tau_ = 4)
        : GlobalSearchConstrained(createProblem(type, constrType, function_number, dim, fraction), epsilon_, r_, tightness, tau_) {}


    GlobalSearchConstrained(IConstrainedOptProblem* problem_, double epsilon_ = 0.01, double r_ = 4.5, int tightness = 10, int tau_ = 1)
    {
        r = r_;
        epsilon = epsilon_;
        N = 0;
        tau = tau_;
        criterion_index = 0;
        problem = problem_;
        evolvent = 0;
        point_count = 0;
        iteration_count = 0;
        max_points = 200000;
        max_iterations = 200000;
        stop_reason = STOP_NOT_STARTED;
        actual_threads = 1;
        true_opt_value = std::numeric_limits<double>::quiet_NaN();

        initializeProblem(tightness);
    }

    ~GlobalSearchConstrained() {
        delete problem;
        delete evolvent;
    }


    void algorithm() {
        // Защита от повторного запуска
        if (stop_reason != STOP_NOT_STARTED) return;

        std::vector<IntervalConstrained> selected; // выбранные интервалы текущей пачки
        int batch_count = 0; // размер пачки
        double min_delta = 0.0;  // мин. ширина среди selected
        std::vector<double> coordinates; // new x
        std::vector<TrialResult> results; // рез испытаний
        bool running = true;  // флаг продолжаем
 
#pragma omp parallel if(tau > 1) num_threads(tau) shared(running, batch_count, min_delta, selected, coordinates, results)
        {
#pragma omp single
            {
#ifdef _OPENMP
                actual_threads = omp_get_num_threads();
#endif
            }

            for (;;) {
                // выбор пачки
#pragma omp single
                {
                    // Проверка лимитов (в начале итерации)
                    if (point_count >= max_points) {
                        stop_reason = STOP_POINT_LIMIT;
                        running = false;
                    }
                    else if (iteration_count >= max_iterations) {
                        stop_reason = STOP_ITERATION_LIMIT;
                        running = false;
                    }
                    else if (ranking.empty()) {
                        stop_reason = STOP_NO_INTERVALS;
                        running = false;
                    }

                    if (running) {
                        size_t room = max_points - point_count; //сколько еще слотов осталось до лимита
                        size_t want = static_cast<size_t>(tau); //сколько хотим
                        if (want > ranking.size()) 
                            want = ranking.size();
                        if (want > room) 
                            want = room;
                        batch_count = static_cast<int>(want); //итоговый размер пачки

                        selected.clear();
                        coordinates.clear();
                        results.resize(batch_count);
                        min_delta = std::numeric_limits<double>::infinity();

                        // Берём batch_count лучших интервалов по R
                        for (int s = 0; s < batch_count; ++s) {
                            Rank candidate = *ranking.begin(); // максимальный R
                            ranking.erase(ranking.begin());    // убираем из ranking
                            selected.push_back(intervals.at(candidate.left_x));
                            const IntervalConstrained& interval = selected.back();

                            // Обновляем min_delta
                            double d = interval.getDelta(N);
                            if (d < min_delta) min_delta = d;

                            // Новая точка внутри интервала
                            int nu = std::max(interval.left.nu, interval.right.nu);
                            coordinates.push_back(interval.getNewPoint(mu_data[nu].getMu(), r, N));
                        }
                    }
                }
                if (!running) break;

                //параллельные испытания
#pragma omp for schedule(static)
                for (int s = 0; s < batch_count; ++s) {
                    results[s] = computeFunctionAndType(coordinates[s]);
                }

#pragma omp single
                {
                    bool changed = false;

                    //для каждого рез-та: вставить в историю и разрезать интервал
                    for (int s = 0; s < batch_count; ++s) {
                        const IntervalConstrained& interval = selected[s]; //выбранный интервал
                        const PointConstrained& point = results[s].point; //новая т
                        //Доб новую точку. 
                        // Внутри insertPointSorted вызывается updateMu.
                        // Если mu изменилась, возвращается true.
                        if (insertPointSorted(results[s])) 
                            changed = true;

                        // Удаляем старый интервал, вставляем два новых
                        intervals.erase(interval.left.x);
                        intervals.insert(std::make_pair(interval.left.x, IntervalConstrained(interval.left, point)));
                        intervals.insert(std::make_pair(point.x, IntervalConstrained(point, interval.right)));
                    }

                    // Пересчёт R:
                    // если mu/best_z изменились — пересчитываем все R
                    // иначе только R двух новых интервалов на каждую пачку
                    if (changed) {
                        ocenkaConstL();
                    }
                    else {
                        for (int s = 0; s < batch_count; ++s) {
                            // Левая половина
                            IntervalConstrained& i1 = intervals.at(selected[s].left.x);
                            i1.R = kharakteristika_R(i1);
                            if (std::isfinite(i1.R)) {
                                Rank rk; 
                                rk.R = i1.R;
                                rk.left_x = i1.left.x;
                                ranking.insert(rk);
                            }
                            // Правая половина
                            IntervalConstrained& i2 = intervals.at(results[s].point.x);
                            i2.R = kharakteristika_R(i2);
                            if (std::isfinite(i2.R)) {
                                Rank rk;
                                rk.R = i2.R;
                                rk.left_x = i2.left.x;
                                ranking.insert(rk);
                            }
                        }
                    }

                    ++iteration_count;

                    if (hasFeasiblePoint() && min_delta <= epsilon) {
                        stop_reason = STOP_ACCURACY;
                        running = false;
                    }
                }
            }
        }
    }



    // Есть ли хоть одна допустимая точка
    bool hasFeasiblePoint() const {
        return std::isfinite(min_point.z);
    }

    const char* getStopReasonName() const {
        switch (stop_reason) {
        case STOP_ACCURACY:        
            return "accuracy";
        case STOP_POINT_LIMIT:   
            return "point_limit";
        case STOP_ITERATION_LIMIT:
            return "iteration_limit";
        case STOP_NO_INTERVALS:   
            return "no_intervals";
        default:                  
            return "not_started";
        }
    }

    int getActualThreads() const { return actual_threads; }
    size_t getIterationCount() const { return iteration_count; }
    size_t getPointCount() const { return point_count; }
    const std::vector<size_t>& getEvaluationCounts() const { return evaluation_counts; }

    // Результаты поиска
    void getResults(std::vector<double>& x_min, double& f_min, size_t& count,
        std::vector<double>& true_x, double& true_f) const {
        x_min = best_point_multidim;
        f_min = min_point.z;
        count = point_count; // всего испытаний
        true_x = true_opt_point;
        true_f = true_opt_value;
    }

    // Значения ограничений в лучшей точке
    void getConstraintsInfo(std::vector<double>& constraints) const {
        constraints.clear();
        if (hasFeasiblePoint()) {
            for (int i = 0; i < criterion_index; ++i)
                constraints.push_back(problem->ComputeConstraint(i, best_point_multidim));
        }
    }
};












