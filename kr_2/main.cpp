#include <iostream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cfloat>
#include <numeric>
#include <regex>

#define FMT_HEADER_ONLY
#include <fmt/core.h>

#include "fraction_number.hpp"
#include "SLE.hpp"

template<typename T>
using Matrix = std::vector<std::vector<T>>;

std::vector<std::string> reader(std::regex &pattern) {
    std::string input;
    std::cin >> std::ws;
    getline(std::cin, input);
    auto words_begin = std::sregex_iterator(input.begin(), input.end(), pattern);
    auto words_end = std::sregex_iterator();
    std::vector<std::string> res;
    for (std::sregex_iterator match = words_begin; match != words_end; ++match) {
        res.push_back(match->str());
    }
    return res;
}

void printMatrix(Matrix<std::string> matrix) {
    std::vector<size_t> indents(matrix[0].size(), 0);
    for (auto &row: matrix) {
        for (int64_t j = 0; j < row.size(); ++j) {
            indents[j] = std::max(indents[j], row[j].length() + 1);
        }
    }
    fmt::print("\n");
    for (auto row: matrix) {
        for (int64_t i = 0; i < row.size(); ++i) {
            std::string indent(indents[i] - row[i].size(), ' ');
            fmt::print("{}", indent + row[i]);
        }
        fmt::print("\n");
    }
}

void printMatrix(Matrix<size_t> matrix) {
    std::vector<size_t> indents(matrix[0].size(), 0);
    for (auto &row: matrix) {
        for (int64_t j = 0; j < row.size(); ++j) {
            indents[j] = std::max(indents[j], std::to_string(row[j]).length() + 1);
        }
    }
    fmt::print("\n");
    for (auto row: matrix) {
        for (int64_t i = 0; i < row.size(); ++i) {
            std::string indent(indents[i] - std::to_string(row[i]).size(), ' ');
            fmt::print("{}", indent + std::to_string(row[i]));
        }
        fmt::print("\n");
    }
}

class SimplexMethodTask {
    bool needs_normalization_{true};
    std::vector<int64_t> f_coeffs;
    SLE<FractionNumber> g_i_coeffs;

public:
    SimplexMethodTask() {
        size_t x_num, g_i_num;
        std::cin >> x_num >> g_i_num >> needs_normalization_;

        f_coeffs = std::vector<int64_t>(x_num, 0);
        std::regex pattern(R"(-?\d+(\/-?\d+)?)");
        auto data = reader(pattern);
        for (size_t i = 0; i < x_num; ++i) {
            f_coeffs[i] = std::stoi(data[i]);
        }
        Matrix<FractionNumber> equations(g_i_num,
                                            std::vector<FractionNumber>(x_num + 1, FractionNumber()));
        for (auto &equality: equations) {
            data = reader(pattern);
            for (size_t j = 0; j < equality.size(); ++j) {
                equality[j] = FractionNumber(data[j]);
            }
        }

        fmt::print("Inputted functions:\n", f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3]);
        fmt::print("f(x) = {}*x_1 + {}*x_2 + {}*x_3 + {}*x_4\n", f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3]);
        for (auto g_i_coeff: equations) {
            fmt::print("{}*x_1 + {}*x_2 + {}*x_3 + {}*x_4 = {}\n", g_i_coeff[0], g_i_coeff[1],
                       g_i_coeff[2].str(),
                       g_i_coeff[3].str(), g_i_coeff[4].str());
        }

        g_i_coeffs = SLE<FractionNumber>(equations);
    }

    void Solve() {
        std::pair<size_t, size_t> vars;
        if (needs_normalization_) {
            vars = g_i_coeffs.nonNegativeDiagonalization();
        } else {
            for (size_t i = 0; i < g_i_coeffs[0].size() - 1; ++i) {
                if (g_i_coeffs[0][i] != 0 && g_i_coeffs[1][i] == 0) {
                    vars.first = i;
                    break;
                }
            }
            for (size_t i = 0; i < g_i_coeffs[1].size() - 1; ++i) {
                if (g_i_coeffs[1][i] != 0 && g_i_coeffs[0][i] == 0) {
                    vars.second = i;
                    break;
                }
            }
        }

        fmt::print("Chosen base vars: x_{}, x_{}\n", vars.first + 1, vars.second + 1);

        std::vector<FractionNumber> base_1 = g_i_coeffs[0].items_;
        base_1.pop_back();
        std::vector<FractionNumber> base_2 = g_i_coeffs[1].items_;
        base_2.pop_back();
        std::vector<FractionNumber> basis_solve = {g_i_coeffs[0][4], g_i_coeffs[1][4]};
        std::vector<FractionNumber> z_i(f_coeffs.size(), FractionNumber());
        std::vector<FractionNumber> deltas(f_coeffs.size(), FractionNumber());
        std::vector<FractionNumber> a_i(2, FractionNumber());
        while (true) {
            FractionNumber max_delta(0);
            int64_t pos;
            for (int i = 0; i < z_i.size(); ++i) {
                z_i[i] = f_coeffs[vars.first] * base_1[i] + f_coeffs[vars.second] * base_2[i];
                deltas[i] = f_coeffs[i] - z_i[i];
                if (max_delta < deltas[i]) {
                    pos = i;
                    max_delta = deltas[i];
                }
            }

            if (max_delta > 0) {
                a_i[0] = basis_solve[0] / base_1[pos];
                a_i[1] = basis_solve[1] / base_2[pos];
            } else {
                a_i[0] = 0;
                a_i[1] = 0;
            }
            fmt::print("\n\n\n");
            printMatrix({
                                {"",                                    "",                       "",                   std::to_string(
                                        f_coeffs[0]),                                                                                    std::to_string(
                                        f_coeffs[1]),                                                                                                     std::to_string(
                                        f_coeffs[2]),                                                                                                                      std::to_string(
                                        f_coeffs[3]),                                                                                                                                       ""},
                                {"c_i",                                 "BP",                     "BR",                 "x_1",           "x_2",           "x_3",           "x_4",           "a_i"},
                                {std::to_string(f_coeffs[vars.first]),  "x_" + std::to_string(vars.first +
                                                                                              1), basis_solve[0].str(), base_1[0].str(), base_1[1].str(), base_1[2].str(), base_1[3].str(), a_i[0].str()},
                                {std::to_string(f_coeffs[vars.second]), "x_" + std::to_string(vars.second +
                                                                                              1), basis_solve[1].str(), base_2[0].str(), base_2[1].str(), base_2[2].str(), base_2[3].str(), a_i[1].str()},
                                {"",                                    "",                       "z_i",                z_i[0].str(),    z_i[1].str(),    z_i[2].str(),    z_i[3].str(),    ""},
                                {"",                                    "",                       "delta",              deltas[0].str(), deltas[1].str(), deltas[2].str(), deltas[3].str(), ""},
                        });
//            fmt::print("             {}   {}   {}   {}\n", f_coeffs[0], f_coeffs[1], f_coeffs[2], f_coeffs[3]);
//            fmt::print(" c  BP   BR  x_1  x_2  x_3  x_4  a_i\n");
//            fmt::print("{}  x_{}   {}   {}    {}    {}    {}    {}\n", f_coeffs[vars.first], vars.first + 1, basis_solve[0].str(), base_1[0].str(), base_1[1].str(), base_1[2].str(), base_1[3].str(), a_i[0].str());
//            fmt::print("{}  x_{}   {}   {}    {}    {}    {}    {}\n", f_coeffs[vars.second], vars.second + 1, basis_solve[1].str(), base_2[0].str(), base_2[1].str(), base_2[2].str(), base_2[3].str(), a_i[1].str());
//            fmt::print("        z_i   {}   {}   {}   {}\n", z_i[0].str(), z_i[1].str(), z_i[2].str(), z_i[3].str());
//            fmt::print("      delta   {}   {}   {}   {}\n", deltas[0].str(), deltas[1].str(), deltas[2].str(), deltas[3].str());

            if (max_delta <= 0) {
                fmt::print("Max delta is non-positive, stop\n");
                break;
            }

            bool first_line = true;
            std::vector<FractionNumber> &stale = base_1;
            std::vector<FractionNumber> &removable = base_2;

            fmt::print("Adding x_{} as new basis var (this line has max(delta | delta > 0))\n", pos + 1);
            if (a_i[0] <= a_i[1] && a_i[0] > 0 || a_i[1] < 0) {
                fmt::print("Instead of x_{} (this line has min(a_i | a_i > 0))\n", vars.first + 1);
            } else {
                first_line = false;
                std::swap(stale, removable);
                fmt::print("Instead of x_{} (this line has min(a_i | a_i > 0))\n", vars.second + 1);
            }

            FractionNumber support_element = stale[pos];
            FractionNumber mul_supp = removable[pos];
            for (int i = 0; i < removable.size(); ++i) {
                removable[i] = removable[i] - (mul_supp * stale[i] / support_element);
            }
            basis_solve[first_line] = basis_solve[first_line] - (mul_supp * basis_solve[!first_line] / support_element);

            for (FractionNumber &elem: stale) {
                elem /= support_element;
            }
            basis_solve[!first_line] /= support_element;

            if (first_line) {
                vars.first = pos;
            } else {
                std::swap(stale, removable);
                vars.second = pos;
            }
        }

        std::vector<FractionNumber> answer(4, FractionNumber(0));
        answer[vars.first] = basis_solve[0];
        answer[vars.second] = basis_solve[1];

        fmt::print("x = ({})\n", fmt::join(answer, ", "));
    }
};

class TransportationMethodTask {
    size_t cols_{};
    size_t rows_{};
    Matrix<int64_t> c_ij_;
    Matrix<size_t> units_;
    std::vector<size_t> needs_;
    std::vector<size_t> stocks_;

public:
    TransportationMethodTask() {
        std::cin >> rows_ >> cols_;
        c_ij_ = Matrix<int64_t>(rows_, std::vector<int64_t>(cols_, 0));
        units_ = Matrix<size_t>(rows_, std::vector<size_t>(cols_, 0));
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                std::cin >> c_ij_[i][j];
            }
        }
        needs_ = std::vector<size_t>(cols_, 0);
        for (size_t i = 0; i < cols_; ++i) {
            std::cin >> needs_[i];
        }
        stocks_ = std::vector<size_t>(rows_, 0);
        for (size_t i = 0; i < rows_; ++i) {
            std::cin >> stocks_[i];
        }
    }

    void calculateCosts() {
        int64_t cost = 0;
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                cost += units_[i][j] * c_ij_[i][j];
            }
        }
        fmt::print("Transportation costs: {}\n", cost);
    }

    void Solve() {
        std::vector<size_t> current_stocks = stocks_;
        std::vector<size_t> current_needs = needs_;
        for (size_t i = 0; i < rows_; ++i) {
            for (size_t j = 0; j < cols_; ++j) {
                size_t fulfillment = std::min(current_stocks[i], current_needs[j]);
                units_[i][j] = fulfillment;
                current_stocks[i] -= fulfillment;
                current_needs[j] -= fulfillment;
            }
        }

        fmt::print("Initial matrix:\n");

        while (true) {
            printMatrix(units_);
            calculateCosts();

            break;
        }

    }



};

class MulticriteriaOptimizationTask {
public:
    struct Point {
        size_t x;
        size_t y;

        Point() : x(0), y(0) {}

        bool operator>=(Point &other) const {
            return x >= other.x && y >= other.y;
        }

        bool operator<=(Point &other) const {
            return x <= other.x && y <= other.y;
        }

        bool operator!=(Point &other) const {
            return x != other.x || y != other.y;
        }

        bool operator==(Point &other) const {
            return !(*this != other);
        }

        bool Comparable(Point &other) const {
            return *this <= other || *this >= other;
        }

        Point &GetMax(Point &other) {
            if (!Comparable(other)) {
                throw std::runtime_error("Points are not comparable");
            }
            if (*this <= other) {
                return other;
            }
            return *this;
        }
    };

    MulticriteriaOptimizationTask() {
        size_t points_num;
        std::cin >> points_num;
        points = std::vector<Point>(points_num);
        for (size_t i = 0; i < points_num; ++i) {
            std::cin >> points[i].x >> points[i].y;
        }
    }

    void Solve() {
        std::vector<std::pair<size_t, Point>> optimal_points;
        for (size_t i = 0; i < points.size(); ++i) {
            bool is_max_point = true;
            for (size_t j = 0; j < points.size(); ++j) {
                if (i == j) {
                    continue;
                }
                if (!points[i].Comparable(points[j])) {
                    continue;
                }
                if (points[i].GetMax(points[j]) != points[i]) {
                    is_max_point = false;
                    break;
                }
            }
            if (is_max_point) {
                optimal_points.emplace_back(i, points[i]);
            }
        }

        fmt::print("Optimal points:\n");
        for (auto point: optimal_points) {
            fmt::print("#{}: ({}, {})\n", point.first, point.second.x, point.second.y);
        }

        Point utopian_point;
        for (Point &point: points) {
            utopian_point.x = std::max(utopian_point.x, point.x);
            utopian_point.y = std::max(utopian_point.y, point.y);
        }

        fmt::print("Utopian point: ({}, {})\n", utopian_point.x, utopian_point.y);

        long double best_distance = DBL_MAX;
        std::vector<std::pair<size_t, Point>> ideal_points;
        for (size_t i = 0; i < points.size(); ++i) {
            long double distance = sqrtl(
                    powl(utopian_point.x - points[i].x, 2) + powl(utopian_point.y - points[i].y, 2));
            if (distance == best_distance) {
                ideal_points.emplace_back(i, points[i]);
            } else if (distance < best_distance) {
                best_distance = distance;
                ideal_points.clear();
                ideal_points.emplace_back(i, points[i]);
            }
        }

        fmt::print("Ideal points:\n");
        for (auto point: ideal_points) {
            fmt::print("#{}: ({}, {})\n", point.first, point.second.x, point.second.y);
        }
    }

    std::vector<Point> points;
};


int main() {
    SimplexMethodTask task;
    task.Solve();
    return 0;
}

