#pragma once

#include <vector>

#define FMT_HEADER_ONLY

#include <fmt/core.h>

template<class T>
class SLE {
    struct Row {
        std::vector<T> items_;

        Row operator+(const Row &other) const {
            Row row = *this;
            for (size_t i = 0; i < row.items_.size(); ++i) {
                row.items_[i] += other.items_[i];
            }
            return row;
        }

        Row operator-(const Row &other) const {
            Row row = *this;
            for (size_t i = 0; i < row.items_.size(); ++i) {
                row.items_[i] -= other.items_[i];
            }
            return row;
        }

        void operator-=(const Row &other) {
            for (size_t i = 0; i < items_.size(); ++i) {
                items_[i] -= other.items_[i];
            }
        }

        Row operator*(const Row &other) const {
            Row row = *this;
            for (size_t i = 0; i < row.items_.size(); ++i) {
                row.items_[i] *= other.items_[i];
            }
            return row;
        }

        Row operator*(const T &other) const {
            Row row = *this;
            for (T &item: row.items_) {
                item *= other;
            }
            return row;
        }

        void operator*=(const Row &other) {
            for (size_t i = 0; i < items_.size(); ++i) {
                items_[i] *= other.items_[i];
            }
        }

        Row operator/(const Row &other) {
            Row row = *this;
            for (size_t i = 0; i < row.items_.size(); ++i) {
                row.items_[i] /= other.items_[i];
            }
            return row;
        }

        Row operator/(const T &other) const {
            Row row = *this;
            for (T &item: row.items_) {
                item /= other;
            }
            return row;
        }

        void operator/=(const Row &other) {
            for (size_t i = 0; i < items_.size(); ++i) {
                items_[i] /= other.items_[i];
            }
        }

        T &operator[](size_t pos) {
            return items_[pos];
        }

        size_t size() {
            return items_.size();
        }

        void normalize(size_t pos) {
            T divisor = items_[pos];
            for (T &item: items_) {
                item /= divisor;
            }
        }
    };

    std::vector<Row> system_;

public:
    SLE() = default;

    explicit SLE(std::vector<std::vector<T>> &matrix) {
        for (std::vector<T> &row: matrix) {
            system_.push_back(Row(row));
        }
    }

    Row &operator[](size_t pos) {
        return system_[pos];
    }

    bool checkSignCondition(int i, int j) {
        return (system_[0][i] < 0 == system_[0][4] < 0) && (system_[1][j] < 0 == system_[1][4] < 0);
    }

    std::pair<size_t, size_t> nonNegativeDiagonalization() {
        for (int i = 0; i < system_[0].size() - 1; ++i) {
            for (int j = 0; j < system_[0].size() - 1; ++j) {
                if (i == j) {
                    continue;
                }
                fmt::print("From system:\n");
                fmt::print("[{} {} {} {} | {}]\n", system_[0][0], system_[0][1], system_[0][2], system_[0][3],
                           system_[0][4]);
                fmt::print("[{} {} {} {} | {}]\n", system_[1][0], system_[1][1], system_[1][2], system_[1][3],
                           system_[1][4]);

                if (system_[0][i] != 0) {
                    system_[1] -= system_[0] * (system_[1][i] / system_[0][i]);
                }
                if (system_[1][j] != 0) {
                    system_[0] -= system_[1] * (system_[0][j] / system_[1][j]);
                }
                fmt::print("Transformed to:\n");
                fmt::print("[{} {} {} {} | {}]\n", system_[0][0], system_[0][1], system_[0][2], system_[0][3],
                           system_[0][4]);
                fmt::print("[{} {} {} {} | {}]\n\n", system_[1][0], system_[1][1], system_[1][2], system_[1][3],
                           system_[1][4]);

                if (system_[0][i] != 0 && system_[1][i] == 0 && system_[1][j] != 0 && system_[0][j] == 0) {
                    if (checkSignCondition(i, j)) {
                        system_[0].normalize(i);
                        system_[1].normalize(j);
                        fmt::print("Found solution, normalized form:\n");
                        fmt::print("[{} {} {} {} | {}]\n", system_[0][0], system_[0][1], system_[0][2], system_[0][3],
                                   system_[0][4]);
                        fmt::print("[{} {} {} {} | {}]\n\n", system_[1][0], system_[1][1], system_[1][2], system_[1][3],
                                   system_[1][4]);
                        return {i, j};
                    }
                }
            }
        }
        throw std::runtime_error("Couldn't find the solution for initial basis variables.");
    }

};