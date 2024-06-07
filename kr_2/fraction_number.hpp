#pragma once

#include <cstdint>
#include <string>
#include <numeric>

#define FMT_HEADER_ONLY

#include <fmt/core.h>
#include <fmt/format.h>

class FractionNumber {
private:
    int64_t numerator;
    int64_t denominator;

public:
    FractionNumber() : numerator(0), denominator(1) {}

    explicit FractionNumber(int64_t seed) : numerator(seed), denominator(1) {}

    explicit FractionNumber(std::string &seed) {
        numerator = std::stoi(seed);
        if (seed.contains("/")) {
            denominator = std::stoi(std::string(seed.begin() + seed.find('/') + 1, seed.end()));
        } else {
            denominator = 1;
        }
    }

    FractionNumber &operator=(int64_t seed) {
        numerator = seed;
        denominator = 1;
        return *this;
    }

    bool operator<(const FractionNumber &other) const {
        int64_t lcm = std::lcm(denominator, other.denominator);
        return numerator * (lcm / denominator) < other.numerator * (lcm / other.denominator);
    }

    bool operator>(const FractionNumber &other) const {
        return other < *this;
    }

    bool operator<=(const FractionNumber &other) const {
        return !(other < *this);
    }

    bool operator>=(const FractionNumber &other) const {
        return !(*this < other);
    }

    bool operator==(const FractionNumber &other) const {
        return numerator == other.numerator && denominator == other.denominator;
    }

    bool operator!=(const FractionNumber &other) const {
        return !(*this == other);
    }

    bool operator<(const int64_t other) const {
        return numerator < other * denominator;
    }

    bool operator>(const int64_t other) const {
        return numerator > other * denominator;
    }

    bool operator<=(const int64_t other) const {
        return !(*this > other);
    }

    bool operator>=(const int64_t other) const {
        return !(*this < other);
    }

    bool operator==(const int64_t other) const {
        return numerator == other && denominator == 1;
    }

    bool operator!=(const int64_t other) const {
        return !(*this == other);
    }

    friend FractionNumber operator+(const FractionNumber &a, const FractionNumber &b);

    friend FractionNumber operator+(const FractionNumber &a, int64_t other);

    friend FractionNumber operator+(int64_t other, const FractionNumber &a);

    friend FractionNumber operator-(const FractionNumber &a, const FractionNumber &b);

    friend FractionNumber operator-(const FractionNumber &a, int64_t other);

    friend FractionNumber operator-(int64_t other, const FractionNumber &a);

    void operator-=(const FractionNumber &other) {
        int64_t lcm = std::lcm(denominator, other.denominator);
        numerator = numerator * (lcm / denominator) - other.numerator * (lcm / other.denominator);
        denominator = lcm;
        reduction();
    }

    friend FractionNumber operator*(const FractionNumber &a, const FractionNumber &b);

    friend FractionNumber operator*(const FractionNumber &a, int64_t other);

    friend FractionNumber operator*(int64_t other, const FractionNumber &a);

    void operator*=(const FractionNumber &other) {
        numerator *= other.numerator;
        denominator *= other.denominator;
        reduction();
    }

    void operator*=(int64_t other) {
        numerator *= other;
        reduction();
    }

    friend FractionNumber operator/(const FractionNumber &a, const FractionNumber &b);

    friend FractionNumber operator/(const FractionNumber &a, int64_t other);

    friend FractionNumber operator/(int64_t other, const FractionNumber &a);

    void operator/=(FractionNumber &other) {
        if (other.numerator == 0) {
            throw std::logic_error("Division by 0 in FractionNumber");
        }
        numerator *= other.denominator;
        denominator *= other.numerator;
        if (denominator < 0) {
            denominator *= -1;
            numerator *= -1;
        }
        reduction();
    }

    void operator/=(int64_t other) {
        if (other == 0) {
            throw std::logic_error("Division by 0 in FractionNumber");
        }
        denominator *= other;
        if (denominator < 0) {
            denominator *= -1;
            numerator *= -1;
        }
        reduction();
    }

    void reduction() {
        int64_t gcd = std::gcd(numerator, denominator);
        if (gcd > 1) {
            numerator /= gcd;
            denominator /= gcd;
        }
    }

    std::string str() const {
        if (denominator == 1) {
            return std::to_string(numerator);
        }
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }
};

FractionNumber
operator+(const FractionNumber &a, const FractionNumber &b) {
    FractionNumber num;
    int64_t lcm = std::lcm(a.denominator, b.denominator);
    num.numerator = a.numerator * (lcm / a.denominator) + b.numerator * (lcm / b.denominator);
    num.denominator = lcm;
    num.reduction();
    return num;
}

FractionNumber operator+(const FractionNumber &a, int64_t other) {
    FractionNumber num = a;
    num.numerator += other * num.denominator;
    num.reduction();
    return num;
}

FractionNumber operator+(int64_t other, const FractionNumber &a) {
    return a + other;
}

FractionNumber
operator-(const FractionNumber &a, const FractionNumber &b) {
    FractionNumber num;
    int64_t lcm = std::lcm(a.denominator, b.denominator);
    num.numerator = a.numerator * (lcm / a.denominator) - b.numerator * (lcm / b.denominator);
    num.denominator = lcm;
    num.reduction();
    return num;
}

FractionNumber operator-(const FractionNumber &a, int64_t other) {
    FractionNumber num = a;
    num.numerator -= other * num.denominator;
    num.reduction();
    return num;
}

FractionNumber operator-(int64_t other, const FractionNumber &a) {
    FractionNumber num(other);
    return num - a;
}

FractionNumber
operator*(const FractionNumber &a, const FractionNumber &b) {
    FractionNumber num = a;
    num.numerator *= b.numerator;
    num.denominator *= b.denominator;
    num.reduction();
    return num;
}

FractionNumber operator*(const FractionNumber &a, int64_t other) {
    FractionNumber num = a;
    num.numerator *= other;
    num.reduction();
    return num;
}

FractionNumber operator*(int64_t other, const FractionNumber &a) {
    return a * other;
}

FractionNumber
operator/(const FractionNumber &a, const FractionNumber &b) {
    if (b.numerator == 0) {
        throw std::logic_error("Division by 0 in FractionNumber");
    }
    FractionNumber num = a;
    num.numerator *= b.denominator;
    num.denominator *= b.numerator;
    if (num.denominator < 0) {
        num.denominator *= -1;
        num.numerator *= -1;
    }
    num.reduction();
    return num;
}

FractionNumber operator/(const FractionNumber &a, int64_t other) {
    if (other == 0) {
        throw std::logic_error("Division by 0 in FractionNumber");
    }
    FractionNumber num = a;
    num.denominator *= other;
    if (num.denominator < 0) {
        num.denominator *= -1;
        num.numerator *= -1;
    }
    num.reduction();
    return num;
}

FractionNumber operator/(int64_t other, const FractionNumber &a) {
    FractionNumber num(other);
    return num / a;
}

template<>
struct fmt::formatter<FractionNumber> {
    template<class ParseContext>
    constexpr auto parse(ParseContext &ctx) {
        return ctx.begin();
    }

    template<class FormatContext>
    auto format(FractionNumber const &number, FormatContext &ctx) {
        return fmt::format_to(ctx.out(), "{0}", number.str());
    }

};
