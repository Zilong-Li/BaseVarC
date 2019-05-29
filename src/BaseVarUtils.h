#ifndef BASEVAR_UTILS_H
#define BASEVAR_UTILS_H

#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <numeric>
#include <algorithm>

namespace BaseVar
{
struct Line
{
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
        return std::getline(is, line.data);
    }
};

template<typename T>
inline std::string tostring(T d) {
    std::stringstream ss;
    ss << d;
    return ss.str();
}

template<typename T>
std::vector<size_t> sortidx(const std::vector<T>& v) {
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}

inline std::tuple<std::string, int32_t, int32_t> splitrg(std::string rg) {
    size_t p;
    std::string chr;
    int32_t s = 0, e = 0;
    if ((p = rg.find(":")) != std::string::npos) {
        chr = rg.substr(0, p);
        rg.erase(0, p + 1);
    }
    if ((p = rg.find("-")) != std::string::npos) {
        s = std::stoi(rg.substr(0, p)) - 1;      // make 0-based
        rg.erase(0, p + 1);
        e = std::stoi(rg) - 1;
    }
    return std::make_tuple(chr, s, e);
}

}

#endif