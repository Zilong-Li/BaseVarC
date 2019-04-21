#ifndef BASEVAR_UTILS_H
#define BASEVAR_UTILS_H

#include <sstream>
#include <string>

namespace BaseVar
{

template<typename T> 
inline std::string tostring(T d) { 
    std::stringstream ss;
    ss << d;
    return ss.str();
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