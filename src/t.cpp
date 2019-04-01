#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>


namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

struct Line {
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
	return std::getline(is, line.data);
    }
};

struct SnpInfo {
    std::string chr;
    std::string pos;
    char ref;
    char alt;
    // provide an overload of operator + to return samtools-like region
    std::string operator+(const SnpInfo& snp) {
	assert(this->chr == snp.chr);
	const char *colon = ":";
	const char *hyphen = "-";
	auto rg_s = patch::to_string(this->pos);
	auto rg_e = patch::to_string(snp.pos);
	return snp.chr + colon + rg_s + hyphen + rg_e;
    }
    // define an overload of operator >>
    friend std::istream& operator>>(std::istream& is, SnpInfo& info) {
	is >> std::ws >> info.chr >> info.pos >> info.ref >> info.alt;
	return is;
    }
};
typedef std::vector<SnpInfo> SnpInfoVector;

int main(int argc, char** argv)
{
    std::string bamlst = argv[1];
    std::ifstream ifs1(bamlst);
    std::vector<std::string> files(std::istream_iterator<Line>{ifs1},
    	                           std::istream_iterator<Line>{});
    for (const auto& f: files) {
	std::cout << f << std::endl;
    }

    std::string input = argv[2];
    std::ifstream ifs(input);
    SnpInfoVector snps;
    for (SnpInfo info; ifs >> info;) {
	snps.push_back(info);
    }
    // SnpInfo s = snps.front();
    // SnpInfo e = snps.back();
    std::string r = snps.front() + snps.back();
    std::cout << snps.size() << "\t" << r << std::endl;
}