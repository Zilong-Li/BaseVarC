#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>

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
	return snp.chr + colon + this->pos + hyphen + snp.pos;
    }
    // define an overload of operator >>
    friend std::istream& operator>>(std::istream& is, SnpInfo& info) {
	is >> std::ws >> info.chr >> info.pos >> info.ref >> info.alt;
	return is;
    }
};

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
    std::vector<SnpInfo> snps;
    for (SnpInfo info; ifs >> info;) {
	snps.push_back(info);
    }
    SnpInfo s = snps.front();
    SnpInfo e = snps.back();
    std::string r = s + e;
    std::cout << snps.size() << "\t" << s.pos << "\t" << e.pos << std::endl;
}