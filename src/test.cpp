#include <string>
#include <vector>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>
#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"

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

int main(int argc, char* argv[])
{

    std::string bamlst = argv[1];
    std::ifstream ifs1(bamlst);
    std::vector<std::string> files(std::istream_iterator<Line>{ifs1},
    	                           std::istream_iterator<Line>{});

    std::string posfile = argv[2];
    std::ifstream ifs2(posfile);
    std::vector<SnpInfo> snps;
    for (SnpInfo info; ifs2 >> info;) {
    	snps.push_back(info);
    }
    snps.shrink_to_fit();        // request for the excess capacity to be released
    SnpInfo s = snps.front();
    SnpInfo e = snps.back();
    std::string rg = s + e;

    SeqLib::BamReader r;
    r.Open(files);
    SeqLib::GenomicRegion gr(rg, r.Header());
    r.SetRegion(gr);

    SeqLib::BamRecord rec;
    while(r.GetNextRecord(rec)) {
    	std::cout << rec.Position() << std::endl;
    }
    
}
