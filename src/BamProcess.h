#ifndef __BASEVARC_BAM_PROCESS_H__
#define __BASEVARC_BAM_PROCESS_H__

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "SeqLib/SeqLibUtils.h"

struct Line {
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
	return std::getline(is, line.data);
    }
};

struct PosInfo {
    std::string chr;
    int pos;
    char ref;
    char alt;
    /* provide an overload of operator + to return samtools-like region */
    std::string operator+(const PosInfo& snp) const {
	assert(this->chr == snp.chr);
	const char *colon = ":";
	const char *hyphen = "-";
	auto rg_s = SeqLib::tostring(this->pos);
	auto rg_e = SeqLib::tostring(snp.pos);
	return snp.chr + colon + rg_s + hyphen + rg_e;
    }
    friend std::istream& operator>>(std::istream& is, PosInfo& info) {
	is >> std::ws >> info.chr >> info.pos >> info.ref >> info.alt;
	return is;
    }
};
typedef std::vector<PosInfo> PosInfoVector;

class BamProcess: public SeqLib:: BamReader {

 public:
    BamProcess(){}
    void findSnpAtPos(const SeqLib::GenomicRegion& gr, const PosInfoVector& pv);
    int32_t mapq = 10;
    std::vector<char> snps;

 private:
    inline char getSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const;
    // std::unordered_map<char, int> snpCode{ {'N', 0},{'A', 1},{'C', 2},{'G', 3},{'T', 4} };
};


#endif