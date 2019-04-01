#ifndef __BAM_PROCESS_H__
#define __BAM_PROCESS_H__

#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/ReadFilter.h"

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

struct PosInfo {
    std::string chr;
    int pos;
    char ref;
    char alt;
    // provide an overload of operator + to return samtools-like region
    std::string operator+(const PosInfo& snp) const {
	assert(this->chr == snp.chr);
	const char *colon = ":";
	const char *hyphen = "-";
	auto rg_s = patch::to_string(this->pos);
	auto rg_e = patch::to_string(snp.pos);
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
    void findSnpAtPos(const PosInfoVector& pv);
    int32_t mapq = 10;
    std::vector<char> snps;

 private:
    inline char getSnpCode(const SeqLib::BamRecord& r, const PosInfo& s);
};

// class PopMatrix {

//  public:
//     void subPopMatrix(const std::vector<std::string>& bams, const PosInfoVector& pv);

// };

// void PopMatrix::subPopMatrix(const std::vector<std::string>& bams, const PosInfoVector& pv) {
    
// }

#endif