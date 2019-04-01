#include "BamProcess.h"

void BamProcess::findSnpAtPos(const PosInfoVector& pv) {
    SeqLib::BamRecord r;
    // check if the BAM is sorted
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
	std::cerr << "       Sorted BAMs are required." << std::endl;
	exit(EXIT_FAILURE);
    }

    while (GetNextRecord(r), r.MapQuality() < mapq) {};  // need check out
    int32_t pos = r.Position() + 1; // make 1-based
    bool flag = true;

    snps.reserve(pv.size());       // best practice;
    for (auto const& s: pv) {
	// select the first record covering this position
	if (s.pos < pos) {
	    snps.push_back('N'); 
	} else {
	    while (s.pos > r.PositionEnd() + 1 && (flag = GetNextRecord(r), flag)) {
		pos = r.Position() + 1;
	    }
	    if (s.pos < pos || !flag) {
		snps.push_back('N'); 
	    } else {
		snps.push_back(getSnpCode(r, s));
	    }
	}
    }
}

inline char BamProcess::getSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) {
    int32_t offset = s.pos - (r.Position() + 1);
    const char* seq = r.Sequence().c_str();
    return seq[offset - 1];
}
