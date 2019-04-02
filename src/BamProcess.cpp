#include "BamProcess.h"

void BamProcess::findSnpAtPos(const PosInfoVector& pv) {
    SeqLib::BamRecord r;
    /* check if the BAM is sorted */
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }

    while (GetNextRecord(r), r.MapQuality() < mapq);  // need check out
    int32_t pos = r.Position() + 1; // make 1-based
    bool flag = true;

    snps.reserve(pv.size());       // best practice;
    for (auto const& s: pv) {
	// select the first record covering this position
	if (s.pos < pos) {
	    snps.push_back('N'); 
	} else {
	    // r.PositionEnd is 1-based
	    while (s.pos > r.PositionEnd() && (flag = GetNextRecord(r))) {
		while (r.MapQuality() < mapq && flag) flag = GetNextRecord(r);
		SeqLib::Cigar c = r.GetCigar();
		// only choose the BAM_CMATCH record
		while ((c.size()!=1 || c[0].RawType() != BAM_CMATCH) && flag) {
		    flag = GetNextRecord(r);
		    c=r.GetCigar();
		}
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

inline char BamProcess::getSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const {
    int32_t offset = s.pos - (r.Position() + 1);
    // must copy r.r.Sequence();
    char seq[r.Sequence().length() + 1];
    std::strcpy(seq, r.Sequence().c_str());
    return seq[offset];
}
