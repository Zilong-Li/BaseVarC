#include "BamProcess.h"

void BamProcess::FindSnpAtPos(const SeqLib::GenomicRegion& gr, const PosInfoVector& pv) {

    SetRegion(gr);
    SeqLib::BamRecord r;
    /* check if the BAM is sorted */
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }

    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
	// filter reads by mapq and cigar type (only BAM_CMATCH);
	if (r.MapQuality() < mapq) continue;
	SeqLib::Cigar c = r.GetCigar();
	if (c.size()!=1 || c[0].RawType() != BAM_CMATCH) continue;
	rv.push_back(r);
    } 

    bool flag = false;
    unsigned int i = 0;
    r = rv[i];
    snps.reserve(pv.size());       // best practice;
    for (auto const& s: pv) {
	/* select the first record covering this position */
	if (s.pos < r.Position() + 1) { // make 1-based
	    snps.push_back('N');
	}else{
	     // r.PositionEnd is 1-based
	    while (s.pos > r.PositionEnd()) {
		if (i == rv.size() - 1) {flag = true; break;}
		i++;
		r = rv[i];
	    } 
	    if (s.pos < r.Position() + 1 || flag) {
		snps.push_back('N');
	    } else {
	 	snps.push_back(GetSnpCode(r, s));
	    }
	}
    }
}

inline char BamProcess::GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const {
    int32_t offset = s.pos - (r.Position() + 1);
    char seq[r.Sequence().length() + 1]; // must copy r.r.Sequence();
    std::strcpy(seq, r.Sequence().c_str());
    return seq[offset];
}
