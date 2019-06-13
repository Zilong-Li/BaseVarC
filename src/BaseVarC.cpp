#include <getopt.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>

#include "htslib/bgzf.h"
#include "RefReader.h"
#include "BamProcess.h"
#include "BaseType.h"
#include "ThreadPool.h"

static const char* BASEVARC_USAGE_MESSAGE = 
"Program: BaseVarC -- A c version of BaseVar\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC <command> [options]\n\n"
"Commands:\n"
"           basetype       Variants Caller\n"
"           popmatrix      Create population matrix at specific positions.\n" 
"           concat         Concat popmatrix.\n" 
"           merge          Merge vcf/cvg files.\n" 
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* BASETYPE_MESSAGE = 
"Program: BaseVarC basetype\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC basetype [options]\n\n"
"Commands:\n"
"  --input,      -i        BAM/CRAM files list, one file per row.\n"
"  --output,     -o        Output filename prefix\n"
"  --reference,  -r        Reference file\n"
"  --region,     -s        Samtools-like region\n"
"  --group,      -g        Population group information\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --thread,     -t <INT>  Number of thread\n"
"  --batch,      -b <INT>  Number of samples each batch\n"
"  --load,                 Load data only\n"
"  --rerun,                Read previous loaded data and rerun\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --input,      -i        BAM/CRAM files list, one file per row.\n"
"  --posfile,    -p        Position file <CHRID POS REF ALT>\n"
"  --output,     -o        Output filename prefix(.mat.gz will be added auto)\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* CONCAT_MESSAGE =
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC concat [options]\n\n"
"Commands:\n"
"  --input,      -i       List of matrix files for concat, one file per row.\n"
"  --output,     -o       Output filename prefix(.gz will be added auto)\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* CVG_HEADER =
"##fileformat=CVGv1.0\n"
"##Group information is the depth of A:C:G:T\n"
"#CHROM\tPOS\tREF\tDepth\tA\tC\tG\tT\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)";

static const char* VCF_HEADER =
"##fileformat=VCFv4.2\n"
"##FILTER=<ID=LowQual,Description=\"Low quality (QUAL < 60)\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"##FORMAT=<ID=AB,Number=1,Type=String,Description=\"Allele Base\">\n"
"##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Strand orientation of the mapping base. Marked as + or -\">\n"
"##FORMAT=<ID=BP,Number=1,Type=String,Description=\"Base Probability which calculate by base quality\">\n"
"##INFO=<ID=CM_AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies base on LRT algorithm\">\n"
"##INFO=<ID=CM_CAF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies just base on read count\">\n"
"##INFO=<ID=CM_AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited allele depth in CMDB\">\n"
"##INFO=<ID=CM_DP,Number=A,Type=Integer,Description=\"Total Depth\">\n"
"##INFO=<ID=SB_REF,Number=A,Type=Integer,Description=\"Read number support REF: Forward,Reverse\">\n"
"##INFO=<ID=SB_ALT,Number=A,Type=Integer,Description=\"Read number support ALT: Forward,Reverse\">\n"
"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
"##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">\n"
"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Phred-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n"
"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence Quality by Depth\">\n";

typedef std::string String;
typedef std::vector<String> StringV;
typedef std::vector<int32_t> IntV;
typedef std::vector<PosAlleleMap> PosAlleleMapVec;
typedef std::map<String, IntV> GroupIdx;

struct BtRes
{
    String  cvg;
    String  vcf;
};

void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void runConcat(int argc, char **argv);
void parseOptions(int argc, char **argv, const char* msg);

void bt_read(StringV bams, const String& region, const IntV& pv, String fout);
BtRes bt_f(int32_t p, const GroupIdx& popg_idx, const AlleleInfoVector& aiv, const DepM& idx, int32_t N, const String& chr, int32_t rg_s, const String& seq);
// String bt_group(const GroupIdx& popg_idx, const AlleleInfoVector& aiv, const DepM& idx, const BaseV& base);

namespace opt {
    static bool verbose = false;
    static bool rerun   = false;
    static bool load    = false;
    static int mapq;
    static int thread;
    static int batch;
    static std::string input;
    static std::string reference;
    static std::string posfile;
    static std::string group;
    static std::string region;
    static std::string output;
}

static const char* shortopts = "hvi:r:p:s:o:q:t:b:g:";

static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "rerun",                   no_argument, NULL,  8  },
  { "load",                    no_argument, NULL,  7  },
  { "input",                   required_argument, NULL, 'i' },
  { "reference",               required_argument, NULL, 'r' },
  { "posfile",                 required_argument, NULL, 'p' },
  { "group",                   required_argument, NULL, 'g' },
  { "region",                  required_argument, NULL, 's' },
  { "output",                  required_argument, NULL, 'o' },
  { "batch",                   required_argument, NULL, 'b' },
  { "thread",                  required_argument, NULL, 't' },
  { "mapq",                    required_argument, NULL, 'q' },
  { NULL, 0, NULL, 0 }
};


int main(int argc, char** argv)
{
    if (argc <= 1 ) {
        std::cerr << BASEVARC_USAGE_MESSAGE;
        return 0;
    } else {
        std::string command(argv[1]);
        if (command == "help" || command == "--help") {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        } else if (command == "basetype") {
            runBaseType(argc - 1, argv + 1);
        } else if (command == "popmatrix") {
            runPopMatrix(argc - 1, argv + 1);
        } else if (command == "concat") {
            runConcat(argc - 1, argv + 1);
        } else {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        }
    }

    return 0;
}

void runBaseType(int argc, char **argv)
{
    parseOptions(argc, argv, BASETYPE_MESSAGE);
    std::cerr << "basetype start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    StringV bams(std::istream_iterator<BaseVar::Line>{ibam},
    	         std::istream_iterator<BaseVar::Line>{});
    const int32_t N = bams.size();
    // fetch ref bases;
    if (opt::region.empty()) {
        std::cerr << "Error: region must be specified" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (opt::reference.empty()) {
        std::cerr << "Error: reference must be specified" << std::endl;
        exit(EXIT_FAILURE);
    }
    RefReader fa;
    String seq = fa.GetTargetBase(opt::region, opt::reference);
    String chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(opt::region);
    IntV pv;
    for (size_t i = 0; i < seq.length(); ++i) {
        if (seq[i] == 'N') continue;
        pv.push_back(i + rg_s);       // 1-based
    }
    // begin to read bams
    String headcvg = String(CVG_HEADER);
    String headvcf = String(VCF_HEADER);
    String tmp;
    StringV ftmp_v;
    int bc = opt::batch;
    int bt = 1 + (N - 1) / bc;    // ceiling
    if (!opt::rerun) {
        int thread;
        if (opt::thread > 1) thread = opt::thread;
        else {
            std::cerr << "threads must be larger than 1" << std::endl;
            exit(EXIT_FAILURE);
        }
        // @todo: optimize threadpool. current cpu usage less than 50%
        ThreadPool pool(thread);
        std::vector<std::future<void>> res;
        String region = opt::region;
        StringV bams_t;
        std::cerr << "begin to read bams and save as tmp file" << std::endl;
        for (int i = 0; i < bt; ++i) {
            tmp = opt::output + ".batch." + BaseVar::tostring(i) + ".tmp";
            ftmp_v.push_back(tmp);
            if (i == bt - 1) {
                StringV t(bams.begin() + i * bc, bams.end());
                bams_t = t;
            } else {
                StringV t(bams.begin() + i * bc, bams.begin() + (i + 1)*bc);
                bams_t = t;
            }
            res.emplace_back(pool.enqueue(bt_read, bams_t, std::cref(region), std::cref(pv), tmp));
        }
        for (auto && r: res) {
            r.get();
        }
        res.clear();
        if (opt::load) exit(EXIT_SUCCESS);
    } else {
        for (int i = 0; i < bt; ++i) {
            tmp = opt::output + ".batch." + BaseVar::tostring(i) + ".tmp";
            ftmp_v.push_back(tmp);
        }
    }
    // begin to call basetype and output
    String vcfout = opt::output + ".vcf.gz";
    String cvgout = opt::output + ".cvg.gz";
    BGZF* fpv = bgzf_open(vcfout.c_str(), "w");
    BGZF* fpc = bgzf_open(cvgout.c_str(), "w");
    BGZF* fpi;
    std::vector<BGZF*> fpiv;
    for (auto & f: ftmp_v) {
        fpi = bgzf_open(f.c_str(), "r");
        fpiv.push_back(fpi);
    }
    kstring_t ks = {0, 0, NULL};
    String sams;
    for (auto & fp: fpiv) {
        if (bgzf_getline(fp, '\n', &ks) >= 0) {
            sams += (String)ks.s;
        }
    }
    sams.pop_back();
    // fetch popgroup information
    GroupIdx popg_idx;
    if (!opt::group.empty()) {
        std::ifstream ifg(opt::group);
        String id, grp;
        std::unordered_map<String, String> popg_m;
        while (ifg >> id >> grp) {
            popg_m.insert({id, grp});
        }
        std::istringstream iss(sams);
        int i = -1;
        while (std::getline(iss, id, '\t')) {
            i++;
            if (popg_m.count(id)) {
                grp = popg_m.at(id);
            } else {
                continue;
            }
            if (popg_idx.count(grp)) {
                popg_idx[grp].push_back(i);
            } else {
                IntV v{i};
                popg_idx.insert({grp, v});
            }
        }
        if (!popg_idx.empty()) {
            for (GroupIdx::iterator it = popg_idx.begin(); it != popg_idx.end(); ++it) {
                headcvg += "\t" + it->first;
                headvcf += "##INFO=<ID=" + it->first + "_AF,Number=A,Type=Float,Description=\"Allele frequency in the " + it->first + " populations calculated based on LRT.[0,1]\">\n";
            }
        }
    }
    // output header
    headcvg += "\n";
    if (bgzf_write(fpc, headcvg.c_str(), headcvg.length()) != headcvg.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // get contig from fai file.
    String fai = opt::reference + ".fai";
    std::ifstream ifai(fai);
    if (ifai.is_open()) {
        String contig, len, t;
        while (ifai >> contig >> len >> t >> t >> t) {
            headvcf += "##contig=<ID=" + contig + ",length=" + len + ",assembly=" + opt::reference + ">\n";
        }
    }
    headvcf += "#reference=file://" + opt::reference + "\n";
    headvcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sams + "\n";
    if (bgzf_write(fpv, headvcf.c_str(), headvcf.length()) != headvcf.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    AlleleInfo ai;
    AlleleInfoVector aiv;
    DepM idx;
    int32_t j, k, i, count=0;
    IntV pv_t;
    std::cerr << "begin to load data and run basetype" << std::endl;
    char *buf=NULL, *str=NULL, *str2=NULL, *pti=NULL, *pto=NULL;
    for (auto & p : pv) {
        j = 0; k = 0;
        for (auto & fp: fpiv) {
            if (bgzf_getline(fp, '\n', &ks) >= 0) {
                buf = ks.s;
                while ((str = strtok_r(buf, " ", &pto)) != NULL) {
                    if (str[0] != '.') {
                        buf = str;
                        for (i = 0; i < 5; ++i) {
                            if ((str2 = strtok_r(buf, ",", &pti)) != NULL) {
                                switch(i){
                                case 0: ai.base = std::atoi(str2);break;
                                case 1: ai.mapq = std::atoi(str2);break;
                                case 2: ai.qual = std::atoi(str2);break;
                                case 3: ai.rpr = std::atoi(str2);break;
                                case 4: ai.strand = std::atoi(str2);break;
                                }
                                buf = NULL;
                            }
                        }
                        // skip N base
                        if (ai.base != 4) {
                            aiv.push_back(ai);
                            idx.insert({j, k++});
                        }
                    }
                    buf = NULL;
                    j++;
                }
            }
        }
        if (aiv.size() > 0) {
            auto btr = bt_f(p, popg_idx, aiv, idx, N, chr, rg_s, seq);
            if (btr.vcf != "NA" && bgzf_write(fpv, btr.vcf.c_str(), btr.vcf.length()) != btr.vcf.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (bgzf_write(fpc, btr.cvg.c_str(), btr.cvg.length()) != btr.cvg.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
            aiv.clear();
            idx.clear();
            if (!(++count % 1000)) std::cerr << "basetype completed " << count << " sites" << std::endl;
        }
    }
    if (bgzf_close(fpv) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    if (bgzf_close(fpc) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    // remove tmp file
    for (auto & f: ftmp_v) {
        std::remove(f.c_str());
    }
    // done
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cerr << "basetype done" << std::endl;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
}

void bt_read(StringV bams, const String& region, const IntV& pv, String fout)
{
    char space = ' ';
    char comma = ',';
    size_t n = bams.size();
    PosAlleleMapVec allele_mv;
    String names;
    allele_mv.reserve(n);
    int32_t count = 0;
    for (auto const& bam: bams) {
        BamProcess reader;
        if (!(++count % 100)) std::cerr << "reading number " << count << " bam -- " << fout << std::endl;
        if (!reader.Open(bam)) {
            std::cerr << "ERROR: could not open file " << bam << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.FindSnpAtPos(region, pv)) {
            std::cerr << "Warning: " << reader.sm << " region " << region << " is empty." << std::endl;
        }
        allele_mv.push_back(reader.allele_m);
        names += reader.sm + '\t';
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bam << std::endl;
        }
    }
    names += "\n";
    std::ostringstream out;
    out << names;
    BGZF* fp = bgzf_open(fout.c_str(), "w");
    if (bgzf_write(fp, out.str().c_str(), out.str().length()) != out.str().length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    out.str("");
    out.clear();
    for (auto & p : pv) {
        for (auto const& m : allele_mv) {
            if (m.count(p)) {
                auto const& a = m.at(p);
                out << a.base << comma << a.mapq << comma << a.qual << comma << a.rpr << comma << a.strand << space;
            } else {
                out << ". ";
            }
        }
        out << "\n";
        if (bgzf_write(fp, out.str().c_str(), out.str().length()) != out.str().length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        out.str("");
        out.clear();
    }
    if (bgzf_close(fp) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
}

BtRes bt_f(int32_t p, const GroupIdx& popg_idx, const AlleleInfoVector& aiv, const DepM& idx, int32_t N, const String& chr, int32_t rg_s, const String& seq)
{
    int8_t alt_base, ref_base;
    int32_t dep, na, nc, ng, nt, ref_fwd, ref_rev, alt_fwd, alt_rev;
    double fs, sor, left_p, right_p, twoside_p;
    double min_af = 0.001;
    BaseV bases, quals;
    std::ostringstream oss;
    BtRes res;
    // output cvg;
    ref_base = BASE_INT8_TABLE[static_cast<int>(seq[p - rg_s])];
    na = 0; nc = 0; ng = 0; nt = 0;
    for (auto const& a: aiv) {
        switch (a.base) {
        case 0 : na += 1; break;
        case 1 : nc += 1; break;
        case 2 : ng += 1; break;
        case 3 : nt += 1; break;
        }
        bases.push_back(a.base);
        quals.push_back(a.qual);
    }
    IntV tmp = {na, nc, ng, nt};
    std::vector<size_t> didx = BaseVar::sortidx(tmp);
    if (didx[0] != ref_base) {
        alt_base = didx[0];
    } else{
        alt_base = didx[1];
    }
    ref_fwd = 0, ref_rev = 0, alt_fwd = 0, alt_rev = 0;
    for (auto const& a: aiv) {
        if (a.strand == 1) {
            if (a.base == ref_base) {
                ref_fwd += 1;
            } else if (a.base == alt_base) {
                alt_fwd += 1;
            }
        } else if (a.strand == 0) {
            if (a.base == ref_base) {
                ref_rev += 1;
            } else if (a.base == alt_base) {
                alt_rev += 1;
            }
        }
    }
    kt_fisher_exact(ref_fwd, ref_rev, alt_fwd, alt_rev, &left_p, &right_p, &twoside_p);
    fs = -10 * log10(twoside_p);
    if (std::isinf(fs)) fs = 10000.0;
    else if (fs == 0) fs = 0.0;
    if (alt_fwd * ref_rev > 0) {
        sor = static_cast<double>(ref_fwd * alt_rev) / (ref_rev * alt_fwd);
    } else {
        sor = 10000.0;
    }
    dep = na + nc + ng + nt;
    oss << chr << '\t' << p << '\t' << BASE2CHAR[ref_base] << '\t' << dep << '\t' << na << '\t' << nc << '\t' << ng << '\t' << nt << '\t' << fs << '\t' << sor << '\t' << ref_fwd << ',' << ref_rev << ',' << alt_fwd << ',' << alt_rev << '\t';
    // basetype caller;
    BaseType bt(bases, quals, ref_base, min_af);
    const bool bt_success = bt.LRT();
    BaseV base_comb{ref_base};
    base_comb.insert(base_comb.end(), bt.alt_bases.begin(), bt.alt_bases.end());
    // popgroup depth
    InfoM info;
    if (!popg_idx.empty()) {
        BaseV gr_bases, gr_quals;
        String gr_af;
        for (GroupIdx::const_iterator it = popg_idx.begin(); it != popg_idx.end(); ++it) {
            na = 0; nc = 0; ng = 0; nt = 0;
            for (auto i : it->second) {
                if (idx.count(i)) {
                    auto & a = aiv[idx.at(i)];
                    switch (a.base) {
                    case 0 : na += 1; break;
                    case 1 : nc += 1; break;
                    case 2 : ng += 1; break;
                    case 3 : nt += 1; break;
                    }
                    if (bt_success) {
                        gr_bases.push_back(a.base);
                        gr_quals.push_back(a.qual);
                    }
                }
            }
            oss << na << ':' << nc << ':' << ng << ':' << nt << '\t';
            if (bt_success) {
                // todo: check out the result
                BaseType gr_bt(gr_bases, gr_quals, ref_base, min_af);
                gr_bt.SetBase(base_comb);
                gr_bt.LRT();
                gr_af = "";
                for (auto b : bt.alt_bases) {
                    if (gr_bt.af_lrt.count(b)) {
                        gr_af += BaseVar::tostring(gr_bt.af_lrt[b]) + ",";
                    } else {
                        gr_af += "0,";
                    }
                }
                gr_af.pop_back();
                gr_bases.clear();
                gr_quals.clear();
                info.insert({it->first + "_AF", gr_af});
            }
        }
        res.cvg = oss.str(); res.cvg.pop_back(); res.cvg += "\n";
    } else {
        res.cvg = oss.str(); res.cvg += "\n";
    }
    if (bt_success) {
        res.vcf = bt.WriteVcf(bt, chr, p, ref_base, aiv, idx, info, N);
    } else {
        res.vcf = "NA";
    }

    return res;
}

void runPopMatrix (int argc, char **argv)
{
    parseOptions(argc, argv, POPMATRIX_MESSAGE);
    std::cerr << "popmatrix start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    std::ifstream ipos(opt::posfile);
    if (!ibam.is_open() || !ipos.is_open()) {
        std::cerr << "bamlist or posifle can not be opend" << std::endl;
        exit(EXIT_FAILURE);
    }
    StringV bams(std::istream_iterator<BaseVar::Line>{ibam},
    	         std::istream_iterator<BaseVar::Line>{});
    PosInfoVector pv;
    for (PosInfo p; ipos >> p;) pv.push_back(p);
    pv.shrink_to_fit();        // request for the excess capacity to be released

    const int32_t N = bams.size();
    const int32_t M = pv.size();
    String out = BaseVar::tostring(N) + "\t" + BaseVar::tostring(M) + "\n";
    BGZF* fp = bgzf_open(opt::output.c_str(), "w");
    if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // ready for run
    String rg = pv.front() + pv.back();
    int32_t count = 0;
    for (int32_t i = 0; i < N; i++) {
        BamProcess reader;
        if (!(++count % 1000)) std::cerr << "Processing the number " << count / 1000 << "k bam" << std::endl;
        if (!reader.Open(bams[i])) {
            std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
            exit(EXIT_FAILURE);
        }
        reader.FindSnpAtPos(rg, pv);
        String out(reader.snps.begin(), reader.snps.end());
        out += "\n";
        if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bams[i] << std::endl;
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "failed to close \n";
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cerr << "popmatrix done" << std::endl;
}

void runConcat(int argc, char **argv)
{
    parseOptions(argc, argv, CONCAT_MESSAGE);
    String fm = opt::input;
    String fo = opt::output;
    clock_t ctb = clock();
    std::ifstream ifm(fm);
    StringV fm_v(std::istream_iterator<BaseVar::Line>{ifm},
	             std::istream_iterator<BaseVar::Line>{});
    int32_t i, k, m, count, n = 0, mt = 0;
    int32_t nm = fm_v.size();
    BGZF* fp = NULL;
    kstring_t ks = {0, 0, NULL};
    String ss;
    StringV Drg;
    std::vector<StringV> D;
    D.reserve(nm);
    for (k = 0; k < nm ; ++k) {
        std::cout << "reading file : " << fm_v[k] << std::endl;
        fp = bgzf_open(fm_v[k].c_str(), "r");
        if (bgzf_getline(fp, '\n', &ks) >= 0) {
            ss = ks.s;
            StringV tokens;
            String token;
            std::istringstream ts(ss);
            while (std::getline(ts, token, '\t')) tokens.push_back(token);
            if (n > 0 && n != std::stol(tokens[0])) {
                std::cerr << "error: the number of samples among the inputs are different!" << std::endl;
                exit(EXIT_FAILURE);
            } else {
                n = std::stol(tokens[0]);
                m = std::stol(tokens[1]);
                mt += m;
            }
        }
        count = 0;
        Drg.reserve(n);
        while (bgzf_getline(fp, '\n', &ks) >= 0) {
            ++count;
            Drg.push_back(ks.s);
        }
        if (count != n) {
            std::cerr << "error: the number of samples dont match the header!" << std::endl;
            exit(EXIT_FAILURE);
        }
        D.push_back(Drg);
        Drg.clear();
        if (bgzf_close(fp) < 0) {
            std::cerr << "warning: file cannot be closed!" << std::endl;
        }
    }

    fp = bgzf_open(fo.c_str(), "w");
    ss = BaseVar::tostring(n) + "\t" + BaseVar::tostring(mt) + "\n";
    if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // fast concatenation consider using rope data structure
    // @see https://brianbondy.com/blog/tagged/data-structure
    for (i = 0; i < n; ++i){
        ss = "";
        ss.reserve(mt+1);
        for (k = 0; k < nm; ++k) {
            ss += D[k][i];
        }
        ss += "\n";
        if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
            std::cerr << "fail to write" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cout << "concat done" << std::endl;
}

void parseOptions(int argc, char **argv, const char* msg)
{
    bool die = false;
    bool help = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case 'q': arg >> opt::mapq; break;
        case 'b': arg >> opt::batch; break;
        case 't': arg >> opt::thread; break;
        case 'i': arg >> opt::input; break;
        case 'r': arg >> opt::reference; break;
        case 'p': arg >> opt::posfile; break;
        case 's': arg >> opt::region; break;
        case 'g': arg >> opt::group; break;
        case 'o': arg >> opt::output; break;
        case  8 : opt::rerun   = true; break;
        case  7 : opt::load    = true; break;
        case 'v': opt::verbose = true; break;
        default: die = true;
        }
    }
    // todo : need more check
    if (die || help || opt::input.empty() || opt::output.empty()) {
        std::cerr << msg;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}
