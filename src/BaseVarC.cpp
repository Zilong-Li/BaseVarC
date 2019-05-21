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
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
"  --output,     -o        Output filename prefix\n"
"  --reference,  -r        Reference file\n"
"  --region,     -g        Samtools-like region\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --thread,     -t <INT>  Number of thread\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
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
"  --input,      -l       List of matrix files for concat, one file per row.\n"
"  --output,     -o       Output filename prefix(.gz will be added auto)\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* CVG_HEADER =
"##fileformat=CVGv1.0\n"
"##Group information is the depth of A:C:G:T\n"
"#CHROM\tPOS\tREF\tDepth\tA\tC\tG\tT\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)\n";

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

void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void runConcat(int argc, char **argv);
void parseOptions(int argc, char **argv, const char* msg);

void bt_f(std::shared_ptr<std::ofstream>& fpc, std::shared_ptr<std::ofstream>& fpv, const IntV& pv, const std::vector<PosAlleleMap>& allele_mv, int32_t N, const String& chr, int32_t rg_s, const String& seq);
void bt_read(const std::vector<String>& bams, PosAlleleMapVec& allele_mv, const String& region, const IntV& pv, String& headvcf);

namespace opt {
    static bool verbose = false;
    static int mapq;
    static int thread;
    static std::string input;
    static std::string reference;
    static std::string posfile;
    static std::string region;
    static std::string output;
}

static const char* shortopts = "hvl:r:p:g:o:q:t:";

static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "input",                   required_argument, NULL, 'l' },
  { "reference",               required_argument, NULL, 'r' },
  { "posfile",                 required_argument, NULL, 'p' },
  { "region",                  required_argument, NULL, 'g' },
  { "output",                  required_argument, NULL, 'o' },
  { "thread",                  required_argument, NULL, 't' },
  { "mapq",                    required_argument, NULL, 'q' },
  { NULL, 0, NULL, 0 }
};

std::mutex mut;

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
    kstring_t ks;
    ks.s = 0; ks.l = 0; ks.m = 0;
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
    ss = std::to_string(n) + "\t" + std::to_string(mt) + "\n";
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

void runBaseType(int argc, char **argv)
{
    parseOptions(argc, argv, BASETYPE_MESSAGE);
    std::cerr << "basetype start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    StringV bams(std::istream_iterator<BaseVar::Line>{ibam},
    	         std::istream_iterator<BaseVar::Line>{});
    RefReader fa;
    String seq = fa.GetTargetBase(opt::region, opt::reference);
    String chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(opt::region);
    IntV pv;
    for (size_t i = 0; i < seq.length(); i++) {
        if (seq[i] == 'N') continue;
        pv.push_back(i + rg_s);       // 1-based
    }
    String headcvg = String(CVG_HEADER);
    String headvcf = String(VCF_HEADER) + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    int nt = opt::thread;
    const int32_t N = bams.size();
    int32_t step = N / nt;
    std::vector<StringV> bams_v;
    for (int i = 0; i < nt; ++i) {
        if (i == nt - 1){
            StringV t(bams.begin() + i * step, bams.end());
            bams_v.push_back(t);
        } else {
            StringV t(bams.begin() + i * step, bams.begin() + (i + 1)*step);
            bams_v.push_back(t);
        }
    }
    PosAlleleMapVec allele_mv;
    allele_mv.reserve(N);
    ThreadPool pool(nt);
    std::vector<std::future<void>> res(nt);
    for (int i = 0; i < nt; ++i) {
        res[i] = pool.enqueue(bt_read, std::cref(bams_v[i]), std::ref(allele_mv), std::cref(opt::region), std::cref(pv), std::ref(headvcf));
    }
    for (int i = 0; i < nt; ++i) {
        res[i].get();
    }
    bams_v.clear();

    headvcf += "\n";
    assert(allele_mv.size() == N);
    std::vector<std::shared_ptr<std::ofstream> > fpvv;
    std::vector<std::shared_ptr<std::ofstream> > fpcv;
    StringV out_cv;
    StringV out_vv;
    for (int i = 0; i < nt; ++i) {
        String vcfout = opt::output + ".tmp." + std::to_string(i) + ".vcf";
        String cvgout = opt::output + ".tmp." + std::to_string(i) + ".cvg";
        out_vv.push_back(vcfout);
        out_cv.push_back(cvgout);
        fpvv.push_back(std::make_shared<std::ofstream>(vcfout));
        fpcv.push_back(std::make_shared<std::ofstream>(cvgout));
    }
    step = pv.size() / nt;
    std::vector<IntV> pp;
    pp.reserve(nt);
    for (int i = 0; i < nt; ++i) {
        if (i == nt - 1){
            IntV t(pv.begin() + i * step, pv.end());
            pp.push_back(t);
        } else {
            IntV t(pv.begin() + i * step, pv.begin() + (i + 1)*step);
            pp.push_back(t);
        }
    }
    pv.clear();
    for (int i = 0; i < nt; ++i) {
        res[i] = pool.enqueue(bt_f, std::ref(fpcv[i]), std::ref(fpvv[i]), std::cref(pp[i]), std::cref(allele_mv), N, std::cref(chr), rg_s, std::cref(seq));
    }
    for (int i = 0; i < nt; ++i) {
        res[i].get();
    }
    fpvv.clear();
    fpcv.clear();
    // merge output files
    String vcfout = opt::output + ".vcf.gz";
    String cvgout = opt::output + ".cvg.gz";
    BGZF* fpv = bgzf_open(vcfout.c_str(), "w");
    BGZF* fpc = bgzf_open(cvgout.c_str(), "w");
    String line;
    for (int i = 0; i < nt; ++i) {
        std::ifstream fv(out_vv[i]);
        while (std::getline(fv, line)) {
            line += "\n";
            if (bgzf_write(fpv, line.c_str(), line.length()) != line.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::remove(out_vv[i].c_str());
        std::ifstream fc(out_cv[i]);
        while (std::getline(fc, line)) {
            line += "\n";
            if (bgzf_write(fpc, line.c_str(), line.length()) != line.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::remove(out_cv[i].c_str());
    }
    if (bgzf_close(fpv) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    if (bgzf_close(fpc) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    // done
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cerr << "basetype done" << std::endl;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
}

void bt_read(const std::vector<String>& bams, PosAlleleMapVec& allele_mv, const String& region, const IntV& pv, String& headvcf)
{
    for (size_t i = 0; i < bams.size(); ++i) {
        BamProcess reader;
        if (!reader.Open(bams[i])) {
            std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.FindSnpAtPos(region, pv)) {
            std::cerr << "Warning: " << reader.sm << " region " << region << " is empty." << std::endl;
        }
        std::unique_lock<std::mutex> guard(mut);
        headvcf += "\t" + reader.sm;
        allele_mv.push_back(reader.allele_m);
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bams[i] << std::endl;
        }
    }
}

void bt_f(std::shared_ptr<std::ofstream>& fpc, std::shared_ptr<std::ofstream>& fpv, const IntV& pv, const std::vector<PosAlleleMap>& allele_mv, int32_t N, const String& chr, int32_t rg_s, const String& seq)
{
    char col = ';';
    char tab = '\t';
    int8_t alt_base, ref_base;
    int32_t j = 0, na, nc, ng, nt, ref_fwd, ref_rev, alt_fwd, alt_rev;
    double fs, sor, left_p, right_p, twoside_p;
    double min_af = 0.001;
    BaseV bases, quals;
    AlleleInfoVector aiv;
    DepM idx;
    IntV tmp;
    std::ostringstream sout;
    String out;
    for (auto const& p : pv) {
        for (int32_t i = 0; i < N; ++i) {
            auto& m = allele_mv[i];
            if (m.count(p) == 0) {
                continue;
            } else if (m.at(p).base != 4) {
                // skip N base
                aiv.push_back(m.at(p));
                idx.insert({i, j++});
            }
        }
        // skip coverage==0
        if (aiv.size() > 0) {
            // here call BaseType
            ref_base = BASE_INT8_TABLE[static_cast<int>(seq[p - rg_s])];
            na = 0, nc = 0, ng = 0, nt = 0;
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
            // output cvg;
            tmp = {na, nc, ng, nt};
            std::sort(tmp.begin(), tmp.end(), std::greater<int>());
            if (tmp[0] != ref_base) {
                if (tmp[0] == na) alt_base = 0;
                else if (tmp[0] == nc) alt_base = 1;
                else if (tmp[0] == ng) alt_base = 2;
                else alt_base = 3;
            } else{
                if (tmp[1] == na) alt_base = 0;
                else if (tmp[1] == nc) alt_base = 1;
                else if (tmp[1] == ng) alt_base = 2;
                else alt_base = 3;
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
            if (isinf(fs)) fs = 10000.0;
            else if (fs == 0) fs = 0.0;
            if (alt_fwd * ref_rev > 0) {
                sor = static_cast<double>(ref_fwd * alt_rev) / (ref_rev * alt_fwd);
            } else {
                sor = 10000.0;
            }
            sout << chr << tab << p << tab << BASE2CHAR[ref_base] << tab << aiv.size() << tab << na << tab << nc << tab << ng << tab << nt << tab << fs << tab << sor << tab << ref_fwd << col << ref_rev << col << alt_fwd << col << alt_rev << "\n";
            out = sout.str();
            sout.str("");
            sout.clear();
            *fpc << out;
            // basetype caller;
            BaseType bt(bases, quals, ref_base, min_af);
            if (bt.LRT()) {
                String outvcf = bt.WriteVcf(bt, chr, p, ref_base, aiv, idx, N);
                *fpv << outvcf;
            }
            bases.clear();
            quals.clear();
        }
        aiv.clear();
        idx.clear();
    }
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
    String out = std::to_string(N) + "\t" + std::to_string(M) + "\n";
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

void parseOptions(int argc, char **argv, const char* msg)
{
    bool die = false;
    bool help = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case 'v': opt::verbose = true; break;
        case 'q': arg >> opt::mapq; break;
        case 't': arg >> opt::thread; break;
        case 'l': arg >> opt::input; break;
        case 'r': arg >> opt::reference; break;
        case 'p': arg >> opt::posfile; break;
        case 'g': arg >> opt::region; break;
        case 'o': arg >> opt::output; break;
        default: die = true;
        }
    }
    // todo : need more check
    if (die || help || (opt::input.empty() && opt::output.empty())) {
        std::cerr << msg;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}
