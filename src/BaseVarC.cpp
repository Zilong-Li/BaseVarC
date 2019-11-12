#include <getopt.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>
#include <unistd.h>

#include "htslib/bgzf.h"
#include "RefReader.h"
#include "BamProcess.h"
#include "BaseType.h"
#include "ThreadPool.h"
#define FMT_HEADER_ONLY
#include "fmt/format.h"

#define AUTHOR "Zilong Li"
#define EMAIL "[zimusen94@gmail.com]"

static const char* BASEVARC_USAGE_MESSAGE = 
"Program: BaseVarC -- C++ Version of BaseVar\n"
"Contact: " AUTHOR " " EMAIL "\n"
"Usage  : BaseVarC <command> [options]\n\n"
"Commands:\n"
"           basetype       Variants Caller\n"
"           popmatrix      Create population matrix\n" 
"           concat         Concat popmatrix\n";

static const char* BASETYPE_MESSAGE = 
"Program: BaseVarC basetype\n"
"Contact: " AUTHOR " " EMAIL "\n"
"Usage  : BaseVarC basetype [options]\n\n"
"Commands:\n"
"  --input,      -i        BAM/CRAM file list, one file per row\n"
"  --output,     -o        Output file prefix\n"
"  --reference,  -r        Reference file\n"
"  --region,     -s        Samtools-like region <chr:start-end>\n"
"  --group,      -g        Population group information <SampleID Group>\n"
"  --mapq,       -q <INT>  Mapping quality >= INT [10]\n"
"  --thread,     -t <INT>  Number of threads\n"
"  --batch,      -b <INT>  Number of samples each batch\n"
"  --maf,                  Minimum allele count frequency [min(0.001, 100/N, maf)]\n"
"  --load,                 Load data only\n"
"  --rerun,                Read previous loaded data and rerun\n"
"  --keep_tmp              Don't remove tmp files when basetype finished\n"
"  --verbose,    -v        Set verbose output\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: " AUTHOR " " EMAIL "\n"
"Usage  : BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --input,      -i        BAM/CRAM files list, one file per row.\n"
"  --posfile,    -p        Position file <CHRID POS REF ALT>\n"
"  --output,     -o        Output filename prefix(.mat.gz will be added auto)\n"
"  --mapq,       -q <INT>  Mapping quality >= INT [10]\n"
"  --verbose,    -v        Set verbose output\n";

static const char* CONCAT_MESSAGE =
"Program: BaseVarC concat\n"
"Contact: " AUTHOR " " EMAIL "\n"
"Usage  : BaseVarC concat [options]\n\n"
"Commands:\n"
"  --input,      -i       List of matrix files for concat, one file per row.\n"
"  --output,     -o       Output filename prefix(.gz will be added auto)\n";

static const char* CVG_HEADER =
"##fileformat=CVGv1.0\n"
"##Group information is the depth of A:C:G:T:Indel\n"
"#CHROM\tPOS\tREF\tDepth\tA\tC\tG\tT\tIndels\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)";

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
typedef std::unordered_map<String, int> IndelMap;

struct BtRes
{
    String  cvg;
    String  vcf;
};

void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void runConcat(int argc, char **argv);
void parseOptions(int argc, char **argv, const char* msg);

void bt_r(const StringV& bams, const IntV& pv, const String& refseq, const String& region, const String& fout, int nb, int bc, int ib, int32_t rg_s, int thread);
void bt_s(const StringV& ftmp_v, const IntV& pv, const String& refseq, const String& chr, int32_t rg_s, int32_t N, int thread, int ithread);
BtRes bt_f(int32_t p, const GroupIdx& popg_idx, const AlleleInfoVector& aiv, const DepM& idx, int32_t N, const String& chr, int32_t rg_s, const String& refseq);

namespace opt {
    static bool verbose = false;
    static bool rerun   = false;
    static bool load    = false;
    static bool keep_tmp= false;
    static uint8_t mapq = 10;
    static int thread = 1;
    static int batch  = 10;
    static double maf = 0.001;
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
  { "keep_tmp",                no_argument, NULL,  6  },
  { "load",                    no_argument, NULL,  7  },
  { "rerun",                   no_argument, NULL,  8  },
  { "maf",                     required_argument, NULL,  9  },
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
    time_t tim = time(0);
    clock_t ctb = clock();
    std::cout << "basetype start -- " << ctime(&tim);
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
    String refseq = fa.GetTargetBase(opt::region, opt::reference);
    String chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(opt::region);
    IntV pv;
    String acgt = "ACGT";
    for (size_t i = 0; i < refseq.length(); ++i) {
        // skip non-acgt character
        if (acgt.find(refseq[i]) != std::string::npos) {
            pv.push_back(i + rg_s);       // 1-based
        }
    }
    // begin to read bams
    String tmp;
    int thread = opt::thread;
    for (int i = 0; i < thread; ++i) {
        tmp = fmt::format("mkdir -p {}.tmp.thread.{}", opt::output, i);
        if (!std::system(tmp.c_str())) continue;
        else { std::cerr << "Error: couldn't mkdir " << std::endl; exit(EXIT_FAILURE);}
    }
    std::vector<StringV> ftmp_vv(thread);
    int bc = opt::batch;
    int ngz = 0, nb = 1 + (N - 1) / bc;    // ceiling
    int bk = nb - 1;
    for (int j = 0; j < nb; ++j) {
        int k = 0;
        for (int i = 0; i < thread; ++i) {
            tmp = fmt::format("{}.tmp.thread.{}/batch.{}", opt::output, i, j);
            ftmp_vv[i].push_back(tmp);
            if (bgzf_is_bgzf(tmp.c_str())) {
                BGZF* fp = bgzf_open(tmp.c_str(), "r");
                if (bgzf_check_EOF(fp) == 1) { k +=1; ngz += 1; }
            }
        }
        if (k != thread) bk = bk > j ? j : bk;
    }
    if (opt::rerun && ngz > 0) {
        if (ngz != thread * nb) {
            BaseVar::ThreadPool pool(thread);
            std::vector<std::future<void>> res;
            std::cerr << "begin to read bams and save as tmp file" << std::endl;
            for (int i = bk; i < nb; ++i) {
                res.emplace_back(pool.enqueue(bt_r, std::cref(bams), std::cref(pv), std::cref(refseq), std::cref(opt::region), std::cref(opt::output), nb, bc, i, rg_s, thread));
            }
            for (auto && r: res) {
                r.get();
            }
            res.clear();
        }
    } else {
        BaseVar::ThreadPool pool(thread);
        std::vector<std::future<void>> res;
        std::cerr << "begin to read bams and save as tmp file" << std::endl;
        for (int i = 0; i < nb; ++i) {
            res.emplace_back(pool.enqueue(bt_r, std::cref(bams), std::cref(pv), std::cref(refseq), std::cref(opt::region), std::cref(opt::output), nb, bc, i, rg_s, thread));
        }
        for (auto && r: res) {
            r.get();
        }
        res.clear();
    }
    time_t tim1 = time(0);
    std::cout << "basetype loading done -- " << ctime(&tim1);
    if (opt::load) exit(EXIT_SUCCESS);
    // begin to call basetype
    std::vector<std::thread> workers;
    for (int i = 0; i < thread; ++i) {
        workers.push_back(std::thread(bt_s, std::cref(ftmp_vv[i]), std::cref(pv), std::cref(refseq), std::cref(chr), rg_s, N, thread, i));
    }
    // merge all subfile
    String vcfout = opt::output + ".vcf.gz", subvcf;
    String cvgout = opt::output + ".cvg.gz", subcvg;
    BGZF* fov = bgzf_open(vcfout.c_str(), "w");
    BGZF* foc = bgzf_open(cvgout.c_str(), "w");
    BGZF* fiv = NULL; BGZF* fic = NULL;
    kstring_t ks = {0, 0, NULL};
    for (int i = 0; i < thread; ++i) {
        auto & t = workers[i];
        if (t.joinable()) t.join();
        subvcf = fmt::format("{}.{}.vcf.gz", opt::output, i);
        subcvg = fmt::format("{}.{}.cvg.gz", opt::output, i);
        fiv = bgzf_open(subvcf.c_str(), "r");
        fic = bgzf_open(subcvg.c_str(), "r");
        while (bgzf_getline(fiv, '\n', &ks) >= 0) {
            tmp = (String)ks.s + '\n';
            if (bgzf_write(fov, tmp.c_str(), tmp.length()) != tmp.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        while (bgzf_getline(fic, '\n', &ks) >= 0) {
            tmp = (String)ks.s + '\n';
            if (bgzf_write(foc, tmp.c_str(), tmp.length()) != tmp.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::remove(subvcf.c_str());
        std::remove(subcvg.c_str());
    }
    std::cout << "merge subfiles done" << std::endl;
    if (bgzf_close(fov) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    if (bgzf_close(foc) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    for (int i = 0; i < thread; ++i) {
        tmp = fmt::format("{}.tmp.thread.{}", opt::output, i);
        // for unix-system;
        rmdir(tmp.c_str());
    }

    // done
    time_t tim2 = time(0);
    std::cout << "basetype computing done -- " << ctime(&tim2);
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "basetype elapsed cpu secs : " << elapsed_secs << std::endl;
    std::cout << "basetype done" << std::endl;
}

void bt_s(const StringV& ftmp_v, const IntV& pv, const String& refseq, const String& chr, int32_t rg_s, int32_t N, int thread, int ithread)
{
    // hold all tmp file pointers
    String headcvg = String(CVG_HEADER);
    String headvcf = String(VCF_HEADER);
    String vcfout = fmt::format("{}.{}.vcf.gz", opt::output, ithread);
    String cvgout = fmt::format("{}.{}.cvg.gz", opt::output, ithread);
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
    // get contig from fai file.
    String fai = opt::reference + ".fai";
    std::ifstream ifai(fai);
    if (ifai.is_open()) {
        String contig, len, t;
        while (ifai >> contig >> len >> t >> t >> t) {
            headvcf += "##contig=<ID=" + contig + ",length=" + len + ">\n";
        }
    }
    headvcf += "##reference=file://" + opt::reference + "\n";
    headvcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sams + "\n";
    headcvg += "\n";
    // output header
    if (ithread == 0) {
        if (bgzf_write(fpc, headcvg.c_str(), headcvg.length()) != headcvg.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (bgzf_write(fpv, headvcf.c_str(), headvcf.length()) != headvcf.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    // begin to call basetype and output
    AlleleInfo ai;
    AlleleInfoVector aiv;
    DepM idx;
    int32_t j, k, i, count=0;
    std::cerr << "begin to load data and run basetype" << std::endl;
    char *buf=NULL, *str=NULL, *str2=NULL, *pti=NULL, *pto=NULL;
    int32_t window = pv.size() % thread + pv.size() / thread;
    IntV::const_iterator itp, itp2;
    if (ithread == thread - 1) itp2 = pv.end();
    else itp2 = pv.begin() + (ithread + 1) * window;
    for (itp = pv.begin() + ithread * window; itp != itp2; ++itp) {
        auto & p = *itp;
        j = 0; k = 0;
        // merge all data together from tmp files
        for (auto & fp: fpiv) {
            if (bgzf_getline(fp, '\n', &ks) >= 0) {
                buf = ks.s;
                while ((str = strtok_r(buf, " ", &pto)) != NULL) {
                    if (str[0] != '+' && str[0] != '-' && str[0] != 'N' && str[0] != '.') {
                        buf = str;
                        ai.is_indel = 0;
                        for (i = 0; i < 5; ++i) {
                            if ((str2 = strtok_r(buf, ",", &pti)) != NULL) {
                                switch(i){
                                case 0: ai.base = std::atoi(str2);break;
                                case 1: ai.mapq = std::atoi(str2);break;
                                case 2: ai.qual = std::atoi(str2);break;
                                case 3: ai.rpr  = std::atoi(str2);break;
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
                    } else if (str[0] != '.') {
                        ai.is_indel = 1;
                        ai.indel = str;
                        aiv.push_back(ai);
                        idx.insert({j, k++});
                    }
                    buf = NULL;
                    j++;
                }
            }
        }
        if (!aiv.empty()) {
            auto btr = bt_f(p, popg_idx, aiv, idx, N, chr, rg_s, refseq);
            if (!btr.vcf.empty() && bgzf_write(fpv, btr.vcf.c_str(), btr.vcf.length()) != btr.vcf.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (bgzf_write(fpc, btr.cvg.c_str(), btr.cvg.length()) != btr.cvg.length()) {
                std::cerr << "fail to write - exit" << std::endl;
                exit(EXIT_FAILURE);
            }
            aiv.clear();
            idx.clear();
            if (!(++count % 1000)) std::cerr << "basetype completed " << count << " sites -- thread" << ithread << std::endl;
        }
    }
    if (bgzf_close(fpv) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    if (bgzf_close(fpc) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    // whether remove tmp file or not
    if (!opt::keep_tmp) {
        for (auto & f: ftmp_v) {
            std::remove(f.c_str());
        }
    }
}

void bt_r(const StringV& bams, const IntV& pv, const String& refseq, const String& region, const String& fout, int nb, int bc, int ib, int32_t rg_s, int thread)
{
    PosAlleleMapVec allele_mv;
    String names, fw;
    StringV::const_iterator itb = bams.begin() + ib * bc, itb2;
    if (ib == nb - 1) {
        itb2 = bams.end();
    } else {
        itb2 = bams.begin() + (ib + 1) * bc;
    }
    int32_t size = itb2 - itb, count = 0;
    allele_mv.reserve(size);
    for (; itb != itb2; ++itb) {
        auto bam = *itb;
        try {
        BamProcess reader(opt::mapq);
        if (!(++count % 100)) std::cerr << "reading number " << count << " bam -- " << fout << ".tmp.batch." << ib << std::endl;
        if (!reader.Open(bam)) {
            std::cerr << "ERROR: could not open file " << bam << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.FindSnpAtPos(rg_s, refseq, region, pv)) {
            std::cerr << "Warning: " << reader.sm << " region " << region << " is empty." << std::endl;
        }
        allele_mv.push_back(reader.allele_m);
        names += reader.sm + '\t';
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bam << std::endl;
        }
        } catch (std::out_of_range e) {
            std::cout <<  bam << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    names += "\n";    // we keep '\t' ahead of '\n' in order to connect different batches' names directly
    int32_t psize = pv.size();
    int32_t window = psize % thread + psize / thread;
    std::vector<BGZF*> fpv;
    BGZF* fp;
    for (int i = 0; i < thread; ++i) {
        fw = fmt::format("{}.tmp.thread.{}/batch.{}", fout, i, ib);
        fp = bgzf_open(fw.c_str(), "w");
        if (bgzf_write(fp, names.c_str(), names.length()) != names.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        fpv.push_back(fp);
    }
    String out;
    for (int i = 0, j = 0; i < psize; ++i) {
        auto p = pv[i];
        out = "";
        for (auto const& m : allele_mv) {
            if (m.count(p)) {
                auto const& a = m.at(p);
                if (a.is_indel == 1) out += fmt::format("{} ", a.indel);
                else out += fmt::format("{},{},{},{},{} ", a.base, a.mapq, a.qual, a.rpr, a.strand);
            } else {
                out += ". ";
            }
        }
        out += "\n";
        if (i == (j + 1) * window) ++j;
        if (bgzf_write(fpv[j], out.c_str(), out.length()) != out.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    for (auto & fp : fpv) {
        if (bgzf_close(fp) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    }
}

BtRes bt_f(int32_t p, const GroupIdx& popg_idx, const AlleleInfoVector& aiv, const DepM& idx, int32_t N, const String& chr, int32_t rg_s, const String& refseq)
{
    int8_t alt_base, ref_base;
    int32_t dep, na, nc, ng, nt, ref_fwd, ref_rev, alt_fwd, alt_rev;
    double fs, sor;
    double min_af = 100.0 / N;
    if (min_af > 0.001) min_af = 0.001;
    if (opt::maf < min_af ) min_af = opt::maf;
    BaseV bases, quals;
    BtRes res;
    IndelMap indel_m;
    // output cvg;
    ref_base = BASE_INT8_TABLE[static_cast<int>(refseq[p - rg_s])];
    na = 0; nc = 0; ng = 0; nt = 0;
    for (auto const& a: aiv) {
        if (a.is_indel == 0) {
            switch (a.base) {
            case 0 : na += 1; break;
            case 1 : nc += 1; break;
            case 2 : ng += 1; break;
            case 3 : nt += 1; break;
            }
            bases.push_back(a.base);
            quals.push_back(a.qual);
        } else if (a.is_indel == 1) {
            if (indel_m.count(a.indel)) {
                indel_m[a.indel] += 1;
            } else {
                indel_m.insert({a.indel, 1});
            }
        }
    }
    String oss, indels = ".";
    if (!indel_m.empty()) {
        indels = "";
        for (IndelMap::iterator it = indel_m.begin(); it != indel_m.end(); ++it) {
            indels += fmt::format("{}|{},", it->first, it->second);
        }
        indels.pop_back();
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
    fs = bt_fisher_exact(ref_fwd, ref_rev, alt_fwd, alt_rev);
    if (alt_fwd * ref_rev > 0) {
        sor = static_cast<double>(ref_fwd * alt_rev) / (ref_rev * alt_fwd);
    } else {
        sor = 10000.0;
    }
    dep = na + nc + ng + nt;
    oss = fmt::format("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}\t{},{},{},{}\t", chr, p, BASE2CHAR[ref_base], dep, na, nc, ng, nt, indels, fs, sor, ref_fwd, ref_rev, alt_fwd, alt_rev);
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
                    if (a.is_indel == 0) {
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
            }
            oss += fmt::format("{}:{}:{}:{}\t", na, nc, ng, nt);
            if (!gr_bases.empty()) {
                BaseType gr_bt(gr_bases, gr_quals, ref_base, min_af);
                gr_bt.SetBase(base_comb);
                gr_bt.LRT();
                gr_af = "";
                for (auto b : bt.alt_bases) {
                    if (gr_bt.af_lrt.count(b)) {
                        gr_af += fmt::format("{:.6f},", gr_bt.af_lrt[b]);
                    } else {
                        gr_af += "0,";
                    }
                }
                gr_af.pop_back();
                gr_bases.clear();
                gr_quals.clear();
                info.insert({it->first + "_AF", gr_af});
            } else {
                info.insert({it->first + "_AF", "0"});
            }
        }
    }
    oss.pop_back(); oss += "\n";
    res.cvg = oss;
    if (bt_success) {
        res.vcf = WriteVcf(bt, chr, p, ref_base, aiv, idx, info, N);
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
    String out = fmt::format("{}\t{}\n", N, M);
    BGZF* fp = bgzf_open(opt::output.c_str(), "w");
    if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // ready for run
    String rg = pv.front() + pv.back();
    int32_t count = 0;
    for (int32_t i = 0; i < N; i++) {
        BamProcess reader(opt::mapq);
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
    ss = fmt::format("{}\t{}\n", n, mt);
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
        case  9 : arg >> opt::maf; break;
        case  8 : opt::rerun   = true; break;
        case  7 : opt::load    = true; break;
        case  6 : opt::keep_tmp= true; break;
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
