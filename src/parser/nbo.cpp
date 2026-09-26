#include "cov/nbo.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cov { namespace {
const std::string num=R"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?)";
std::string trim(std::string s) { auto a=s.find_first_not_of(" \t\r\n"); return a==s.npos?"":s.substr(a,s.find_last_not_of(" \t\r\n")-a+1); }
std::string upper(std::string s) { for(auto& c:s)c=static_cast<char>(std::toupper(static_cast<unsigned char>(c))); return s; }
double number(std::string s) { for(auto& c:s)if(c=='D'||c=='d')c='E'; double x=std::stod(s); if(!std::isfinite(x))throw std::runtime_error("NBO nonfinite number"); return x; }
std::vector<double> numbers(const std::string& s) { static const std::regex r(num); std::vector<double> v; for(std::sregex_iterator i(s.begin(),s.end(),r),e;i!=e;++i)v.push_back(number(i->str())); return v; }
std::vector<std::string> lines(const std::filesystem::path& p) { std::ifstream f(p); if(!f)throw std::runtime_error("Cannot read NBO file: "+p.string()); std::vector<std::string> l; std::string s; while(std::getline(f,s))l.push_back(s); return l; }
NboSource source(const std::filesystem::path& p,std::string block,std::size_t line,std::string raw={},std::string version={},std::size_t seg=0) { return {p.string(),std::move(block),std::move(version),std::move(raw),line,line,seg}; }
const NboMatrix* matrix(const std::vector<NboMatrix>& a,const std::string& kind,NboSpin spin) { for(auto& m:a)if(m.kind==kind&&m.spin==spin)return &m; return nullptr; }
std::vector<double> multiply(const std::vector<double>& a,const std::vector<double>& b,std::size_t n) { std::vector<double> c(n*n); for(std::size_t i=0;i<n;++i)for(std::size_t k=0;k<n;++k)for(std::size_t j=0;j<n;++j)c[i*n+j]+=a[i*n+k]*b[k*n+j]; return c; }
std::string quote(const std::string& s) { std::ostringstream o;o<<'"'; for(unsigned char c:s) { if(c=='"'||c=='\\')o<<'\\'<<c;else if(c=='\n')o<<"\\n";else if(c=='\r')o<<"\\r";else if(c=='\t')o<<"\\t";else if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else o<<c; }o<<'"';return o.str(); }
void json_source(std::ostream& o,const NboSource& s) { o<<"{\"path\":"<<quote(s.path)<<",\"block\":"<<quote(s.block)<<",\"producer_version\":"<<quote(s.producer_version)<<",\"line_begin\":"<<s.line_begin<<",\"line_end\":"<<s.line_end<<",\"analysis_segment\":"<<s.analysis_segment<<",\"raw\":"<<quote(s.raw)<<"}"; }
void optional_number(std::ostream& o,const std::optional<double>& v) { if(v)o<<*v;else o<<"null"; }
template<class T,class F> void array(std::ostream& o,const std::vector<T>& a,F f) { o<<'[';bool first=true;for(auto& x:a){if(!first)o<<',';first=false;f(x);}o<<']'; }
}
const char* nbo_spin_name(NboSpin s) noexcept { return s==NboSpin::Alpha?"alpha":s==NboSpin::Beta?"beta":"total"; }

NboArchive read_nbo_archive(const std::filesystem::path& path) {
    auto l=lines(path); std::map<std::string,std::string> blocks; std::map<std::string,std::size_t> starts,ends;
    std::string key;
    for(std::size_t i=0;i<l.size();++i) { std::string s=l[i]; std::size_t pos=0;
        while(pos<s.size()) { auto dollar=s.find('$',pos); if(dollar==s.npos){if(!key.empty())blocks[key]+=s.substr(pos)+"\n";break;}
            if(!key.empty())blocks[key]+=s.substr(pos,dollar-pos)+"\n";
            auto end=dollar+1;while(end<s.size()&&std::isalnum(static_cast<unsigned char>(s[end])))++end;
            auto name=upper(s.substr(dollar+1,end-dollar-1));if(name=="END"){if(!key.empty())ends[key]=i+1;key.clear();}else{if(!key.empty())throw std::runtime_error("Nested .47 blocks are unsupported");if(blocks.count(name))throw std::runtime_error("Duplicate .47 block: "+name);key=name;blocks[key]="";starts[key]=i+1;} pos=end;
        }
    }
    if(!key.empty())throw std::runtime_error("Unterminated .47 block");
    auto required=[&](const std::string& k)->const std::string& {auto i=blocks.find(k);if(i==blocks.end())throw std::runtime_error("Missing .47 $"+k);return i->second;};
    const auto gen=upper(required("GENNBO")); std::smatch m; NboArchive a;a.source=source(path,"GENNBO",starts["GENNBO"]);a.source.line_end=l.size();
    if(!std::regex_search(gen,m,std::regex(R"(NBAS\s*=\s*(\d+))")))throw std::runtime_error("Missing NBAS");a.basis_count=std::stoul(m[1]);
    if(a.basis_count==0||a.basis_count>10000)throw std::runtime_error("Unsupported .47 NBAS");
    if(!std::regex_search(gen,m,std::regex(R"(NATOMS\s*=\s*(\d+))")))throw std::runtime_error("Missing NATOMS");auto nat=std::stoul(m[1]);
    a.open_shell=std::regex_search(gen,std::regex(R"(\bOPEN\b)"));a.density_is_bond_order=std::regex_search(gen,std::regex(R"(\bBODM\b)"));bool packed=std::regex_search(gen,std::regex(R"(\bUPPER\b)"));
    bool bohr=std::regex_search(gen,std::regex(R"(\bBOHR\b)"));
    std::istringstream coords(required("COORD"));std::string row;
    const std::regex cr("^\\s*(\\d+)\\s+("+num+")\\s+("+num+")\\s+("+num+")\\s+("+num+")\\s*$");
    while(std::getline(coords,row))if(std::regex_match(row,m,cr)){Atom at;at.atomic_number=std::stoi(m[1]);at.nuclear_charge=number(m[2]);const double u=bohr?1:kAngstromToBohr;at.x=number(m[3])*u;at.y=number(m[4])*u;at.z=number(m[5])*u;a.atoms.push_back(at);}
    if(a.atoms.size()!=nat)throw std::runtime_error(".47 COORD atom count mismatch");
    auto fields=[](const std::string& block){std::map<std::string,std::vector<double>> f;const std::regex r(R"(([A-Za-z][A-Za-z0-9]*)\s*=)");std::vector<std::pair<std::string,std::size_t>> keys;std::vector<std::size_t> ends;for(std::sregex_iterator i(block.begin(),block.end(),r),e;i!=e;++i){keys.push_back({upper((*i)[1]),static_cast<std::size_t>(i->position()+i->length())});ends.push_back(i->position());}for(std::size_t i=0;i<keys.size();++i)f[keys[i].first]=numbers(block.substr(keys[i].second,(i+1<keys.size()?ends[i+1]:block.size())-keys[i].second));return f;};
    auto b=fields(required("BASIS"));auto c=fields(required("CONTRACT"));
    auto ints=[](const std::vector<double>& x){std::vector<int> y;for(double v:x){if(v!=std::floor(v))throw std::runtime_error("Noninteger .47 metadata");y.push_back(static_cast<int>(v));}return y;};
    a.centers=ints(b["CENTER"]);a.labels=ints(b["LABEL"]);a.ncomp=ints(c["NCOMP"]);a.nprim=ints(c["NPRIM"]);a.nptr=ints(c["NPTR"]);a.exponents=c["EXP"];a.cs=c["CS"];a.cp=c["CP"];a.cd=c["CD"];a.cf=c["CF"];a.cg=c["CG"];
    if(a.centers.size()!=a.basis_count||a.labels.size()!=a.basis_count||a.ncomp.size()!=a.nprim.size()||a.nprim.size()!=a.nptr.size()||std::accumulate(a.ncomp.begin(),a.ncomp.end(),0)!=static_cast<int>(a.basis_count))throw std::runtime_error(".47 basis/contract dimensions inconsistent");
    for(std::size_t i=0;i<a.nprim.size();++i)if(a.nprim[i]<=0||a.nptr[i]<=0||static_cast<std::size_t>(a.nptr[i]-1+a.nprim[i])>a.exponents.size())throw std::runtime_error(".47 primitive pointer outside exponent array");
    for(auto& k:{"OVERLAP","DENSITY","FOCK","LCAOMO"}) { auto it=blocks.find(k);if(it==blocks.end())continue; auto v=numbers(it->second);const std::size_t n=a.basis_count, count=(packed&&std::string(k)!="LCAOMO")?n*(n+1)/2:n*n;const std::size_t spins=a.open_shell&&std::string(k)!="OVERLAP"?2:1;
        if(v.size()!=count*spins)throw std::runtime_error(std::string(".47 matrix size mismatch: ")+k);
        for(std::size_t s=0;s<spins;++s){NboMatrix mat;mat.kind=k;mat.spin=spins==1?NboSpin::Total:s==0?NboSpin::Alpha:NboSpin::Beta;mat.rows=mat.columns=n;mat.values.resize(n*n);mat.source=source(path,k,starts[k]);mat.source.line_end=ends[k];std::size_t p=s*count;for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<(count==n*n?n:j+1);++i){mat.values[i*n+j]=v[p++];if(count!=n*n)mat.values[j*n+i]=mat.values[i*n+j];}a.matrices.push_back(std::move(mat));}
    }
    return a;
}

namespace {
std::vector<NboMatrix> read_matrix_file(const std::filesystem::path& p,const std::string& kind,std::size_t n,bool open) {
    auto ls=lines(p);
    // W37/PLOT appends occupancies and orbital/atom descriptors to each spin
    // block. Only the first NBAS*NBAS real values belong to the matrix.
    // OPEN files delimit their independent payloads with ALPHA/BETA SPIN.
    const std::regex numeric("^\\s*"+num+"(?:\\s+"+num+")*\\s*$");
    const std::string heading=kind=="AONBO"?"NBOs in the AO basis:":"MOs in the NBO basis:";
    std::size_t header=ls.size();for(std::size_t i=0;i<ls.size();++i)if(ls[i].find(heading)!=ls[i].npos){if(header!=ls.size())throw std::runtime_error("Repeated matrix header");header=i;}
    if(header==ls.size())throw std::runtime_error("NBO matrix type header absent or mismatched: "+kind);
    std::vector<std::size_t> starts;if(open){for(std::size_t i=header+1;i<ls.size();++i){auto u=upper(trim(ls[i]));if(u=="ALPHA SPIN"||u=="BETA  SPIN"||u=="BETA SPIN")starts.push_back(i+1);}if(starts.size()!=2||upper(ls[starts[0]-1]).find("ALPHA")==std::string::npos||upper(ls[starts[1]-1]).find("BETA")==std::string::npos)throw std::runtime_error("OPEN matrix requires exactly one alpha and beta block");}else starts={header+1};
    const std::size_t count=n*n;std::vector<NboMatrix> out;
    for(std::size_t s=0;s<starts.size();++s){std::vector<double> values;std::size_t first=0,last=0;const auto end=s+1<starts.size()?starts[s+1]-1:ls.size();for(std::size_t i=starts[s];i<end&&values.size()<count;++i){if(!std::regex_match(ls[i],numeric)){if(!values.empty())throw std::runtime_error("Interrupted NBO matrix payload");continue;}if(ls[i].find_first_of(".EeDd")==std::string::npos){if(values.empty())continue;throw std::runtime_error("Integer metadata encountered before complete NBO matrix");}if(!first)first=i+1;auto row=numbers(ls[i]);if(values.size()+row.size()>count)throw std::runtime_error("NBO matrix payload boundary not aligned");values.insert(values.end(),row.begin(),row.end());last=i+1;}if(values.size()!=count)throw std::runtime_error(kind+" incomplete matrix payload");NboMatrix m;m.kind=kind;m.spin=open?(s==0?NboSpin::Alpha:NboSpin::Beta):NboSpin::Total;m.rows=m.columns=n;m.values.resize(count);m.source=source(p,kind,first);m.source.line_end=last;for(std::size_t j=0;j<n;++j)for(std::size_t i=0;i<n;++i)m.values[i*n+j]=values[j*n+i];out.push_back(std::move(m));}return out;
}
}

NboDataset read_nbo(const std::filesystem::path& path,const NboReadOptions& options) {
    auto ls=lines(path);NboDataset d;std::vector<std::size_t> segments;
    const std::regex banner(R"(\*+\s+NBO\s+[0-9])");for(std::size_t i=0;i<ls.size();++i)if(std::regex_search(ls[i],banner))segments.push_back(i);
    if(segments.empty())segments.push_back(0);if(segments.size()>1&&!options.analysis_segment)throw std::runtime_error("Multiple NBO analyses: explicit analysis_segment is required");
    auto seg=options.analysis_segment.value_or(0);if(seg>=segments.size())throw std::runtime_error("NBO analysis_segment outside available segments");const auto begin=segments[seg],end=seg+1<segments.size()?segments[seg+1]:ls.size();
    d.source=source(path,"NBO",begin+1,{}, {},seg);d.source.line_end=end;
    std::smatch m;for(std::size_t i=begin;i<end;++i)if(std::regex_search(ls[i],m,std::regex(R"(\[NBO\s+([^\]]+)\])"))){d.producer_version=m[1];break;}if(d.producer_version.empty())for(std::size_t i=begin;i<end;++i)if(std::regex_search(ls[i],m,std::regex(R"(NBO\s+([0-9]+\.[0-9]+))"))){d.producer_version=m[1];break;}d.source.producer_version=d.producer_version;
    auto src=[&](const char* block,std::size_t i){return source(path,block,i+1,ls[i],d.producer_version,seg);};
    NboSpin spin=NboSpin::Total;std::string section;std::size_t current=static_cast<std::size_t>(-1);std::vector<std::size_t> columns;bool nao_spin_column=false;
    const std::regex nao("^\\s*(\\d+)\\s+([A-Za-z]+)\\s+(\\d+)\\s+(\\S+)\\s+([A-Za-z]+\\([^)]*\\))\\s+("+num+")(?:\\s+("+num+"))?\\s*$");
    const std::regex pop("^\\s*([A-Za-z]+)\\s+(\\d+)\\s+("+num+")\\s+("+num+")\\s+("+num+")\\s+("+num+")\\s+("+num+")(?:\\s+("+num+"))?\\s*$");
    const std::regex orb("^\\s*(\\d+)\\.\\s*\\(("+num+")\\)\\s*([A-Za-z0-9]+\\*?)\\s*\\(\\s*(\\d+)\\)\\s*(.*)$");
    const std::regex atom(R"(([A-Z][a-z]?)\s+(\d+))");
    const std::regex comp("^\\s*\\(\\s*("+num+")%\\)\\s*("+num+")\\*\\s*([A-Z][a-z]?)\\s+(\\d+)\\s*(.*)$");
    const std::regex e2("^\\s*(\\d+)\\.\\s+.*?\\s+(\\d+)\\.\\s+.*?\\s+("+num+")\\s+("+num+")\\s+("+num+")\\s*$");
    // Printed occupations have a decimal point; this prevents the final atom
    // index from being mistaken for occupation when summary annotations follow.
    const std::regex summ("^\\s*(\\d+)\\.\\s+([A-Za-z0-9]+\\*?)\\s*\\(\\s*\\d+\\).*?\\s+([0-9]+\\.[0-9]+)\\s+("+num+")(?:\\s+.*)?$");
    for(std::size_t i=begin;i<end;++i){const auto& l=ls[i];auto u=upper(l);
        if(u.find("ALPHA SPIN ORBITALS")!=u.npos){spin=NboSpin::Alpha;section.clear();current=-1;}if(u.find("BETA  SPIN ORBITALS")!=u.npos||u.find("BETA SPIN ORBITALS")!=u.npos){spin=NboSpin::Beta;section.clear();current=-1;}
        if(l.find("NATURAL POPULATIONS:")!=l.npos){section="nao";nao_spin_column=false;}
        if(section=="nao"&&l.find("NAO Atom No")!=l.npos)nao_spin_column=l.find("Spin")!=l.npos;
        if(l.find("Summary of Natural Population Analysis:")!=l.npos)section="npa";
        if(l.find("(Occupancy)")!=l.npos&&l.find("Bond orbital")!=l.npos)section="nbo";
        if(l.find("SECOND ORDER PERTURBATION")!=l.npos){section="e2";NboE2Section s;s.spin=spin;s.source=src("E2",i);s.missing_reason="Unprinted pairs are below threshold or omitted; never interpreted as zero";d.e2_sections.push_back(s);}
        if(l.find("NATURAL BOND ORBITALS (Summary)")!=l.npos)section="summary";
        if(section=="summary"&&l.find("Total Lewis")!=l.npos)section.clear();
        if(l.find("Wiberg bond index matrix")!=l.npos){section="wiberg";columns.clear();}
        if(section=="wiberg"&&l.find("Wiberg bond index, Totals")!=l.npos)section.clear();
        if(l.find("CMO: NBO Analysis")!=l.npos){section="cmo";d.cmo_summaries.push_back(src("CMO thresholded summary",i));}
        if(section=="cmo"){if(!d.cmo_summaries.empty()&&d.cmo_summaries.back().line_begin!=i+1){d.cmo_summaries.back().raw+="\n"+l;d.cmo_summaries.back().line_end=i+1;}continue;}
        if(section=="nao"&&std::regex_match(l,m,nao)){NboNao x;x.id=std::stoul(m[1]);x.symbol=m[2];x.atom=std::stoul(m[3]);x.angular=m[4];x.type=m[5];x.occupation=number(m[6]);if(m[7].matched){if(nao_spin_column)x.spin_density=number(m[7]);else x.energy_hartree=number(m[7]);}x.spin=spin;x.source=src(nao_spin_column?"NAO occupation and spin density":"NAO occupation and energy",i);d.naos.push_back(x);}
        if(section=="npa"&&std::regex_match(l,m,pop)){NboPopulation x;x.symbol=m[1];x.atom=std::stoul(m[2]);x.charge=number(m[3]);x.core=number(m[4]);x.valence=number(m[5]);x.rydberg=number(m[6]);x.total=number(m[7]);if(m[8].matched)x.spin_density=number(m[8]);x.spin=spin;x.source=src("NPA",i);d.populations.push_back(x);}
        if(section=="nbo"&&std::regex_match(l,m,orb)){NboOrbital x;x.id=std::stoul(m[1]);x.occupation=number(m[2]);x.kind=m[3];x.ordinal=std::stoul(m[4]);x.label=x.kind+"("+m[4].str()+") "+trim(m[5]);x.spin=spin;x.source=src("NBO",i);auto tail=m[5].str();for(std::sregex_iterator ai(tail.begin(),tail.end(),atom),ae;ai!=ae;++ai)x.atoms.push_back(std::stoul((*ai)[2]));if(x.atoms.size()==1){NboLocalComponent c;c.atom=x.atoms[0];c.percent=100;c.coefficient=1;c.hybrid=trim(std::regex_replace(tail,atom,"",std::regex_constants::format_first_only));c.source=x.source;x.components.push_back(c);}d.orbitals.push_back(x);current=d.orbitals.size()-1;}
        else if(section=="nbo"&&current<d.orbitals.size()&&std::regex_match(l,m,comp)){NboLocalComponent c;c.percent=number(m[1]);c.coefficient=number(m[2]);c.atom=std::stoul(m[4]);c.hybrid=trim(m[5]);c.source=src("NBO hybrid",i);d.orbitals[current].components.push_back(c);}
        if(section=="summary"&&std::regex_match(l,m,summ)){auto id=std::stoul(m[1]);for(auto& o:d.orbitals)if(o.id==id&&o.spin==spin){if(o.energy_source)throw std::runtime_error("Duplicate NBO energy identity inside summary");if(o.kind!=m[2].str()||std::abs(o.occupation-number(m[3]))>1e-5)throw std::runtime_error("NBO summary identity/occupation mismatch");o.diagonal_fock_hartree=number(m[4]);o.energy_source=src("NBO diagonal Fock summary",i);break;}}
        if(section=="e2") { if(std::regex_search(l,m,std::regex("Threshold for printing:\\s*("+num+")\\s*(\\S+)"))){d.e2_sections.back().printing_threshold=number(m[1]);d.e2_sections.back().units=m[2];}if(std::regex_match(l,m,e2)){NboE2 x;x.donor=std::stoul(m[1]);x.acceptor=std::stoul(m[2]);x.value=number(m[3]);x.energy_gap_hartree=number(m[4]);x.fock_hartree=number(m[5]);x.spin=spin;x.printing_threshold=d.e2_sections.back().printing_threshold;x.units=d.e2_sections.back().units;x.source=src("E2",i);d.e2.push_back(x);}d.e2_sections.back().source.line_end=i+1; }
        if(section=="wiberg") {if(std::regex_search(l,std::regex(R"(^\s*Atom\s+\d)"))){columns.clear();for(double v:numbers(l))columns.push_back(static_cast<std::size_t>(v));}else if(std::regex_match(l,m,std::regex(R"(^\s*(\d+)\.\s+([A-Za-z]+)\s+(.*)$)"))){auto v=numbers(m[3]);if(v.size()!=columns.size())throw std::runtime_error("Incomplete Wiberg row");for(std::size_t j=0;j<v.size();++j){NboWiberg x;x.atom_a=std::stoul(m[1]);x.atom_b=columns[j];x.value=v[j];x.spin=spin;x.source=src("Wiberg NAO",i);d.wiberg.push_back(x);}}}
    }
    if(d.e2_sections.empty()){NboE2Section s;s.source=d.source;s.missing_reason="E2 analysis not printed or not requested; unavailable, not zero";d.e2_sections.push_back(s);}
    if(!options.archive47.empty())d.archive=read_nbo_archive(options.archive47);
    if(d.archive)for(auto& p:d.populations){if(p.atom==0||p.atom>d.archive->atoms.size())throw std::runtime_error("NPA atom outside archive");const auto& a=d.archive->atoms[p.atom-1];p.effective_core_electrons=(a.atomic_number-a.nuclear_charge)*(p.spin==NboSpin::Total?1.0:0.5);p.explicit_population=p.total-*p.effective_core_electrons;}
    std::size_t n=d.archive?d.archive->basis_count:0;bool open=d.archive?d.archive->open_shell:false;
    if(!options.aonbo.empty()||!options.nbomo.empty()){if(!n)throw std::runtime_error("Complete matrix import requires .47 dimensions/spin metadata");for(auto& item:std::vector<std::pair<std::filesystem::path,std::string>>{{options.aonbo,"AONBO"},{options.nbomo,"NBOMO"}})if(!item.first.empty()){auto v=read_matrix_file(item.first,item.second,n,open);d.matrices.insert(d.matrices.end(),v.begin(),v.end());}}
    std::set<std::pair<NboSpin,std::size_t>> ids;for(auto& o:d.orbitals)if(!ids.insert({o.spin,o.id}).second)throw std::runtime_error("Ambiguous repeated NBO orbital identity in selected analysis");
    if(!d.cmo_summaries.empty())d.warnings.push_back("CMO printed percentages are thresholded summaries; complete decomposition requires NBOMO matrix");
    return d;
}

NboAssociation associate_nbo(NboDataset& d,const Wavefunction& w) {
    NboAssociation r;
    auto fail=[&](const std::string& status,const std::string& why){r.status=status;r.detail=why;d.association=r;return r;};
    if(!d.archive)return fail("missing_archive","Same-source .47 archive is required");const auto& a=*d.archive;const auto n=a.basis_count;
    if(!n||a.centers.size()!=n||a.labels.size()!=n||a.ncomp.size()!=a.nprim.size()||a.nprim.size()!=a.nptr.size()||std::accumulate(a.ncomp.begin(),a.ncomp.end(),0)!=static_cast<int>(n))return fail("invalid_archive","Archive basis metadata is incomplete");
    const std::array<const std::vector<NboMatrix>*,2> matrix_lists{{&a.matrices,&d.matrices}};
    for(const auto* list:matrix_lists)for(const auto& m:*list)if(m.rows!=n||m.columns!=n||m.values.size()!=n*n||!std::all_of(m.values.begin(),m.values.end(),[](double x){return std::isfinite(x);}))return fail("invalid_matrix","Matrix dimensions or finite-value contract violated");
    if(w.source!=WavefunctionSource::Fchk)return fail("unsupported_source","Strict association currently requires Gaussian FCHK source coefficients");
    if(w.basis_count!=n||w.atoms.size()!=a.atoms.size())return fail("incompatible_dimensions","Atom/basis count differs");
    if(w.gaussian_ao_transform.size()!=n||w.ao_overlap.size()!=n*n)return fail("missing_metric","FCHK Gaussian AO transform and complete AO overlap are required");
    for(std::size_t i=0;i<a.atoms.size();++i){const auto& x=a.atoms[i];const auto& y=w.atoms[i];if(x.atomic_number!=y.atomic_number||std::abs(x.nuclear_charge-y.nuclear_charge)>1e-6)return fail("incompatible_atoms","Ordered atomic numbers or ECP effective nuclear charges differ");for(double e:{std::abs(x.x-y.x),std::abs(x.y-y.y),std::abs(x.z-y.z)})r.geometry_max_error_bohr=std::max(r.geometry_max_error_bohr,e);}
    if(r.geometry_max_error_bohr>3e-6)return fail("incompatible_geometry","Ordered geometry differs beyond .47 printed-coordinate precision; reordering unsupported");
    // Compare actual primitive contractions, respecting Gaussian's SP split in COV.
    std::size_t sh=0,ao=0;for(std::size_t k=0;k<a.ncomp.size();++k){const bool sp=a.ncomp[k]==4&&a.labels[ao]==1;const int parts=sp?2:1;for(int part=0;part<parts;++part){if(sh>=w.shells.size())return fail("incompatible_basis","Shell sequence differs");const auto& s=w.shells[sh++];const std::size_t nc=sp?(part==0?1:3):a.ncomp[k];const auto count=s.pure?2*s.angular_momentum+1:(s.angular_momentum+1)*(s.angular_momentum+2)/2;const auto l=sp?part:a.labels[ao]/100;if(count!=nc||s.angular_momentum!=l||s.atom_index+1!=static_cast<unsigned>(a.centers[ao])||s.primitive_count!=static_cast<unsigned>(a.nprim[k]))return fail("incompatible_basis","Ordered shell angular/center/primitive metadata differs");const std::vector<double>* coeff=l==0?&a.cs:l==1?&a.cp:l==2?&a.cd:l==3?&a.cf:l==4?&a.cg:nullptr;if(!coeff||coeff->size()!=a.exponents.size())return fail("unsupported_basis","Missing complete angular contraction array or angular momentum above g");for(std::size_t p=0;p<s.primitive_count;++p){const auto q=a.nptr[k]-1+p;const auto& prim=w.primitives.at(s.primitive_offset+p);if(std::abs(prim.exponent-a.exponents[q])>1e-7*std::max(1.0,std::abs(prim.exponent))||std::abs(prim.coefficient-(*coeff)[q])>2e-7*std::max(1.0,std::abs(prim.coefficient)))return fail("incompatible_basis","Primitive exponent or contraction coefficient differs");}if(sp&&part==0)++ao;else ao+=nc;}}
    if(sh!=w.shells.size()||ao!=n)return fail("incompatible_basis","Extra shell or AO entries");
    if(w.electron_counts_provenance==DataProvenance::Unavailable)return fail("missing_electronic_state","Explicit FCHK electron counts are required");
    const auto packed_count=n*(n+1)/2;
    const bool producer_total=w.total_density_provenance==DataProvenance::Producer&&w.total_density_packed.size()==packed_count;
    const bool producer_spin=w.spin_density_provenance==DataProvenance::Producer&&w.spin_density_packed.size()==packed_count;
    const bool density_supported_beta=a.open_shell&&w.orbital_occupation_model!=OrbitalOccupationModel::ExplicitSpin&&producer_total&&producer_spin;
    if(a.open_shell!=(w.orbital_occupation_model==OrbitalOccupationModel::ExplicitSpin)&&!density_supported_beta)return fail("incompatible_spin","Archive/FCHK spin representations differ and producer total/spin densities are unavailable");
    const auto* overlap=matrix(a.matrices,"OVERLAP",NboSpin::Total);if(!overlap)return fail("missing_overlap",".47 overlap required to verify AO convention");
    std::vector<std::size_t> src_to_internal(n),gcenter(n);std::vector<double> sg(n*n);for(std::size_t i=0;i<n;++i)src_to_internal[w.gaussian_ao_transform[i].source_index]=i;
    for(auto& shell:w.shells){const std::size_t count=shell.pure?2*shell.angular_momentum+1:(shell.angular_momentum+1)*(shell.angular_momentum+2)/2;for(std::size_t j=0;j<count;++j)gcenter[w.gaussian_ao_transform[shell.basis_offset+j].source_index]=shell.atom_index+1;}
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){auto ii=src_to_internal[i],jj=src_to_internal[j];sg[i*n+j]=w.ao_overlap[ii*n+jj]/(w.gaussian_ao_transform[ii].basis_scale*w.gaussian_ao_transform[jj].basis_scale);}
    const std::vector<NboSpin> spins=a.open_shell?std::vector<NboSpin>{NboSpin::Alpha,NboSpin::Beta}:std::vector<NboSpin>{NboSpin::Total};
    bool first=true;
    for(auto spin:spins){const auto* lc=matrix(a.matrices,"LCAOMO",spin);const auto* pd=matrix(a.matrices,"DENSITY",spin);if(!lc||!pd)return fail("missing_electronic_evidence","Complete .47 LCAOMO and DENSITY are required");
        std::vector<const MolecularOrbital*> mos(n,nullptr);for(auto& mo:w.orbitals)if((spin!=NboSpin::Beta&&mo.spin==Spin::Alpha)||(spin==NboSpin::Beta&&mo.spin==Spin::Beta)){if(mo.source_orbital_index>=n||mos[mo.source_orbital_index])return fail("ambiguous_canonical_identity","Canonical MO identity missing or repeated");mos[mo.source_orbital_index]=&mo;}
        const bool direct_coefficients=std::all_of(mos.begin(),mos.end(),[&](const MolecularOrbital* mo){return mo&&mo->gaussian_source_coefficients.size()==n;});
        const bool absent_coefficients=std::all_of(mos.begin(),mos.end(),[](const MolecularOrbital* mo){return !mo;});
        if(!direct_coefficients&&!(spin==NboSpin::Beta&&absent_coefficients&&density_supported_beta))return fail("missing_canonical_coefficients","All available canonical source MO columns are required; absent beta requires independent producer total/spin densities");
        NboCanonicalEvidence evidence;evidence.spin=spin;evidence.direct_fchk_coefficients=direct_coefficients;evidence.coefficient_source=direct_coefficients?"FCHK and archive LCAOMO":"archive LCAOMO only; direct FCHK beta coefficients unavailable";
        if(first){r.gaussian_row.resize(n);r.coefficient_scale.resize(n);std::set<std::size_t> used;
            for(std::size_t i=0;i<n;++i){std::vector<std::pair<std::size_t,double>> candidates;for(std::size_t g=0;g<n;++g){if(gcenter[g]!=static_cast<std::size_t>(a.centers[i]))continue;double xx=0,xy=0;for(std::size_t j=0;j<n;++j){xx+=lc->values[i*n+j]*lc->values[i*n+j];xy+=lc->values[i*n+j]*mos[j]->gaussian_source_coefficients[g];}if(xx<1e-18)continue;double scale=xy/xx,error=0;for(std::size_t j=0;j<n;++j)error=std::max(error,std::abs(scale*lc->values[i*n+j]-mos[j]->gaussian_source_coefficients[g]));if(std::abs(scale)>0.01&&std::abs(scale)<100&&error<3e-6)candidates.push_back({g,scale});}if(candidates.size()!=1||used.count(candidates.front().first))return fail("ambiguous_ao_mapping","Full canonical coefficients do not establish a unique same-source AO row mapping");r.gaussian_row[i]=candidates[0].first;r.coefficient_scale[i]=candidates[0].second;used.insert(candidates[0].first);}first=false;
        }
        std::vector<double> pn(n*n);
        const bool use_producer_density=producer_total&&(spin==NboSpin::Total||producer_spin);
        for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){const auto gi=r.gaussian_row[i],gj=r.gaussian_row[j];const auto si=r.coefficient_scale[i],sj=r.coefficient_scale[j];
            if(use_producer_density){const auto ii=src_to_internal[gi],jj=src_to_internal[gj];const auto hi=std::max(ii,jj),lo=std::min(ii,jj),k=hi*(hi+1)/2+lo;double p=w.total_density_packed[k];if(spin!=NboSpin::Total)p=(p+(spin==NboSpin::Alpha?1:-1)*w.spin_density_packed[k])*0.5;pn[i*n+j]=p/(w.gaussian_ao_transform[ii].coefficient_scale*w.gaussian_ao_transform[jj].coefficient_scale*si*sj);}
            else for(std::size_t k=0;k<n;++k)pn[i*n+j]+=mos[k]->occupation*mos[k]->gaussian_source_coefficients[gi]*mos[k]->gaussian_source_coefficients[gj]/(si*sj);
            r.overlap_max_error=std::max(r.overlap_max_error,std::abs(overlap->values[i*n+j]-sg[gi*n+gj]*si*sj));if(direct_coefficients)r.canonical_max_error=std::max(r.canonical_max_error,std::abs(si*lc->values[i*n+j]-mos[j]->gaussian_source_coefficients[gi]));}
        if(r.overlap_max_error>1e-5||r.canonical_max_error>3e-6)return fail("incompatible_ao_convention","Full canonical coefficients/metric disagree after established AO mapping");
        auto expected=a.density_is_bond_order?pn:multiply(multiply(overlap->values,pn,n),overlap->values,n);for(std::size_t i=0;i<n*n;++i)r.density_max_error=std::max(r.density_max_error,std::abs(expected[i]-pd->values[i]));if(r.density_max_error>2e-5)return fail("incompatible_density","Archive density does not reproduce FCHK canonical electron state");
        double electrons=0;for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j)electrons+=pn[i*n+j]*overlap->values[j*n+i];const double ne=spin==NboSpin::Alpha?w.alpha_electrons:spin==NboSpin::Beta?w.beta_electrons:w.alpha_electrons+w.beta_electrons;if(std::abs(electrons-ne)>2e-4)return fail("incompatible_electrons","Explicit electron trace differs");
        evidence.density_verified=true;evidence.detail=use_producer_density?"Independent FCHK producer density and electron trace agree with archive":"Density reconstructed from complete FCHK canonical coefficients agrees with archive";if(!direct_coefficients)evidence.detail+="; beta canonical coefficient comparison unavailable; NBOMO beta refers only to archive canonical columns";r.canonical_evidence.push_back(evidence);
        const auto* b=matrix(d.matrices,"AONBO",spin);const auto* t=matrix(d.matrices,"NBOMO",spin);
        if(b){auto sb=multiply(overlap->values,b->values,n);double orth=0;for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j){double v=0;for(std::size_t k=0;k<n;++k)v+=b->values[k*n+i]*sb[k*n+j];orth=std::max(orth,std::abs(v-(i==j?1:0)));}if(orth>2e-5)return fail("invalid_aonbo_metric","AONBO columns are not orthonormal in archive AO metric");
            auto psb=multiply(pn,sb,n);std::size_t count=0;for(auto& o:d.orbitals)if(o.spin==spin){if(o.id==0||o.id>n)return fail("invalid_nbo_identity","NBO column index outside matrix");double occ=0;for(std::size_t k=0;k<n;++k)occ+=sb[k*n+o.id-1]*psb[k*n+o.id-1];if(std::abs(occ-o.occupation)>4e-5)return fail("incompatible_nbo_output","Printed NBO occupation does not match complete AONBO/source density");++count;}if(count!=n)return fail("incomplete_nbo_output","One orbital record per full AONBO column is required");
            if(t){auto bt=multiply(b->values,t->values,n);double err=0;for(std::size_t i=0;i<n*n;++i)err=std::max(err,std::abs(bt[i]-lc->values[i]));if(err>2e-5)return fail("incompatible_nbomo","AONBO * NBOMO differs from archive canonical coefficients");}
        }else if(t)return fail("missing_aonbo","NBOMO alone cannot establish the NBO/source basis identity");
    }
    r.compatible=true;r.status=density_supported_beta?"strict_same_source_density_supported_beta":"strict_same_source";r.detail="Ordered atoms/effective charges/geometry, primitive basis, unique AO mapping, available complete canonical coefficients, overlap and explicit-spin density agree";if(density_supported_beta)r.detail+="; direct FCHK beta canonical coefficient verification is unavailable (archive beta canonical columns only)";d.association=r;return r;
}

Wavefunction make_nbo_wavefunction(const NboDataset& d,const Wavefunction& canonical) {
    if(!d.association.compatible||!d.archive)throw std::runtime_error("NBO rendering requires successful strict same-source association");
    // Recheck against the supplied canonical object: association is not a reusable
    // permission token for a subsequently changed wavefunction or dataset.
    NboDataset checked=d;if(!associate_nbo(checked,canonical).compatible)throw std::runtime_error("NBO association no longer matches canonical wavefunction");
    Wavefunction out;out.atoms=canonical.atoms;out.primitives=canonical.primitives;out.shells=canonical.shells;out.basis_count=canonical.basis_count;out.pure_d=canonical.pure_d;out.pure_f=canonical.pure_f;out.pure_g=canonical.pure_g;out.source_title="NBO localized orbitals (independent dataset)";out.ao_overlap=canonical.ao_overlap;out.ao_overlap_provenance=canonical.ao_overlap_provenance;
    const auto n=out.basis_count;std::vector<std::size_t> gaussian_to_internal(n);for(std::size_t i=0;i<n;++i)gaussian_to_internal[canonical.gaussian_ao_transform[i].source_index]=i;
    for(auto& x:d.orbitals){const auto* m=matrix(d.matrices,"AONBO",x.spin);if(!m||m->rows!=n||m->columns!=n||x.id==0||x.id>n)throw std::runtime_error("Complete AONBO matrix for every orbital/spin is required");MolecularOrbital mo;mo.energy_hartree=std::numeric_limits<double>::quiet_NaN();mo.occupation=x.occupation;mo.occupation_provenance=DataProvenance::Producer;mo.spin=x.spin==NboSpin::Beta?Spin::Beta:Spin::Alpha;mo.spin_provenance=DataProvenance::Producer;mo.spin_source_text=nbo_spin_name(x.spin);mo.source_orbital_index=x.id-1;mo.coefficients.resize(n);for(std::size_t i=0;i<n;++i){const auto internal=gaussian_to_internal[checked.association.gaussian_row[i]];mo.coefficients[internal]=m->values[i*n+x.id-1]*checked.association.coefficient_scale[i]*canonical.gaussian_ao_transform[internal].coefficient_scale;}out.orbitals.push_back(std::move(mo));}
    if(out.orbitals.empty())throw std::runtime_error("No NBO orbitals available");return out;
}

std::string serialize_nbo_json(const NboDataset& d) {
    std::ostringstream o;o<<std::setprecision(16);o<<"{\"schema\":\"cov.nbo.v1\",\"producer_version\":"<<quote(d.producer_version)<<",\"source\":";json_source(o,d.source);
    o<<",\"energy_semantics\":\"NBO and NAO energies are diagonal Fock elements, not canonical eigenvalues\",\"association\":{\"compatible\":"<<(d.association.compatible?"true":"false")<<",\"status\":"<<quote(d.association.status)<<",\"detail\":"<<quote(d.association.detail)<<",\"geometry_max_error_bohr\":"<<d.association.geometry_max_error_bohr<<",\"overlap_max_error\":"<<d.association.overlap_max_error<<",\"density_max_error\":"<<d.association.density_max_error<<",\"canonical_max_error\":"<<d.association.canonical_max_error<<",\"canonical_evidence\":";array(o,d.association.canonical_evidence,[&](const NboCanonicalEvidence& e){o<<"{\"spin\":"<<quote(nbo_spin_name(e.spin))<<",\"coefficient_source\":"<<quote(e.coefficient_source)<<",\"direct_fchk_coefficients\":"<<(e.direct_fchk_coefficients?"true":"false")<<",\"density_verified\":"<<(e.density_verified?"true":"false")<<",\"detail\":"<<quote(e.detail)<<'}';});o<<",\"gaussian_row_zero_based\":";array(o,d.association.gaussian_row,[&](std::size_t v){o<<v;});o<<",\"coefficient_scale\":";array(o,d.association.coefficient_scale,[&](double v){o<<v;});o<<"}";
    o<<",\"populations\":";array(o,d.populations,[&](const NboPopulation& x){o<<"{\"atom\":"<<x.atom<<",\"symbol\":"<<quote(x.symbol)<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"charge\":"<<x.charge<<",\"core\":"<<x.core<<",\"valence\":"<<x.valence<<",\"rydberg\":"<<x.rydberg<<",\"total\":"<<x.total<<",\"effective_core_electrons_derived\":";optional_number(o,x.effective_core_electrons);o<<",\"explicit_population_derived\":";optional_number(o,x.explicit_population);o<<",\"spin_density\":";optional_number(o,x.spin_density);o<<",\"source\":";json_source(o,x.source);o<<'}';});
    o<<",\"naos\":";array(o,d.naos,[&](const NboNao& x){o<<"{\"id\":"<<x.id<<",\"atom\":"<<x.atom<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"angular\":"<<quote(x.angular)<<",\"type\":"<<quote(x.type)<<",\"occupation\":"<<x.occupation<<",\"diagonal_fock_hartree\":";optional_number(o,x.energy_hartree);o<<",\"spin_density\":";optional_number(o,x.spin_density);o<<",\"source\":";json_source(o,x.source);o<<'}';});
    o<<",\"orbitals\":";array(o,d.orbitals,[&](const NboOrbital& x){o<<"{\"id\":"<<x.id<<",\"ordinal\":"<<x.ordinal<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"kind\":"<<quote(x.kind)<<",\"label\":"<<quote(x.label)<<",\"occupation\":"<<x.occupation<<",\"diagonal_fock_hartree\":";optional_number(o,x.diagonal_fock_hartree);o<<",\"energy_source\":";if(x.energy_source)json_source(o,*x.energy_source);else o<<"null";o<<",\"atoms\":";array(o,x.atoms,[&](std::size_t v){o<<v;});o<<",\"components\":";array(o,x.components,[&](const NboLocalComponent& c){o<<"{\"atom\":"<<c.atom<<",\"percent\":"<<c.percent<<",\"coefficient\":"<<c.coefficient<<",\"hybrid\":"<<quote(c.hybrid)<<",\"source\":";json_source(o,c.source);o<<'}';});o<<",\"source\":";json_source(o,x.source);o<<'}';});
    o<<",\"e2_sections\":";array(o,d.e2_sections,[&](const NboE2Section& x){o<<"{\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"printing_threshold\":";optional_number(o,x.printing_threshold);o<<",\"units\":"<<quote(x.units)<<",\"missing_reason\":"<<quote(x.missing_reason)<<",\"source\":";json_source(o,x.source);o<<'}';});
    o<<",\"e2\":";array(o,d.e2,[&](const NboE2& x){o<<"{\"donor\":"<<x.donor<<",\"acceptor\":"<<x.acceptor<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"value\":"<<x.value<<",\"units\":"<<quote(x.units)<<",\"energy_gap_hartree\":"<<x.energy_gap_hartree<<",\"fock_hartree\":"<<x.fock_hartree<<",\"printing_threshold\":";optional_number(o,x.printing_threshold);o<<",\"source\":";json_source(o,x.source);o<<'}';});
    o<<",\"wiberg\":";array(o,d.wiberg,[&](const NboWiberg& x){o<<"{\"atom_a\":"<<x.atom_a<<",\"atom_b\":"<<x.atom_b<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"value\":"<<x.value<<",\"source\":";json_source(o,x.source);o<<'}';});
    auto jm=[&](const NboMatrix& x){o<<"{\"kind\":"<<quote(x.kind)<<",\"spin\":"<<quote(nbo_spin_name(x.spin))<<",\"rows\":"<<x.rows<<",\"columns\":"<<x.columns<<",\"layout\":\"row-major\",\"values\":";array(o,x.values,[&](double v){o<<v;});o<<",\"source\":";json_source(o,x.source);o<<'}';};
    o<<",\"matrices\":";array(o,d.matrices,jm);o<<",\"archive\":";if(d.archive){const auto& a=*d.archive;o<<"{\"basis_count\":"<<a.basis_count<<",\"open_shell\":"<<(a.open_shell?"true":"false")<<",\"density_is_bond_order\":"<<(a.density_is_bond_order?"true":"false")<<",\"atoms\":";array(o,a.atoms,[&](const Atom& x){o<<"{\"atomic_number\":"<<x.atomic_number<<",\"effective_nuclear_charge\":"<<x.nuclear_charge<<",\"x_bohr\":"<<x.x<<",\"y_bohr\":"<<x.y<<",\"z_bohr\":"<<x.z<<'}';});o<<",\"source\":";json_source(o,a.source);o<<",\"centers\":";array(o,a.centers,[&](int v){o<<v;});o<<",\"labels\":";array(o,a.labels,[&](int v){o<<v;});o<<",\"ncomp\":";array(o,a.ncomp,[&](int v){o<<v;});o<<",\"nprim\":";array(o,a.nprim,[&](int v){o<<v;});o<<",\"nptr\":";array(o,a.nptr,[&](int v){o<<v;});o<<",\"exponents\":";array(o,a.exponents,[&](double v){o<<v;});o<<",\"contractions\":{\"s\":";array(o,a.cs,[&](double v){o<<v;});o<<",\"p\":";array(o,a.cp,[&](double v){o<<v;});o<<",\"d\":";array(o,a.cd,[&](double v){o<<v;});o<<",\"f\":";array(o,a.cf,[&](double v){o<<v;});o<<",\"g\":";array(o,a.cg,[&](double v){o<<v;});o<<"},\"matrices\":";array(o,a.matrices,jm);o<<'}';}else o<<"null";
    o<<",\"cmo_summaries\":";array(o,d.cmo_summaries,[&](const NboSource& x){json_source(o,x);});o<<",\"warnings\":";array(o,d.warnings,[&](const std::string& x){o<<quote(x);});o<<'}';return o.str();
}
} // namespace cov
