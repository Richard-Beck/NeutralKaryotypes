// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include <numeric>
#include <algorithm>

struct VecHash {
  size_t operator()(const std::vector<int>& v) const {
    size_t h = v.size();
    for (int x : v) h ^= std::hash<int>{}(x) + 0x9e3779b9 + (h<<6) + (h>>2);
    return h;
  }
};
using Pop = std::unordered_map<std::vector<int>, long long, VecHash>;

static inline std::vector<int> parse_kt(const std::string& s, int expect=-1){
  std::vector<int> v; std::stringstream ss(s); std::string t;
  while(std::getline(ss,t,'.')){ if(t.empty()) return {}; int x=std::stoi(t); if(x<0) return {}; v.push_back(x); }
  if(expect>0 && (int)v.size()!=expect) return {}; return v;
}
static inline std::string kt_str(const std::vector<int>& v){
  std::string s; s.reserve(v.size()*2);
  for(size_t i=0;i<v.size();++i){ s+=std::to_string(v[i]); if(i+1<v.size()) s+='.'; }
  return s;
}

// one division with k missegregations (drawn outside)
// returns {d1,d2}; empty vector means invalid daughter (<=0 copy of a type)
// nTot = sum(parent)
// rng provided
static inline std::pair<std::vector<int>,std::vector<int>>
  daughters_misseg(int k, const std::vector<int>& p, int nTot, std::mt19937& rng){
    if (k == 0) return {p, p};
    
    std::vector<int> pool; pool.reserve(nTot);
    for (size_t i = 0; i < p.size(); ++i)
      for (int c = 0; c < p[i]; ++c) pool.push_back((int)i);
    
    std::vector<int> d1 = p, d2 = p, pick; pick.reserve(k);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    
    std::sample(pool.begin(), pool.end(), std::back_inserter(pick), k, rng);
    for (int pos : pick) {
      if (U(rng) < 0.5) { ++d1[pos]; --d2[pos]; } else { --d1[pos]; ++d2[pos]; }
    }
    for (size_t i = 0; i < p.size(); ++i) {
      if (d1[i] < 0) d1.clear();
      if (d2[i] < 0) d2.clear();
      if (d1.empty() && d2.empty()) break;
    }
    return {d1, d2};
  }


// [[Rcpp::export]]
Rcpp::List run_karyotype_neutral(Rcpp::NumericVector initial_counts_named,
                                 double rate,                // per-cell division rate
                                 double p_misseg,            // per-chromosome misseg prob
                                 double dt,
                                 int    n_steps,
                                 long long max_pop = 0,      // 0 = no cap
                                 double cull_keep = 0.1,     // fraction to keep if capped
                                 int    record_every = 1,
                                 int    seed = -1) {
  // build initial population & dimension
  Pop pop; pop.reserve(initial_counts_named.size());
  int K = -1;
  Rcpp::CharacterVector nm = initial_counts_named.names();
  for(int i=0;i<initial_counts_named.size();++i){
    std::vector<int> cn = parse_kt(Rcpp::as<std::string>(nm[i]));
    if(cn.empty()) continue;
    if(K==-1) K=(int)cn.size(); else if((int)cn.size()!=K) Rcpp::stop("Mixed karyotype lengths.");
    long long c = std::llround((double)initial_counts_named[i]); if(c>0) pop[cn]+=c;
  }
  if(pop.empty()) return Rcpp::List::create();
  
  std::mt19937 rng; if(seed<0){ std::random_device rd; rng.seed(rd()); } else rng.seed((unsigned)seed);
  std::poisson_distribution<long long> Pois;
  Rcpp::List out;
  
  auto record = [&](const std::string& key){
    Rcpp::NumericVector cnt; Rcpp::CharacterVector nmv;
    cnt.attr("names") = nmv;
    for(auto &kv: pop){ cnt.push_back((double)kv.second); nmv.push_back(kt_str(kv.first)); }
    if(cnt.size()>0) cnt.names()=nmv; out.push_back(cnt, key);
  };
  if(record_every>=1) record("0");
  
  for(int t=1;t<=n_steps;++t){
    if(pop.empty()) break;
    
    Pop delta; delta.reserve(pop.size()*2);
    std::vector<std::pair<std::vector<int>, long long>> now; now.reserve(pop.size());
    for(auto &kv: pop) if(kv.second>0) now.emplace_back(kv.first, kv.second);
    
    for(auto &pc: now){
      const auto& parent = pc.first; long long n = pc.second;
      int nTot = std::accumulate(parent.begin(), parent.end(), 0);
      if(nTot==0 || n==0) continue;
      
      // Poisson divisions
      double mu = rate * dt * n; if (mu < 0) mu = 0;
      Pois.param(std::poisson_distribution<long long>::param_type(mu));
      
      long long n_div = std::min(n, Pois(rng));        // cannot divide more than present
      if(n_div==0) continue;
      
      delta[parent] -= n_div;                          // parents consumed
      std::binomial_distribution<int> Bin(nTot, p_misseg);
      
      for(long long i=0;i<n_div;++i){
        int kerr = (p_misseg>0 && nTot>0) ? Bin(rng) : 0;
        auto ds = daughters_misseg(kerr, parent, nTot, rng);
        if(!ds.first.empty())  delta[ds.first]  += 1;
        if(!ds.second.empty()) delta[ds.second] += 1;
      }
    }
    
    // apply changes
    long long tot=0;
    for(auto &kv: delta) pop[kv.first] += kv.second;
    for(auto it=pop.begin(); it!=pop.end();){
      if(it->second<=0) it = pop.erase(it); else { tot+=it->second; ++it; }
    }
    
    // cull if needed (neutral, uniform binomial downsample)
    if(max_pop>0 && tot>max_pop){
      std::binomial_distribution<long long> B;
      Pop kept; kept.reserve(pop.size());
      for(auto &kv: pop){
        // Binomial downsampling during cull
        double keep = std::max(0.0, std::min(1.0, cull_keep));
        B.param(std::binomial_distribution<long long>::param_type(kv.second, keep));
        long long k=B(rng); if(k>0) kept[kv.first]=k;
      }
      pop.swap(kept);
    }
    
    if(record_every>=1 && (t%record_every==0 || t==n_steps)) record(std::to_string(t));
    Rcpp::checkUserInterrupt();
  }
  return out;
}
