// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <cstring>

// -----------------------------
// Helpers to parse/format karyotypes
// -----------------------------
static std::vector<int> parse_kt(const std::string &s){
  std::vector<int> kt; kt.reserve(32);
  std::stringstream ss(s); std::string tok;
  while (std::getline(ss, tok, '.')) kt.push_back(std::stoi(tok));
  return kt;
}
static std::string kt_str(const std::vector<int> &kt){
  std::ostringstream oss;
  for (size_t i=0;i<kt.size();++i){ if(i) oss<<'.'; oss<<kt[i]; }
  return oss.str();
}

// -----------------------------
// Main
// -----------------------------
// [[Rcpp::export]]
Rcpp::List run_karyotype_neutral(Rcpp::NumericVector initial_counts_named,
                                 double rate,                // per-cell division rate
                                 double p_misseg,            // per-chromosome misseg prob
                                 double dt,
                                 int    n_steps,
                                 long long max_pop = 0,      // must be >0 here (fixed pool)
                                 double cull_keep = 0.1,     // keep fraction on cull (0..1]
                                 int    record_every = 1,
                                 int    seed = -1)
{
  if (max_pop <= 0) Rcpp::stop("max_pop must be > 0 (fixed preallocation only).");
  if (p_misseg < 0.0 || p_misseg > 1.0) Rcpp::stop("p_misseg must be in [0,1].");
  cull_keep = std::max(0.0, std::min(1.0, cull_keep));
  
  // RNG
  std::mt19937_64 rng(seed >= 0 ? (uint64_t)seed : std::random_device{}());
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::uniform_int_distribution<int> sign01(0,1);
  
  // Infer chromosome dimension K from the first named karyotype
  Rcpp::CharacterVector nm = initial_counts_named.names();
  if (nm.size()==0) Rcpp::stop("initial_counts_named must have names (karyotype strings).");
  std::vector<int> kt0 = parse_kt(Rcpp::as<std::string>(nm[0]));
  const int K = (int)kt0.size();
  if (K <= 0) Rcpp::stop("Could not infer karyotype length.");
  
  // Flat uint8_t population storage: rows [0..Npop-1] are live; rows >= Npop are garbage
  std::vector<unsigned char> pop((size_t)max_pop * (size_t)K, 0u);
  auto row_ptr = [&](long long r)->unsigned char* { return &pop[(size_t)r * (size_t)K]; };
  auto write_row_vec = [&](long long dst, const std::vector<int>& v){
    unsigned char* p = row_ptr(dst);
    for (int j=0;j<K;++j) p[j] = (unsigned char)v[j];
  };
  auto copy_row = [&](long long dst, long long src){
    std::memcpy(row_ptr(dst), row_ptr(src), (size_t)K);
  };
  auto swap_rows = [&](long long a, long long b){
    if (a==b) return;
    unsigned char* pa = row_ptr(a);
    unsigned char* pb = row_ptr(b);
    for (int j=0;j<K;++j) std::swap(pa[j], pb[j]);
  };
  
  // Build initial population up to max_pop (truncate extras if needed)
  long long Npop = 0;
  for (int i=0;i<initial_counts_named.size();++i){
    std::vector<int> kt = parse_kt(Rcpp::as<std::string>(nm[i]));
    if ((int)kt.size()!=K) Rcpp::stop("All karyotypes must have the same length.");
    long long c = (long long)std::llround((double)initial_counts_named[i]);
    if (c <= 0) continue;
    for (long long k=0; k<c && Npop<max_pop; ++k){
      write_row_vec(Npop, kt);
      ++Npop;
    }
    if (Npop >= max_pop) break;
  }
  if (Npop == 0) Rcpp::stop("Initial population is empty after truncation to max_pop.");
  
  // Recording helper: histogram of current live rows
  // Build histogram (just a named NumericVector)
  auto record_hist = [&](){
    std::unordered_map<std::string, long long> H; H.reserve((size_t)Npop*2+1);
    std::vector<int> tmp(K);
    for (long long r=0;r<Npop;++r){
      unsigned char* p = row_ptr(r);
      for (int j=0;j<K;++j) tmp[j] = (int)p[j];
      ++H[kt_str(tmp)];
    }
    Rcpp::NumericVector v(H.size());
    Rcpp::CharacterVector nms(H.size());
    size_t idx=0;
    for (auto &kv : H){
      v[idx] = (double)kv.second;
      nms[idx] = kv.first;
      ++idx;
    }
    v.attr("names") = nms;
    return v; // << just return the named vector
  };
  
  
  // Cull: shuffle all live rows, then truncate to floor(Npop * cull_keep) (at least 1 if Npop>0)
  auto cull_shuffle_truncate = [&](){
    if (Npop <= 0) return;
    // Fisher–Yates on [0..Npop-1]
    for (long long i=Npop-1; i>0; --i){
      std::uniform_int_distribution<long long> J(0, i);
      long long j = J(rng);
      swap_rows(i, j);
    }
    long long keep = (long long)std::floor((double)Npop * cull_keep);
    if (keep < 1) keep = 1;
    Npop = std::min<long long>(keep, max_pop);
  };
  
  // Division with balanced misseg (per-chromosome) and 3-case outcome
  enum class DivResult { AppendedNew, OverwriteOnly, RemovedParent };
  auto divide_balanced_3case = [&](long long parent_row, long long append_row) -> DivResult {
    unsigned char* p = row_ptr(parent_row);
    
    // build daughters as int scratch (to allow negatives during proposal)
    std::array<int, 64> d1, d2; // 64 > 22
    for (int j=0;j<K;++j){ d1[j] = (int)p[j]; d2[j] = (int)p[j]; }
    
    // balanced misseg per chromosome (no clamping)
    for (int j=0;j<K;++j){
      if (U(rng) < p_misseg){
        int s = sign01(rng) ? +1 : -1;
        d1[j] += s;
        d2[j] -= s;
      }
    }
    
    auto valid = [&](const std::array<int,64>& d){
      for (int j=0;j<K;++j) if (d[j] < 0 || d[j] > 8) return false;
      return true;
    };
    bool v1 = valid(d1), v2 = valid(d2);
    
    if (v1 && v2){
      // both valid: parent <- d1, append d2
      unsigned char* a = row_ptr(parent_row);
      unsigned char* b = row_ptr(append_row);
      for (int j=0;j<K;++j){ a[j] = (unsigned char)d1[j]; b[j] = (unsigned char)d2[j]; }
      return DivResult::AppendedNew;
    } else if (v1 ^ v2){
      // one valid: overwrite parent with valid, no append
      const auto& dv = v1 ? d1 : d2;
      unsigned char* a = row_ptr(parent_row);
      for (int j=0;j<K;++j) a[j] = (unsigned char)dv[j];
      return DivResult::OverwriteOnly;
    } else {
      // both invalid: remove the parent
      return DivResult::RemovedParent;
    }
  };
  
  // Output container
  Rcpp::List out;
  if (record_every>=1) out.push_back(record_hist(), "0");
  
  
  // Simulation loop
  for (int t=1; t<=n_steps; ++t){
    if (Npop == 0) break;
    
    // Poisson number of division events this step
    double lambda = rate * dt * (double)Npop;
    long long n_events = 0;
    if (lambda <= 0.0) {
      n_events = 0;
    } else if (lambda < 20.0){
      // Knuth
      double L = std::exp(-lambda), p = 1.0; n_events = -1;
      do { ++n_events; p *= U(rng); } while (p > L);
    } else {
      // Normal approximation (simple & fast)
      std::normal_distribution<double> N(lambda, std::sqrt(lambda));
      double x = std::max(0.0, N(rng));
      n_events = (long long)std::llround(x);
    }
    
    // If full and we might need to append, cull ahead of time
    if (Npop >= max_pop && cull_keep > 0.0) cull_shuffle_truncate();
    
    // Execute events
    for (long long e=0; e<n_events; ++e){
      if (Npop == 0) break;
      
      // Ensure space to append if the event yields a second daughter
      if (Npop >= max_pop && cull_keep > 0.0) {
        cull_shuffle_truncate();
        if (Npop == 0) break;
      }
      
      std::uniform_int_distribution<long long> I(0, Npop-1);
      long long parent = I(rng);
      
      // Attempt division with 3-case rule
      DivResult res = divide_balanced_3case(parent, Npop); // append target is Npop
      if (res == DivResult::AppendedNew){
        if (Npop < max_pop){
          ++Npop; // appended the second daughter
        } else {
          // No room to append; treat as overwrite-only fallback
          // (this occurs only if cull_keep==0 or max_pop==Npop post-cull)
        }
      } else if (res == DivResult::RemovedParent){
        // Parent dies: overwrite with last live row and shrink
        if (Npop == 1){
          Npop = 0;
          break;
        } else {
          swap_rows(parent, Npop-1);
          --Npop;
        }
      }
      // OverwriteOnly: nothing else to do
    }
    
    // Optional: if you want to cull exactly at/above cap after events
    if (Npop > max_pop && cull_keep > 0.0) cull_shuffle_truncate();
    
    // Record
    if (record_every>=1 && (t % record_every == 0 || t == n_steps)){
      out.push_back(record_hist(), std::to_string(t));
    }
    
    
    Rcpp::checkUserInterrupt();
  }
  
  return out;
}

