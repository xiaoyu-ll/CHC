#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<algorithm>
#include<cmath>
#include<map>
#include<vector>
#include<cstdlib>
#include <chrono>
#include <iomanip>
#include<stack>
#include<queue>
#include <cstring>
#include <thread>
#include <chrono>
#include <functional>  // std::greater
#include <unordered_set>
#include <ctime>
#include <numeric>   // for std::iota
#include <unordered_map>


static inline long long peak_rss_kb() {
    rusage r{};
    getrusage(RUSAGE_SELF, &r);
#if defined(__APPLE__)
    // macOS: ru_maxrss 是 bytes
    return (long long)r.ru_maxrss / 1024;
#else
    // Linux: ru_maxrss 是 KB
    return (long long)r.ru_maxrss;
#endif
}

static inline void print_peak(const char* tag) {
    std::cerr << "[MEM] " << tag << " peakRSS=" << peak_rss_kb() << " KB\n";
}

//#define maxn 2450000 //hb
//#define maxn 40000 //hc
//#define maxn 40000 //ch
//#define maxn 450000 //sb
//#define maxn 30000 //sc
//#define maxn 900 //cp
using namespace std;
const int vm=16000000;
const int em=4300000;
const int INF = 1000000000; 

struct Vertex
{
    bool flag=false;
    vector<int> inc;          // incident hyperedge ids
    vector<int> pos;          // pos[k] = vertex在 hyperedge[inc[k]].varr 中的位置（动态维护）
    int sigma=0;
};

struct Hyperedge
{
    vector<int> varr;         // 顶点数组（会 swap）
    vector<int> inc_idx;      // inc_idx[i] = 该 varr[i] 在 vertex[varr[i]].inc 中对应下标
    int alive = 0;            // 当前 active 前缀长度
    bool flag;
};



Vertex vertex[vm];
Hyperedge hyperedge[em];
int edgenum=0;
int nodenum=0;
int vmin=INF;
vector<int>indexv;
vector<int>indexc;
int maxeta=-1;



static inline uint64_t pack_uv_key(int u, int v) {
    if (u > v) std::swap(u, v);
    return (uint64_t)(uint32_t)u << 32 | (uint32_t)v;
}
static inline uint32_t cache_pack(uint32_t cnt, bool exact) {
    if (cnt > 0x7fffffffU) cnt = 0x7fffffffU;
    return (exact ? 0x80000000U : 0U) | cnt;
}
static inline bool cache_exact(uint32_t v) { return (v & 0x80000000U) != 0; }
static inline uint32_t cache_cnt(uint32_t v) { return v & 0x7fffffffU; }

static std::unordered_map<uint64_t, uint32_t> co_cache;
static const size_t MAX_CACHE = 5'000'000; // 按内存调：5M 一般够用





void readedge()//读数据
{
    //ifstream rda("dataset/hb8.txt");
    //ifstream rda("dataset/hb8.txt");
    ifstream rda("dataset/mathoverflow.txt");
    //ifstream rda("dataset/walmart.txt");
    //ifstream rda("dataset/trivago.txt");
    //ifstream rda("dataset/CA.txt");
    //ifstream rda("dataset/senate-bills.txt");
    //ifstream rda("dataset/house-bills.txt");
    //ifstream rda("dataset/stackoverflow.txt");
    //ifstream rda("dataset/amazon.txt");
    //ifstream rda("dataset/house-committees.txt");
    //ifstream rda("dataset/senate-committees.txt");
    //ifstream rda("dataset/contact-high.txt");
    //ifstream rda("dataset/contact-primary.txt");
    //ifstream rda("dataset/gptgene3.txt");
    //ifstream rda("dataset/test4.txt");
    //ifstream rda("dataset/tmathoverflow3.txt");
    //ifstream rda("dataset/email_eu3.txt");
    if(!rda)
    {
        cout<<"error!"<<endl;
        //std::exit(1);
    }
    string strline;
    int ei=0;
    int wmax=1;
    while(rda>>strline)//读取每一条超边并把顶点分隔开来
    {
        int ps=0,pt=0,i=0,tempv;
        hyperedge[ei].flag=true;
        while(i<strline.size())
        {
            if(strline[i]==',')
            {
                pt=i;
                string temps=strline.substr(ps,pt-ps);
                tempv=stoi(temps);
                // 原来：
// hyperedge[ei].varr.push_back(tempv);
// vertex[tempv].inc.push_back(ei);

int idx_in_v = (int)vertex[tempv].inc.size();              // 该 eid 在 v.inc 里的位置
vertex[tempv].inc.push_back(ei);
vertex[tempv].pos.push_back((int)hyperedge[ei].varr.size());// v 在该边 varr 中的位置（push 前 size）

hyperedge[ei].varr.push_back(tempv);
hyperedge[ei].inc_idx.push_back(idx_in_v);
                ps=i+1;
                if(tempv>nodenum)
                nodenum=tempv;
                if(tempv<vmin)
                vmin=tempv;
            }
            i++;
        }
        string temps=strline.substr(ps,i-ps);
        tempv=stoi(temps);
        // 原来：
// hyperedge[ei].varr.push_back(tempv);
// vertex[tempv].inc.push_back(ei);

int idx_in_v = (int)vertex[tempv].inc.size();              // 该 eid 在 v.inc 里的位置
vertex[tempv].inc.push_back(ei);
vertex[tempv].pos.push_back((int)hyperedge[ei].varr.size());// v 在该边 varr 中的位置（push 前 size）

hyperedge[ei].varr.push_back(tempv);
hyperedge[ei].inc_idx.push_back(idx_in_v);
        if(tempv>nodenum)
        nodenum=tempv;
        if(tempv<vmin)
        vmin=tempv;
        hyperedge[ei].alive = (int)hyperedge[ei].varr.size();
        ei++;
    }
    edgenum=ei;
    cout<<"read test edge successful!"<<endl;
    rda.close();
}



bool hasEdge(int u, int v, int eta)
{
    if (u == v) return true;
    if (u > v) std::swap(u, v);

    uint64_t key = pack_uv_key(u, v);

    auto it = co_cache.find(key);
    if (it != co_cache.end()) {
        uint32_t cv = it->second;
        uint32_t known = cache_cnt(cv);
        if (known >= (uint32_t)eta) return true;
        if (cache_exact(cv)) return false;  // 精确值仍 < eta
        // 否则是下界，不够，继续算
    }

    const auto &A = vertex[u].inc;
    const auto &B = vertex[v].inc;
    if (A.empty() || B.empty()) {
        if (co_cache.size() > MAX_CACHE) { co_cache.clear(); co_cache.reserve(MAX_CACHE); co_cache.max_load_factor(0.7f); }
        co_cache[key] = cache_pack(0, true);
        return false;
    }

    // 两指针求交：达到 eta 立即停止（此时得到的是下界）
    int i = 0, j = 0;
    uint32_t cnt = 0;
    while (i < (int)A.size() && j < (int)B.size()) {
        int x = A[i], y = B[j];
        if (x == y) {
            cnt++;
            if (cnt >= (uint32_t)eta) break; // ★ 提前停止
            i++; j++;
        } else if (x < y) i++;
        else j++;
    }

    bool exact = (i >= (int)A.size() || j >= (int)B.size()); // 扫到头才是精确

    if (co_cache.size() > MAX_CACHE) { co_cache.clear(); co_cache.reserve(MAX_CACHE); co_cache.max_load_factor(0.7f); }
    co_cache[key] = cache_pack(cnt, exact);

    return cnt >= (uint32_t)eta;
}




vector<int> collect_etas_noindex_all() {
    int maxDeg = 0;
    for (int v = 1; v <= nodenum; ++v) {
        maxDeg = max(maxDeg, (int)vertex[v].inc.size());
    }
    vector<int> etas;
    etas.reserve(maxDeg);
    for (int eta = maxDeg; eta >= 1; --eta) etas.push_back(eta);
    return etas;
}



static inline void build_candidates_noindex(int eta, vector<int>& cand) {
    cand.clear();
    // 你也可以 reserve 一个粗略值，比如 nodenum/10
    for (int v = 1; v <= nodenum; ++v) {
        if ((int)vertex[v].inc.size() >= eta) cand.push_back(v);
    }
}



static inline void color_sort_mcq(const vector<int>& cand, vector<int>& order, vector<int>& colors, int eta)
{
    // 启发式：先按 |inc| 降序（更容易先找到大团，提高剪枝）
    order = cand;
    sort(order.begin(), order.end(), [&](int a, int b){
        int da = (int)vertex[a].inc.size();
        int db = (int)vertex[b].inc.size();
        if (da != db) return da > db;
        return a < b;
    });

    vector<vector<int>> classes; // color classes
    classes.reserve(order.size());
    colors.assign(order.size(), 0);

    for (int idx = 0; idx < (int)order.size(); ++idx) 
    {
        int v = order[idx];

        int c = 0;
        for (; c < (int)classes.size(); ++c) 
        {
            bool conflict = false;
            for (int u : classes[c]) 
            {
                if (hasEdge(u, v, eta)) 
                { 
                    conflict = true; 
                    break; 
                }
            }
            if (!conflict) break;
        }
        if (c == (int)classes.size()) 
        classes.push_back({});
        classes[c].push_back(v);
        colors[idx] = c + 1; // 颜色号从1开始
    }

    // 按颜色类拼起来（MCQ惯例）
    vector<int> newOrder;
    vector<int> newColors;
    newOrder.reserve(order.size());
    newColors.reserve(order.size());

    for (int c = 0; c < (int)classes.size(); ++c) 
    {
        for (int v : classes[c]) 
        {
            newOrder.push_back(v);
            newColors.push_back(c + 1);
        }
    }
    order.swap(newOrder);
    colors.swap(newColors);
}

static void expand_all_max(vector<int>& cur, const vector<int>& cand, int eta, int& bestSize, vector<vector<int>>& all_max_cliques)
{
    if (cand.empty()) 
    {
        int sz = (int)cur.size();
        if (sz > bestSize) 
        {
            bestSize = sz;
            all_max_cliques.clear();
            all_max_cliques.push_back(cur);
        } 
        else if (sz == bestSize && sz > 0) 
        {
            all_max_cliques.push_back(cur);
        }
        return;
    }

    vector<int> order, colors;
    color_sort_mcq(cand, order, colors, eta);

    // 关键：这里用 "< bestSize" 才不会漏掉等于 bestSize 的其他解
    for (int i = (int)order.size() - 1; i >= 0; --i) {
        if ((int)cur.size() + colors[i] < bestSize) return;

        int v = order[i];

        // newCand = order[0..i-1] ∩ N(v)
        vector<int> newCand;
        newCand.reserve(i);
        for (int j = 0; j < i; ++j) {
            int u = order[j];
            if (hasEdge(u, v, eta)) newCand.push_back(u);
        }

        cur.push_back(v);

        if (newCand.empty()) {
            int sz = (int)cur.size();
            if (sz > bestSize) {
                bestSize = sz;
                all_max_cliques.clear();
                all_max_cliques.push_back(cur);
            } else if (sz == bestSize) {
                all_max_cliques.push_back(cur);
            }
        } else {
            // 一个弱剪枝：就算把 newCand 全加上也达不到 bestSize 就不递归
            if ((int)cur.size() + (int)newCand.size() >= bestSize)
                expand_all_max(cur, newCand, eta, bestSize, all_max_cliques);
        }

        cur.pop_back();
    }
}

static inline int upper_bound_by_coloring(const vector<int>& cand, int eta) 
{
    if (cand.empty()) return 0;
    vector<int> order, colors;
    color_sort_mcq(cand, order, colors, eta);
    int ub = 0;
    for (int c : colors) ub = max(ub, c);
    return ub;
}

static void expand_size_only(vector<int>& cur, const vector<int>& cand, int eta,int& bestSize)
{
    if (cand.empty()) 
    {
        bestSize = max(bestSize, (int)cur.size());
        return;
    }
    vector<int> order, colors;
    color_sort_mcq(cand, order, colors, eta);

    for (int i = (int)order.size() - 1; i >= 0; --i) {
        if ((int)cur.size() + colors[i] <= bestSize) return; // size剪枝(求size时用<=更狠)

        int v = order[i];
        vector<int> newCand;
        newCand.reserve(i);
        for (int j = 0; j < i; ++j) 
        {
            int u = order[j];
            // 便宜必败剪枝（可选但便宜）
            if (min((int)vertex[u].inc.size(), (int)vertex[v].inc.size()) < eta) 
            continue;
            if (hasEdge(u, v, eta)) newCand.push_back(u);
        }

        cur.push_back(v);
        if (newCand.empty()) 
        {
            bestSize = max(bestSize, (int)cur.size());
        } else 
        {
            // size下界剪枝
            if ((int)cur.size() + (int)newCand.size() > bestSize)
                expand_size_only(cur, newCand, eta, bestSize);
        }
        cur.pop_back();
    }
}

static void expand_collect_only(vector<int>& cur, const vector<int>& cand, int eta, int omega, vector<vector<int>>& all_max)
{
    if ((int)cur.size() + (int)cand.size() < omega) return; // 强 size 剪枝

    if (cand.empty()) 
    {
        if ((int)cur.size() == omega) all_max.push_back(cur);
        return;
    }

    vector<int> order, colors;
    color_sort_mcq(cand, order, colors, eta);

    for (int i = (int)order.size() - 1; i >= 0; --i) {
        if ((int)cur.size() + colors[i] < omega) return; // 关键：< 才不漏解

        int v = order[i];
        vector<int> newCand;
        newCand.reserve(i);
        for (int j = 0; j < i; ++j) {
            int u = order[j];
            if (min((int)vertex[u].inc.size(), (int)vertex[v].inc.size()) < eta) continue;
            if (hasEdge(u, v, eta)) newCand.push_back(u);
        }

        cur.push_back(v);
        expand_collect_only(cur, newCand, eta, omega, all_max);
        cur.pop_back();
    }
}



vector<vector<int>> findclique_allmax_two_phase(vector<int> vec, int eta, int curBestSize, long long curBestScore)
{
    // 必要过滤
    vec.erase(remove_if(vec.begin(), vec.end(), [&](int v){
        return (int)vertex[v].inc.size() < eta;
    }), vec.end());
    if (vec.empty()) return {};

    // 先用 coloring ub 做“双剪枝”
    int ub = upper_bound_by_coloring(vec, eta);
    if (ub < curBestSize) return {};
    if (1LL * ub * eta < curBestScore) return {};

    // Phase 1: 求 omega
    int omega = 0;
    vector<int> cur;
    expand_size_only(cur, vec, eta, omega);

    // 双剪枝：omega 自己也得过当前下界
    if (omega < curBestSize) return {};
    if (1LL * omega * eta < curBestScore) return {};

    // Phase 2: 枚举所有 size=omega 的最大团
    vector<vector<int>> all_max;
    cur.clear();
    expand_collect_only(cur, vec, eta, omega, all_max);

    // 去重（可选但建议）
    for (auto &cl : all_max) sort(cl.begin(), cl.end());
    sort(all_max.begin(), all_max.end());
    all_max.erase(unique(all_max.begin(), all_max.end()), all_max.end());

    return all_max;
}

vector<vector<int>> findmchc_noindex(
    const vector<int>& etas,
    vector<int>& etavec
){
    long long bestScore = 0;
    int bestSize = 0;
    vector<vector<int>> ans;
    etavec.clear();

    vector<int> cand;

    for (int eta : etas) {
        if(eta<3)
        break;
        // 一个很便宜的上界：如果 cand 最大可能 size * eta <= bestScore 就可跳
        build_candidates_noindex(eta, cand);
        if(cand.size()<3)
        continue;
        int ub = (int)cand.size();
        if (ub < 2) continue;
        if (1LL * ub * eta < bestScore) continue;

        // 精确最大团（两阶段：omega + all max）
        auto cliques = findclique_allmax_two_phase(cand, eta, bestSize, bestScore);
        if (cliques.empty()) continue;

        int omega = (int)cliques[0].size();
        if(omega<3)
        continue;
        long long score = 1LL * omega * eta;

        if (score > bestScore || (score == bestScore && omega > bestSize)) {
            bestScore = score;
            bestSize = omega;
            ans = cliques;
            etavec.assign(ans.size(), eta);
        } else if (score == bestScore && omega == bestSize) {
            // append
            for (auto &cl : cliques) {
                ans.push_back(cl);
                etavec.push_back(eta);
            }
        }

        // 可选：对 ans 去重（和你现在一样）
        for (auto &cl : ans) sort(cl.begin(), cl.end());
        sort(ans.begin(), ans.end());
        // 同步 etavec 去重比较麻烦：baseline里你也可以先不去重，
        // 或者改为 (eta, clique) pair 统一去重
    }

    return ans;
}

int main()
{
    readedge();
    print_peak("after readedge");
    auto start = std::chrono::high_resolution_clock::now();
    auto t0 = std::chrono::high_resolution_clock::now();
    auto etas = collect_etas_noindex_all();
    auto t2 = std::chrono::high_resolution_clock::now();
    double elapsed1 =std::chrono::duration<double>(t2- t0).count();
    std::cout << "colletc eta:"<<elapsed1 << " seconds\n";
    vector<int> etavec;
    auto res = findmchc_noindex(etas, etavec);
    cout<<endl;
    print_peak("after findmchc");
    for(int i=0;i<res.size();i++)
    {
        vector<int>clique=res[i];
        cout<<"clique "<<i<<": ";
        for(int j=0;j<clique.size();j++)
        {
            cout<<clique[j]<<" ";
        }
        cout<<", eta="<<etavec[i]<<", size="<<res[i].size()<<endl;
    }

    //for(int i=0;i<indexc.size();i++)
    //cout<<"("<<i+1<<", "<<indexc[i]<<") ";
    //cout<<endl;

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed2 =std::chrono::duration<double>(end - start).count();
    std::cout << elapsed2 << " seconds\n";
    //std::this_thread::sleep_for(std::chrono::seconds(300)); // 暂停3秒
    return 0;

}
