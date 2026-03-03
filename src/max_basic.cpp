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
const int vm=74000;
const int em=5500;
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


int buildindex()
{
    // P: 候选集合（未 ORDERED，且与 ORDERED 顶点共现过）
    vector<int> P;
    P.reserve(100000);

    // inP: 防止重复入队（你原来用 sort+unique；这里 O(1) 去重）
    static vector<unsigned char> inP;
    if ((int)inP.size() < nodenum + 1) inP.assign(nodenum + 1, 0);
    else fill(inP.begin(), inP.begin() + (nodenum + 1), 0);

    // 临时计数：用于计算某个 x 与每个 ORDERED 顶点的共现次数
    static vector<int> tmpCnt;
    static vector<int> touched;
    if ((int)tmpCnt.size() < nodenum + 1) tmpCnt.assign(nodenum + 1, 0);
    touched.reserve(200000);

    // init
    indexv.clear(); indexc.clear();
    indexv.reserve(nodenum);
    indexc.reserve(nodenum);

    for (int v = 1; v <= nodenum; ++v) {
        vertex[v].flag = false; // flag=true 表示 ORDERED
        vertex[v].sigma = 0;
    }

    int picked = 0;
    int wmax = 0;
    int next_unpicked = 1;

    auto pick_smallest_unordered = [&]() -> int {
        while (next_unpicked <= nodenum && vertex[next_unpicked].flag) ++next_unpicked;
        return (next_unpicked <= nodenum) ? next_unpicked : -1;
    };

    // 计算 sigma(x)=max_{y ORDERED} co(x,y)，不预计算 pair
    auto recompute_sigma_against_ordered = [&](int x) -> int {
        touched.clear();
        int best = 0;

        // 扫 x 的 incident 超边，统计与 ORDERED 顶点的共现次数
        for (int eid : vertex[x].inc) {
            auto &arr = hyperedge[eid].varr;
            for (int y : arr) {
                if (!vertex[y].flag) continue; // 只关心 ORDERED 顶点
                if (y == x) continue;
                if (tmpCnt[y] == 0) touched.push_back(y);
                int c = ++tmpCnt[y];
                if (c > best) best = c;
            }
        }

        // reset
        for (int y : touched) tmpCnt[y] = 0;
        return best;
    };

    while (picked < nodenum)
    {
        int u = -1;
        int bestSig = -1;
        int bestPos = -1;

        if (P.empty())
        {
            u = pick_smallest_unordered();
            if (u == -1) break;

            // 与伪码一致：P 为空时取最小 ID，sigma=0
            bestSig = 0;
        }
        else
        {
            // 每轮对 P 中每个 x 重新算 sigma(x)，然后 argmax
            for (int i = 0; i < (int)P.size(); ++i) {
                int x = P[i];
                if (vertex[x].flag) continue; // 保险：已 ORDERED 的不考虑

                int sx = recompute_sigma_against_ordered(x);
                vertex[x].sigma = sx; // 可存可不存

                // argmax，tie-break: sigma 大优先；sigma 相同选更小 id
                if (sx > bestSig || (sx == bestSig && x < u)) {
                    bestSig = sx;
                    u = x;
                    bestPos = i;
                }
            }

            // 从 P 移除 u
            if (bestPos != -1) {
                int last = P.back();
                P[bestPos] = last;
                P.pop_back();
                inP[u] = 0;
            }
        }

        // commit u
        vertex[u].flag = true;
        indexv.push_back(u);
        indexc.push_back(bestSig);  // 这一轮选中时的 sigma(u)
        wmax = max(wmax, bestSig);
        ++picked;

        // 把与 u 共现的未 ORDERED 顶点加入 P（对应伪码“for v in N(u) add to P”）
        for (int eid : vertex[u].inc) {
            auto &arr = hyperedge[eid].varr;
            for (int v : arr) {
                if (v == u) continue;
                if (vertex[v].flag) continue;
                if (!inP[v]) {
                    inP[v] = 1;
                    P.push_back(v);
                }
            }
        }
    }

    return wmax;
}

vector<int> collect_etas_from_indexc() 
{
    vector<int> etas = indexc;
    etas.push_back(1);
    sort(etas.begin(), etas.end());
    etas.erase(unique(etas.begin(), etas.end()), etas.end());
    sort(etas.begin(), etas.end(), greater<int>());
    return etas;
}

struct DSU {
    vector<int> p, sz;
    vector<char> act;
    DSU(int n=0){ init(n); }
    void init(int n){
        p.resize(n);
        sz.assign(n,1);
        act.assign(n,0);
        iota(p.begin(), p.end(), 0);
    }
    int find(int x){ return p[x]==x? x : p[x]=find(p[x]); }
    int unite(int a,int b){
        a=find(a); b=find(b);
        if(a==b) return sz[a];
        if(sz[a]<sz[b]) swap(a,b);
        p[b]=a;
        sz[a]+=sz[b];
        return sz[a];
    }
};

static inline void precompute_L_and_sufUB(const vector<int>& indexc, const vector<int>& etas, unordered_map<int,int>& L_of_eta, vector<long long>& sufUB)
{
    int n = (int)indexc.size();
    vector<int> pos(n);
    iota(pos.begin(), pos.end(), 0);

    // sort positions by indexc descending
    sort(pos.begin(), pos.end(), [&](int a, int b){
        return indexc[a] > indexc[b];
    });

    DSU dsu(n);
    int maxSeg = 0;
    int pptr = 0;

    // activate positions with indexc[pos] >= eta
    for (int eta : etas) 
    {
        while (pptr < n && indexc[pos[pptr]] >= eta) 
        {
            int x = pos[pptr++];
            dsu.act[x] = 1;
            dsu.p[x] = x;
            dsu.sz[x] = 1;

            if (x-1 >= 0 && dsu.act[x-1]) maxSeg = max(maxSeg, dsu.unite(x, x-1));
            if (x+1 < n && dsu.act[x+1])  maxSeg = max(maxSeg, dsu.unite(x, x+1));
            maxSeg = max(maxSeg, 1);
        }
        L_of_eta[eta] = maxSeg;  // max interval length under threshold eta
    }

    // sufUB[i] = max_{j>=i} etas[j] * L(etas[j])
    int m = (int)etas.size();
    sufUB.assign(m, 0);
    for (int i = m-1; i >= 0; --i) 
    {
        long long val = 1LL * etas[i] * (long long)L_of_eta[etas[i]];
        sufUB[i] = (i == m-1) ? val : max(val, sufUB[i+1]);
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

vector<vector<int>> findmchc(const vector<int> &etas, const vector<long long> &sufUB, const unordered_map<int,int> &L_of_eta, vector<int> &vec)
{
    long long bestScore = 0;
    int bestSize = 0;
    vector<vector<int>> mchc;
    long long best = 0;                          // ★ 用 long long
    int cli = 0;

    int n = (int)indexv.size();                  // indexv/indexc 都是 CHC order, size = nodenum

    for (int idxEta = 0; idxEta < (int)etas.size(); ++idxEta)
    {
        int eta = etas[idxEta];
        if(eta<3)
        break;
        //if(eta%10==1)
        //cout<<"eta="<<eta<<", ";

        // ★ 全局提前终止：后面所有更小 eta 的上界都不可能超过 best
        if (idxEta + 1 < (int)etas.size() && sufUB[idxEta + 1] <= best) break;

        // 生成 [l,r] ranges：indexc[pos] >= eta 的连续段
        vector<int> veclr;
        bool inSeg = false;
        for (int i = 0; i < n; ++i)              // ★ 0-based，别用 1..nodenum 访问 indexc[i]
        {
            if (!inSeg && indexc[i] >= eta) {
                inSeg = true;
                veclr.push_back(i-1);              // start
            }
            if (inSeg && indexc[i] < eta) {
                inSeg = false;
                veclr.push_back(i-1);            // end
            }
        }
        if (inSeg) veclr.push_back(n-1);

        // 遍历所有 ranges
        for (int t = 0; t + 1 < (int)veclr.size(); t += 2)
        {
            int l = veclr[t], r = veclr[t+1];
            int len = r - l + 1;
            if (len < 3) continue;

            // ★ 区间级剪枝不变
            if (1LL * len * eta < bestScore) continue;

            vector<int> ceta;
            ceta.reserve(len);
            //for (int i = l; i <= r; ++i) ceta.push_back(indexv[i]);
            for (int i = l; i <= r; ++i) 
            {
                int v = indexv[i];
                if ((int)vertex[v].inc.size() >= eta) 
                {
                    ceta.push_back(v);
                    //cout<<v<<", ";
                }
            }
            //cout<<ceta.size()<<endl;

            //vector<vector<int>> cliques = findclique(ceta, eta);
            //if (cliques.empty()) continue;
            auto cliques = findclique_allmax_two_phase(ceta, eta, bestSize, bestScore);
            if (cliques.empty()) 
            continue;
            
            int tsize = (int)cliques[0].size();
            long long score = 1LL * tsize * eta;
            if (score > bestScore || (score == bestScore && tsize > bestSize)) 
            {
                bestScore = score;
                bestSize = tsize;
                mchc.clear();
                vec.clear();
                cli = 0;
                for (auto &cl : cliques) mchc.push_back(cl);
            } 
            else if (score == bestScore && tsize == bestSize) 
            {
                for (auto &cl : cliques) mchc.push_back(cl);
            }

        }

        for (auto &cl : mchc) sort(cl.begin(), cl.end());
        sort(mchc.begin(), mchc.end());
        mchc.erase(unique(mchc.begin(), mchc.end()), mchc.end());

        int chcs = (int)mchc.size();
        for (int j = cli; j < chcs; ++j) vec.push_back(eta);
        cli = chcs;
    }
    return mchc;
}

int main()
{
    readedge();
    print_peak("after readedge");
    auto start = std::chrono::high_resolution_clock::now();
    auto t0 = std::chrono::high_resolution_clock::now();
    int comax = buildindex();
    cout<<"eta max="<<comax<<endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed =std::chrono::duration<double>(t1- t0).count();
    std::cout << "build index:"<<elapsed << " seconds\n";
    print_peak("after buildindex");
    vector<int> etas = collect_etas_from_indexc(); // 再从 indexc 提取阈值
    auto t2 = std::chrono::high_resolution_clock::now();
    double elapsed1 =std::chrono::duration<double>(t2- t1).count();
    std::cout << "colletc eta:"<<elapsed1 << " seconds\n";
    // ★ 预处理 L(eta) 和 sufUB
    unordered_map<int,int> L_of_eta;
    vector<long long> sufUB;
    precompute_L_and_sufUB(indexc, etas, L_of_eta, sufUB);
    vector<int> etavec;
    vector<vector<int>> res = findmchc(etas, sufUB, L_of_eta, etavec);
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
