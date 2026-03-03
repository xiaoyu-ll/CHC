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

struct Vertex //1791488
{
    bool flag;
    vector<int>inc;
    int sigma;
};
struct Hyperedge//1735400 
{
    vector<int> varr;
    bool flag;
};

Vertex vertex[vm];
Hyperedge hyperedge[em];
int edgenum=0;
int nodenum=0;
int vmin=INF;
vector<int>indexv;
vector<int>indexc;


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
    //ifstream rda("dataset/walmart.txt");
    //ifstream rda("dataset/trivago.txt");
    //ifstream rda("dataset/CA.txt");
    //ifstream rda("dataset/senate-bills.txt");
    //ifstream rda("dataset/house-bills.txt");
    //ifstream rda("dataset/stackoverflow.txt");
    //ifstream rda("dataset/amazon.txt");
    //ifstream rda("dataset/house-committees.txt");
    ifstream rda("dataset/senate-committees.txt");
    //ifstream rda("dataset/contact-high.txt");
    //ifstream rda("dataset/contact-primary.txt");
    //ifstream rda("dataset/mathoverflow.txt");
    //ifstream rda("dataset/wiki_topcats.txt");
    //ifstream rda("dataset/gptgene3.txt");
    //ifstream rda("dataset/test4.txt");
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
        vector<int> nodes;
        hyperedge[ei].flag=true;
        while(i<strline.size())
        {
            if(strline[i]==',')
            {
                pt=i;
                string temps=strline.substr(ps,pt-ps);
                tempv=stoi(temps);
                hyperedge[ei].varr.push_back(tempv);
                vertex[tempv].flag=false;
                vertex[tempv].sigma=0;
                vertex[tempv].inc.push_back(ei);
                nodes.push_back(tempv);
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
        hyperedge[ei].varr.push_back(tempv);
        vertex[tempv].flag=false;
        vertex[tempv].sigma=0;
        vertex[tempv].inc.push_back(ei);
        nodes.push_back(tempv);
        if(tempv>nodenum)
        nodenum=tempv;
        if(tempv<vmin)
        vmin=tempv;
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


struct NodeKey {
    int sigma;
    int v;
    bool operator<(const NodeKey& o) const {
        if (sigma != o.sigma) return sigma < o.sigma; // max-heap
        return v > o.v; // tie: smaller id first
    }
};

int buildindex_heap_exact()
{
    static vector<int> tmpCnt;
    static vector<int> touched;
    if ((int)tmpCnt.size() < nodenum + 1) tmpCnt.assign(nodenum + 1, 0);
    touched.reserve(100000);

    // init
    indexv.clear(); indexc.clear();
    for (int i = 1; i <= nodenum; ++i) { vertex[i].flag = false; vertex[i].sigma = 0; }

    priority_queue<NodeKey> pq;
    int picked = 0;
    int wmax = 0;

    auto pick_any_unflagged = [&](){
        for (int i = 1; i <= nodenum; ++i) if (!vertex[i].flag) return i;
        return -1;
    };

    while (picked < nodenum)
    {
        int uu = -1;

        // get best from heap
        while (!pq.empty()) {
            auto top = pq.top(); pq.pop();
            int v = top.v;
            if (vertex[v].flag) continue;
            if (top.sigma != vertex[v].sigma) continue; // stale
            uu = v;
            break;
        }
        if (uu == -1) uu = pick_any_unflagged(); // all remaining have sigma=0 and never pushed

        // commit
        vertex[uu].flag = true;
        indexv.push_back(uu);
        indexc.push_back(vertex[uu].sigma);
        wmax = max(wmax, vertex[uu].sigma);

        // count co(uu, v)
        touched.clear();
        for (int eid : vertex[uu].inc) {
            auto &arr = hyperedge[eid].varr;
            for (int v : arr) {
                if (v == uu) continue;
                if (vertex[v].flag) continue;
                if (tmpCnt[v] == 0) touched.push_back(v);
                tmpCnt[v] += 1;
            }
        }

        // update sigma[v] and push to heap
        for (int v : touched) {
            int cand = tmpCnt[v];
            tmpCnt[v] = 0;
            if (cand > vertex[v].sigma) {
                vertex[v].sigma = cand;
                pq.push({vertex[v].sigma, v});
            }
        }

        ++picked;
    }
    return wmax;
}



// cur: 当前 clique 中的点
// cand: 还能尝试加入 clique 的候选点（都在 vec 里面）
void dfs_clique(vector<int> &best_clique, vector<vector<int>> &all_max_cliques, vector<int> &cur, vector<int> cand, int eta)
{
    // 剪枝：就算把 cand 全部加进来，大小也不会超过当前最优
    if (cur.size() + cand.size() < best_clique.size())
        return;

    bool extended = false;

    while (!cand.empty())
    {
        if (cur.size() + cand.size() < best_clique.size())
        return;

        int v = cand.back();
        cand.pop_back();

        // 检查 v 是否和当前 clique 中所有点都连边
        bool ok = true;
        for (int u : cur)
        {
            if (!hasEdge(u, v, eta)) { ok = false; break; }
        }
        if (!ok) continue;

        extended = true;

        // newCand = cand ∩ N(v)
        vector<int> newCand;
        newCand.reserve(cand.size());
        for (int u : cand)
            if (hasEdge(u, v, eta))
                newCand.push_back(u);

        cur.push_back(v);
        dfs_clique(best_clique, all_max_cliques, cur, newCand, eta);
        cur.pop_back();
    }

    // maximal clique
    if (!extended)
    {
        int sz = (int)cur.size();
        int bestsz = (int)best_clique.size();

        if (sz > bestsz)
        {
            best_clique = cur;

            // 发现更大的最大 clique：清空并重置全局答案集
            all_max_cliques.clear();
            all_max_cliques.push_back(cur);
        }
        else if (sz == bestsz && sz > 0)
        {
            // 发现同样大小的最大 clique：加入全局答案集（可选去重）
            all_max_cliques.push_back(cur);
        }
    }
}

vector<vector<int>> findclique(vector<int> vec, int eta)
{
    vector<vector<int>> all_max_cliques;
    vector<int> best;
    vector<int> cur;
    dfs_clique(best, all_max_cliques, cur, vec, eta);
    return all_max_cliques;
}

vector<vector<int>> findmchc(int eta, int comax)
{
    vector<vector<int>> mchc;
    //vector<int>mchc;
    vector<int>veclr;
    int mchcs=-1;
    if(eta>comax)
    return mchc;
    bool flag=false;
    for(int i=0;i<(int)indexc.size();i++)
    {
        if(!flag&&indexc[i]>=eta)
        {
            flag=true;
            veclr.push_back(i-1);
        }
        if(flag&&indexc[i]<eta)
        {
            flag=false;
            veclr.push_back(i-1);
        }
    }
    int veclrs=(int)veclr.size();
    if(veclrs%2==1)
    veclr.push_back(indexv.size()-1);
    while(!veclr.empty())
    {
        int l=veclr[0];
        int r=veclr[1];
        veclr.erase(veclr.begin());
        veclr.erase(veclr.begin());
        if((r-l+1)<3)
        continue;
        if(mchcs!=-1&&(r-l+1)<mchcs)
        continue;
        vector<int>ceta;
        for(int i=l;i<=r;i++)
        {
            ceta.push_back(indexv[i]);
        }
        //ALL_MAX_CLIQUES.clear();                 // ★ 新增：清空本区间结果
        vector<vector<int>> cliques = findclique(ceta, eta); // 保留：仍然触发 dfs
        if(cliques.empty())
        continue;
        int temps = (int)cliques[0].size();
        if (temps > mchcs)
        {
            mchcs = temps;
            mchc.clear();
            // ★ 新增：把该区间所有最大 clique 都放进去
            for (auto &cl : cliques) mchc.push_back(cl);
        }
        else if (temps == mchcs)
        {
            // ★ 新增：同大小也要把所有最大 clique 都并进来
            for (auto &cl : cliques) mchc.push_back(cl);
        }
    }
    for (auto &cl : mchc) sort(cl.begin(), cl.end());
    sort(mchc.begin(), mchc.end());
    mchc.erase(unique(mchc.begin(), mchc.end()), mchc.end());
    return mchc;
}

int main()
{
    readedge();
    auto t0 = chrono::high_resolution_clock::now();
    int comax=buildindex_heap_exact();
    auto t1 = chrono::high_resolution_clock::now();
    cout << "buildindex: " << chrono::duration<double>(t1-t0).count() << " s\n";
    int eta=4;
    cout<<"eta="<<eta<<endl;
    cout<<"comax="<<comax<<endl;
    vector<vector<int>> res=findmchc(eta, comax);
    for(int i=0;i<res.size();i++)
    {
        vector<int>clique=res[i];
        cout<<"clique "<<i<<": "<<endl;
        for(int j=0;j<clique.size();j++)
        {
            cout<<clique[j]<<" ";
        }
        cout<<endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double elapsed =std::chrono::duration<double>(end - t0).count();
    std::cout << elapsed << " seconds\n";
    //std::this_thread::sleep_for(std::chrono::seconds(300)); // 暂停3秒
    return 0;

}
