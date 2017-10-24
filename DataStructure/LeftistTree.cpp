// 左偏树，可合并堆
// Verified!
// O(logn)时间合并
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<algorithm>
#include<vector>
#include<set>
#include<stack>
#include<functional>
using namespace std;
typedef long long LL;
const LL MAXN=1e5+10;
const LL MAXM=2e5+10;
const LL INF=1e9+9;

typedef pair<int,int> query;
vector<query> vec;
int LH[MAXN],RH[MAXN],L[MAXN],R[MAXN];
int X[MAXN],O[MAXN];

// 左偏树相关
int tot, v[MAXM], l[MAXM], r[MAXM], d[MAXM], Heap[MAXN];

// 合并左偏树
int merge(int x, int y) {
    if(x == 0) {
        return y;
    }
    if(y == 0) {
        return x;
    }
    if(v[x] > v[y]) {
        swap(x, y);
    }
    r[x] = merge(r[x], y);
    if(d[l[x]] < d[r[x]]) {
        swap(l[x], r[x]);   
    }
    d[x] = d[r[x]] + 1;
    return x;
}

// 初始化可并堆结点
inline int init(int x) {
    v[++tot] = x;
    l[tot] = r[tot] = d[tot] = 0;
    return tot;
}

// 左偏树的插入操作
inline int insert(int x, int y) {
    return merge(x, init(y));
}

// 取得左偏树中的最小值
inline int top(int x) {
    return v[x];
}

// 弹出左偏树
inline int pop(int x) {
    return merge(l[x], r[x]);
}

// 判断左偏树是否非空
inline bool empty(int x) {
    return x == 0;
}

//  初始化可并堆
void initHeap() {
    memset(Heap, 0, sizeof(Heap));
    tot = 0;
}

// 并查集相关
int p[MAXN];

// 初始化并查集
void initSet(int n) {
    for(int i = 1; i <= n; i++) {
        p[i] = i;
    }
}

// 查找集合的祖先
int find(int x) {
    return x == p[x] ? x : p[x] = find(p[x]);
}

// 合并集合
inline void Union(int x, int y) {
    x = find(x);
    y = find(y);
    if(x == y) {
        return;
    }
    p[y] = x;
    if(x < y) {
        RH[x] = RH[y];
        R[x] = R[y];
        L[R[x]] = x;
    }
    else {
        LH[x] = LH[y];
        L[x] = L[y];
        R[L[x]] = x;
    }
    // 合并可并堆
    Heap[x] = merge(Heap[x], Heap[y]);
    X[x] += X[y];
    O[x] += O[y];
}
int main()
{
#ifdef LOCAL
    freopen("in.txt","r",stdin);
#endif
    int t;
    scanf("%d",&t);
    for(int tt=1;tt<=t;tt++)
    {
        int n,m;
        scanf("%d%d",&n,&m);
        LH[1]=RH[n]=INF;
        L[n]=n-1;
        for(int i=1;i<=n-1;i++) 
        {
            scanf("%d",&RH[i]);
            LH[i+1]=RH[i];
            L[i]=i-1;
            R[i]=i+1;
        }
        initHeap();
        vec.clear();
        int ans=0;
        for(int i=1;i<=m;i++)
        {
            int x,y,z;
            scanf("%d%d%d",&x,&y,&z);
            if(z==1) vec.push_back(query(y+1,x));
            else 
            {
                Heap[x]=Heap[x]?insert(Heap[x],y):init(y);
                ans++;
            }
        }
        initSet(n);
        sort(vec.begin(),vec.end());
        for(int i=1;i<=n;i++) X[i]=O[i]=0;
        for(unsigned int i=0;i<vec.size();i++)
        {
            int x,y;
            x=find(vec[i].second);
            y=vec[i].first;
            while(y>LH[x])
            {
                Union(x,L[x]);
                x=find(x);
            }
            while(y>RH[x])
            {
                Union(x,R[x]);
                x=find(x);
            }
            while(!empty(Heap[x]) && top(Heap[x])<y)
            {
                Heap[x]=pop(Heap[x]);
                X[x]++;
            }
            if(++O[x]>=X[x])
            {
                ans+=(O[x]-X[x]);
                O[x]=X[x]=0;
            }
        }
        printf("Case #%d: %d\n", tt, ans);
    }
    return 0;
}
