// 找到最长的字符串，判断其他字符串是否都在这个字符串里。
// find返回的是有多少单词模板在查询的字符串中
// Verified!
#include<cstdio>
#include<cstdlib>
#include<vector>
#include<set>
#include<cstring>
#include<algorithm>
#include<cmath>
#include<iostream>
#include<queue>
using namespace std;
typedef long long LL;
const LL MOD=998244353;
const LL MAXN=1e5+10;

struct AC
{
    int ch[MAXN][26];
    int val[MAXN];
    int sz;
    AC() { init(); }
    void init() { sz=1; memset(ch[0],0,sizeof(ch[0])); val[0]=0;}
    int idx(char c) { return c-'a'; }
    void insert(char *s)
    {
        int u=0,n=strlen(s);
        for(int i=0;i<n;i++)
        {
            int c=idx(s[i]);
            if(!ch[u][c])
            {
                memset(ch[sz],0,sizeof(ch[sz]));
                val[sz]=0;
                ch[u][c]=sz++;
            }
            u=ch[u][c];
        }
        val[u]++;
    }
    int f[MAXN];
    int last[MAXN];
    int find(char *T)
    {
        int ans=0;
        int n=strlen(T);
        int j=0;
        for(int i=0;i<n;i++)
        {
            int c=idx(T[i]);
            j=ch[j][c];
            if(val[j])
            {
                ans+=val[j];
                val[j]=0;
            }
            int k=j;
            while(last[k])
            {
                k=last[k];
                ans+=val[k];
                val[k]=0;
            }
        }
        return ans;
    }
    void getFail()
    {
        queue<int> q;
        f[0]=0;
        for(int c=0;c<26;c++)
        {
            int u=ch[0][c];
            if(u) {f[u]=0;q.push(u);last[u]=0;}
        }
        while(!q.empty())
        {
            int r=q.front();q.pop();
            for(int c=0;c<26;c++)
            {
                int u=ch[r][c];
                if(!u) 
                {
                    ch[r][c]=ch[f[r]][c];
                    continue;
                }
                q.push(u);
                int v=f[r];
                while(v && !ch[v][c]) v=f[v];
                f[u]=ch[v][c];
                last[u]=val[f[u]]?f[u]:last[f[u]];
            }
        }
    }
};
AC T;
char s[MAXN];
char ms[MAXN];
int main()
{
#ifdef LOCAL
    freopen("in.txt","r",stdin);
    // freopen("out.txt","w",stdout);
#endif
    int t;
    scanf("%d",&t);
    for(int tt=1;tt<=t;tt++)
    {
        T.init();
        int msl=0;
        int n;
        scanf("%d",&n);
        for(int i=1;i<=n;i++)
        {
            scanf("%s",s);
            int l=strlen(s);
            if(l>msl)
            {
                if(msl) T.insert(ms);
                msl=l;
                strcpy(ms,s);
            }
            else T.insert(s);
        }
        T.getFail();
        int tmp=T.find(ms);
        if(tmp==n-1) printf("%s\n",ms);
        else printf("No\n");
    }
    return 0;

}
