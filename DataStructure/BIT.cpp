//树状数组
//可用于快速查询第k大的值
//别忘了清零tree[]
//Verified!
#define MAXN 100010

int n,tree[MAXN];
int lowbit(int x)
{
    return x&(-x);
}
void add(int k,int num)
{
    while(k<=n)
    {
        tree[k]+=num;
        k+=lowbit(k);
    }
}
int sum(int k)
{
    int ans=0;
    while(k)
    {
        ans+=tree[k];
        k-=lowbit(k);
    }
    return ans;
}
