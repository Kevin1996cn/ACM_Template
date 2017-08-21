#include<vector>
#include<map>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<algorithm>
using namespace std;

typedef long long LL;
const LL MAXN=10010;


// 乘积取模
// Verified!
LL mulMod(LL a,LL b,LL mod)
{
    LL ret=0;
    while(b)
    {
        ret=(ret+(b%2*a)) % mod;
        a=a*2%mod;
        b=b/2;
    }
    return ret%mod;
}



// 快速幂
// 求x^n%mod
// Verified!
LL powMod(LL x,LL n,LL mod)
{
    LL res=1;
    while(n>0)
    {
        if(n&1) res=res*x % mod;
        x=x*x % mod;
        n>>=1;
    }
    return res;
}



// 素数筛法
// prime[0]存储素数个数
// 复杂度为O(n)
// Verified
LL prime[MAXN+1];
void getPrime()
{
    memset(prime,0,sizeof(prime));
    for(LL i=2;i<=MAXN;i++)
    {
        if(!prime[i])prime[++prime[0]]=i;
        for(LL j=1;j<=prime[0]&&prime[j]<=MAXN/i;j++)
        {
            prime[prime[j]*i]=1;
            if(i%prime[j]==0) break;
        }
    }
}



// 求最大公约数
// Verified!
LL gcd(LL a,LL b)
{
    if(b==0) return a;
    return gcd(b,a%b);
}


// 拓展欧几里得算法
// 求ax+by=gcd(a,b)的一组解，满足x和y的绝对值之和最小
// 其他解为x=x0+kb',y=y0-ka'
// a'=a/gcd(a,b),b'=b/gcd(a,b)
// Verified!
LL extgcd(LL a,LL b,LL &x,LL &y)
{
    LL d=a;
    if(b!=0)
    {
        d=extgcd(b,a%b,y,x);
        y-=(a/b)*x;
    }
    else { x=1; y=0; }
    return d;
}


// 欧拉函数
// 求不超过n且和n互质的数的个数
// 对于和n互质的x，有(x^(phi(n)) % n == 1
// Verified!
LL eularPhi(LL n)
{
    LL res=n;
    for(LL i=2;i*i<=n;i++)
    {
        if(n%i==0)
        {
            res=res/i*(i-1);
            for(;n%i==0;n/=i);
        }
    }
    if(n!=1) res=res/n*(n-1);
    return res;
}


// O(nloglogn)时间筛出欧拉函数表
// Verified!
void phiTable(LL n,LL phi[])
{
	for(LL i=2;i<=n;i++) phi[i]=0;
	phi[1]=1;
	for(LL i=2;i<=n;i++) if(!phi[i])
		for(LL j=i;j<=n;j+=i)
		{
			if(!phi[j]) phi[j]=j;
			phi[j]=phi[j]/i*(i-1);
		}
}


// 求逆元
// a和ｍ应该互质
// Verified!
LL modInverse(LL a,LL m)
{
    LL x,y;
    extgcd(a,m,x,y);
    return (m+x%m)%m;
}


// 求逆元
// a和ｍ应该互质
// 根据欧拉定理：a的逆即a^(phi(m)-1)
// Verified!
LL inv(LL a,LL m)
{
    return powMod(a,m-2,m);
	// return powMod(a,eularPhi(m)-1,m);
}


// 求解同余方程组
// 返回一个<b,m>的数对
pair<int,int> linearCongruence(const vector<int> &A,const vector<int> &B,
                               const vector<int> &M)
{
    int x=0,m=1;
    for(unsigned int i=0;i<A.size();i++)
    {
        int a=A[i]*m,b=B[i]-A[i]*x,d=gcd(M[i],a);
        if(b%d!=0) return make_pair(0,-1); //无解
        int t=b/d*modInverse(a/d,M[i]/d) % (M[i]/d);
        x=x+m*t;
        m*=M[i]/d;
    }
    return make_pair(x%m,m);
}


// 离散对数
// 求解a^x === b (mod n)
// n是不是质数都可以，无解返回-1
// 时间复杂度为O(n^(1/2)logn)
#define MOD2 76543
int hs[MOD2],head[MOD2],nex[MOD2],id[MOD2],top;
void insert(int x,int y)
{
    int k = x%MOD2;
    hs[top] = x, id[top] = y, nex[top] = head[k], head[k] = top++;
}
int find(int x)
{
    int k = x%MOD2;
    for(int i = head[k]; i != -1; i = nex[i])
        if(hs[i] == x)
            return id[i];
    return -1;
}
int BSGS(int a,int b,int n)
{
    memset(head,-1,sizeof(head));
    top = 1;
    if(b == 1)return 0;
    int m = sqrt(n*1.0), j;
    long long x = 1, p = 1;
    for(int i = 0; i < m; ++i, p = p*a%n)insert(p*b%n,i);
    for(long long i = m; ;i += m)
    {
        if( (j = find(x = x*p%n)) != -1 )return i-j;
        if(i > n)break;
    }
    return -1;
}



// Lucas定理
// 求C(n,m)%p
// O(lognlogn)
int fact[MAXN];//需预处理ｎ! mod p的表
LL Lucas(LL n,LL m,LL p)  
{  
    LL ans = 1;  
    while(n&&m)  
    {  
        LL a = n%p;  
        LL b = m%p;  
        if(a < b)return 0;  
        ans = ans*fact[a]%p*inv(fact[b]*fact[a-b]%p,p)%p;  
        n /= p;  
        m /= p;  
    }  
    return ans;  
}  


// MillerRabin算法
// 随机化检验素数
// 原理是费马小定理和二次检验定理，一次检验正确率为75%
const int checkTime=8; //随机化测试次数
bool check(LL a,LL n,LL x,LL t)
{
    LL ret=powMod(a,x,n);
    LL last=ret;
    for(LL i=1;i<=t;i++)
    {
        ret=mulMod(ret,ret,n);
		// 二次检验定理:如果p是一个素数，0<x<p,则方程x^2≡1(mod p)的解为x=1或x=p-1。
        if(ret==1 && last!=1 && last!=n-1) return true;
        last=ret;
    }
    if(ret!=1) return true;
    else return false;
}
bool millerRabin(LL n)
{
    if(n<2) return false;
    if(n==2) return true;
    if((n&1)==0) return false;
    LL x=n-1;
    LL t=0;
    while((x&1)==0) {x>>=1;t++;}
    srand(time(NULL));
    for(LL i=0;i<checkTime;i++)
    {
         LL a=rand()%(n-1)+1;
         if(check(a,n,x,t))
             return false;
    }
    return true;
}


// 约瑟夫环
// n个人每次第m个人出列,第k次出列的人的编号
int Josephus(int n,int m,int k)
{
    if(k==1) return m%n;
    return (Josephus(n-1,m,k-1)+m) % n;
}

// 递推写法
// 最后剩下的人的编号
int Josephus(int n,int m)
{
    int ans=0;
    for(int i=2;i<=n;i++)
        ans=(ans+m)%i;
    return ans;
}


// 高斯消元
// 求解Ax=b
// 当无解或无穷多解，返回一个长度为0的数组
const double EPS = 1E-8;
typedef vector<double> vec;
typedef vector<vec> mat;

vec gaussJordan(const mat& A,const vec& b)
{
    int n=A.size();
    mat B(n,vec(n+1));
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            B[i][j]=A[i][j];
    for(int i=0;i<n;i++) B[i][n]=b[i];
    for(int i=0;i<n;i++)
    {
		// 找到当前列最大的那行
        int pivot=i;
        for(int j=i;j<n;j++)
            if(abs(B[j][i]) > abs(B[pivot][i]))
                pivot=j;
        swap(B[i],B[pivot]);

        // 无穷解或无解
        if(abs(B[i][i])<EPS) return vec();

		// 消元
        for(int j=i+1;j<=n;j++) B[i][j]/=B[i][i];
        for(int j=0;j<n;j++)
            if(i!=j)
                for(int k=i+1;k<=n;k++)
                    B[j][k]-=B[j][i]*B[i][k];
    }
	vec x(n);
	for(int i=0;i<n;i++) x[i]=B[i][n];
	return x;
}

// 矩阵乘法/快速幂
// 矩阵r,c需要自己读入的时候设置
// Verified!
const int MOD=1000000000+7;
const int MAXR=100;
const int MAXC=100;
struct Matrix
{
    int m[MAXR][MAXC];
    int r,c;
    Matrix()
    {
        r=0,c=0;
        memset(m,0,sizeof(m));
    }
};
//若不能相乘，返回空矩阵
Matrix operator * (const Matrix &a,const Matrix &b)
{
    Matrix ans;
    int r=a.r,c=b.c,x=0;
    if(a.c!=b.r) return Matrix();
    x=a.c;
    ans.r=a.r;
    ans.c=b.c;
    for (int i=1; i<=r; i++)
        for (int j=1; j<=c; j++)
        {
            ans.m[i][j]=0;
            for (int k=1; k<=x; k++)
                ans.m[i][j]=(ans.m[i][j]+a.m[i][k]*b.m[k][j]) % MOD;
        }
    return ans;
}
// 若不能乘方，返回空矩阵
Matrix operator ^ (Matrix &base,int pow)
{
    if(base.r!=base.c) return Matrix();
    int x=base.r;
    Matrix ans;
    ans.r=ans.c=x;
    for (int i=1; i<=x; i++) ans.m[i][i]=1;
    while (pow)
    {
        if (pow&1) ans=ans*base;
        base=base*base;
        pow>>=1;
    }
    return ans;
}

// 求组合数
// Verified
int C[MAXN][MAXN];
void initC()
{
	memset(C,0,sizeof(C));
	for(int i=0;i<MAXN;i++)
		C[i][0]=1;
	for(int i=1;i<MAXN;i++)
		for(int j=1;j<MAXN;j++)
			C[i][j]=(C[i-1][j]+C[i-1][j-1])%MOD;
}


// SG函数
// Verified!
int sg[MAXN];  
bool vis[MAXN];  
void sg_solve(LL *s,int t,int N)   //N求解范围 S[]数组是可以每次取的值，t是s的长度。  
{  
    int i,j;  
    memset(sg,0,sizeof(sg));  
    for(i=1;i<=N;i++)  
    {  
        memset(vis,0,sizeof(vis));  
        for(j=0;j<t;j++)  
            if(i - s[j] >= 0)  
                vis[sg[i-s[j]]] = 1;  
        for(j=0;j<=N;j++)  
            if(!vis[j])  
                break;  
        sg[i] = j;  
    }  
}


// 辛普森法求积分
// 占空位，自由改动
double F(double x)
{
	return x;
}
// 三点Simpson法，这里要求F是一个全局函数。
double simpson(double a,double b)
{
	double c=a+(b-a)/2;
	return (F(a)+4*F(c)+F(b))*(b-a)/6;
}
// 自适应Simpson公式（递归过程）。已知整个区间[a,b]上的三点Simpson值A
double asr(double a,double b,double eps,double A)
{
	double c=a+(b-a)/2;
	double L=simpson(a,c),R=simpson(c,b);
	if(fabs(L+R-A) <= 15*eps) return L+R+(L+R-A)/15.0;
	return asr(a,c,eps/2,L) + asr(c,b,eps/2,R);
}
// 自适应Simpson公式（主过程）
double asr(double a,double b,double eps)
{
	return asr(a,b,eps,simpson(a,b));
}



// FFT
// 时间复杂度是O(nlogn)
// Verified!
const double PI = acos(-1.0);
//复数结构体
struct complex
{
    double r,i;
    complex(double _r = 0.0,double _i = 0.0)
    {
        r = _r; i = _i;
    }
    complex operator +(const complex &b)
    {
        return complex(r+b.r,i+b.i);
    }
    complex operator -(const complex &b)
    {
        return complex(r-b.r,i-b.i);
    }
    complex operator *(const complex &b)
    {
        return complex(r*b.r-i*b.i,r*b.i+i*b.r);
    }
};
/*
 * 进行FFT和IFFT前的反转变换。
 * 位置i和 （i二进制反转后位置）互换
 * len必须为2的幂
 * 雷德算法，这么做的意义是进行蝶形运算，能减少大概一半的运算量。
 */
void change(complex y[],int len)
{
    int i,j,k;
    for(i = 1, j = len/2;i < len-1; i++)
    {
        if(i < j)swap(y[i],y[j]);
        //交换互为小标反转的元素，i<j保证交换一次
        //i做正常的+1，j左反转类型的+1,始终保持i和j是反转的
        k = len/2;
        while( j >= k)
        {
            j -= k;
            k /= 2;
        }
        if(j < k) j += k;
    }
}
/*
 * 做FFT
 * len必须为2^k形式，
 * on==1时是DFT，on==-1时是IDFT
 */
void fft(complex y[],int len,int on)
{
    change(y,len);
    for(int h = 2; h <= len; h <<= 1)
    {
        complex wn(cos(-on*2*PI/h),sin(-on*2*PI/h));
        for(int j = 0;j < len;j+=h)
        {
            complex w(1,0);
            for(int k = j;k < j+h/2;k++)
            {
                complex u = y[k];
                complex t = w*y[k+h/2];
                y[k] = u+t;
                y[k+h/2] = u-t;
                w = w*wn;
            }
        }
    }
    if(on == -1)
        for(int i = 0;i < len;i++)
            y[i].r /= len;
}
complex x1[MAXN],x2[MAXN];
char str1[MAXN/2],str2[MAXN/2];
int sum[MAXN];
void work()
{
    while(scanf("%s%s",str1,str2)==2)
    {
        int len1 = strlen(str1);
        int len2 = strlen(str2);
        int len = 1;
        while(len < len1*2 || len < len2*2)len<<=1;
        for(int i = 0;i < len1;i++)
            x1[i] = complex(str1[len1-1-i]-'0',0);
        for(int i = len1;i < len;i++)
            x1[i] = complex(0,0);
        for(int i = 0;i < len2;i++)
            x2[i] = complex(str2[len2-1-i]-'0',0);
        for(int i = len2;i < len;i++)
            x2[i] = complex(0,0);
        //求DFT
        fft(x1,len,1);
        fft(x2,len,1);
        for(int i = 0;i < len;i++)
            x1[i] = x1[i]*x2[i];
        fft(x1,len,-1);
        for(int i = 0;i < len;i++)
            sum[i] = (int)(x1[i].r+0.5);
        for(int i = 0;i < len;i++)
        {
            sum[i+1]+=sum[i]/10;
            sum[i]%=10;
        }
        len = len1+len2-1;
        while(sum[len] <= 0 && len > 0)len--;
        for(int i = len;i >= 0;i--)
            printf("%c",sum[i]+'0');
        printf("\n");
    }
}


// NTT
// O(nlogn)
// Verified!
// const LL MAXN = 1 << 18;
const LL P = (479 << 21) + 1;  
// P-1要是N的倍数
// const LL P=998244353; 2^23 * 119
// const LL P=1004535809; 2^21 * 479
const LL G = 3;  
const LL NUM = 20;  
  
LL  wn[NUM];  
LL  a[MAXN], b[MAXN];  
char A[MAXN], B[MAXN];  
  
LL quick_mod(LL a, LL b, LL m)  
{  
    LL ans = 1;  
    a %= m;  
    while(b)  
    {  
        if(b & 1)  
        {  
            ans = ans * a % m;  
            b--;  
        }  
        b >>= 1;  
        a = a * a % m;  
    }  
    return ans;  
}  
  
void GetWn()  
{  
    for(LL i=0; i<NUM; i++)  
    {  
        LL t = 1 << i;  
        wn[i] = quick_mod(G, (P - 1) / t, P);  
    }  
}  
  
void Prepare(char A[], char B[], LL a[], LL b[], LL &len)  
{  
    len = 1;  
    LL len_A = strlen(A);  
    LL len_B = strlen(B);  
    while(len <= 2 * len_A || len <= 2 * len_B) len <<= 1;  
    for(LL i=0; i<len_A; i++)  
        A[len - 1 - i] = A[len_A - 1 - i];  
    for(LL i=0; i<len - len_A; i++)  
        A[i] = '0';  
    for(LL i=0; i<len_B; i++)  
        B[len - 1 - i] = B[len_B - 1 - i];  
    for(LL i=0; i<len - len_B; i++)  
        B[i] = '0';  
    for(LL i=0; i<len; i++)  
        a[len - 1 - i] = A[i] - '0';  
    for(LL i=0; i<len; i++)  
        b[len - 1 - i] = B[i] - '0';  
}  
  
void Rader(LL a[], LL len)  
{  
    LL j = len >> 1;  
    for(LL i=1; i<len-1; i++)  
    {  
        if(i < j) swap(a[i], a[j]);  
        LL k = len >> 1;  
        while(j >= k)  
        {  
            j -= k;  
            k >>= 1;  
        }  
        if(j < k) j += k;  
    }  
}  
  
void NTT(LL a[], LL len, LL on)  
{  
    Rader(a, len);  
    LL id = 0;  
    for(LL h = 2; h <= len; h <<= 1)  
    {  
        id++;  
        for(LL j = 0; j < len; j += h)  
        {  
            LL w = 1;  
            for(LL k = j; k < j + h / 2; k++)  
            {  
                LL u = a[k] % P;  
                LL t = w * (a[k + h / 2] % P) % P;  
                a[k] = (u + t) % P;  
                a[k + h / 2] = ((u - t) % P + P) % P;  
                w = w * wn[id] % P;  
            }  
        }  
    }  
    if(on == -1)  
    {  
        for(LL i = 1; i < len / 2; i++)  
            swap(a[i], a[len - i]);  
        LL Inv = quick_mod(len, P - 2, P);  
        for(LL i = 0; i < len; i++)  
            a[i] = a[i] % P * Inv % P;  
    }  
}  
  
void Conv(LL a[], LL b[], LL n)  
{  
    NTT(a, n, 1);  
    NTT(b, n, 1);  
    for(LL i = 0; i < n; i++)  
        a[i] = a[i] * b[i] % P;  
    NTT(a, n, -1);  
}  
  
void Transfer(LL a[], LL n)  
{  
    LL t = 0;  
    for(LL i = 0; i < n; i++)  
    {  
        a[i] += t;  
        if(a[i] > 9)  
        {  
            t = a[i] / 10;  
            a[i] %= 10;  
        }  
        else t = 0;  
    }  
}  
  
void Print(LL a[], LL n)  
{  
    bool flag = 1;  
    for(LL i = n - 1; i >= 0; i--)  
    {  
        if(a[i] != 0 && flag)  
        {  
            printf("%lld", a[i]);  
            flag = 0;  
        }  
        else if(!flag)  
            printf("%lld", a[i]);  
    }  
	if(flag) printf("0");
	printf("\n");
}  
  
LL work2()  
{  
    GetWn();  
    while(scanf("%s%s", A, B)!=EOF)  
    {  
        LL len;  
        Prepare(A, B, a, b, len);  
        Conv(a, b, len);  
        Transfer(a, len);  
        Print(a, len);  
    }  
    return 0;  
}  

// FWT
// 求位运算卷积
// FWT(C)=FWT(A)*FWT(B)
// 取模直接在FWT内部取模即可
void FWT_And(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  FWT_And(A, len2);
  FWT_And(A + len2, len2);
  for (int i = 0; i < len2; ++i)
    A[i] += A[i + len2];
}
void IFWT_And(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  for (int i = 0; i < len2; ++i)
    A[i] -= A[i + len2];
  IFWT_And(A, len2);
  IFWT_And(A + len2, len2);
}

void FWT_Or(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  FWT_Or(A, len2);
  FWT_Or(A + len2, len2);
  for (int i = 0; i < len2; ++i)
    A[i + len2] += A[i];
}
void IFWT_Or(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  for (int i = 0; i < len2; ++i)
    A[i + len2] -= A[i];
  IFWT_Or(A, len2);
  IFWT_Or(A + len2, len2);
}

// Verified!
void FWT_Xor(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  FWT_Xor(A, len2);
  FWT_Xor(A + len2, len2);
  for (int i = 0; i < len2; ++i) {
    int x = A[i], y = A[i + len2];
    A[i] = x + y;
    A[i + len2] = x - y;
  }
}
void IFWT_Xor(int *A, int len) {
  if (len == 1) return;
  int len2 = len >> 1;
  for (int i = 0; i < len2; ++i) {
    int x = A[i], y = A[i + len2];
    A[i] = (x + y) >> 1;
    A[i + len2] = (x - y) >> 1;
  }
  IFWT_Xor(A, len2);
  IFWT_Xor(A + len2, len2);
}


// 母函数
// 对于表达式(1+x+x^2)(x^8+x^10)(x^5+x^10+x^15+x^20)，v[3]={1,2,5}，n1[3]={0,4,1}，n2[3]={2,5,4}。
// K为表达式数量，P为要求的最大次数
// a为结果，b为中间变量
// 其实没什么大用，用dp就能解决的事，我也是闲的￣へ￣
// Verified!
void generateFunc(int a[],int b[],int n1[],int n2[],int v[],int K,int P)
{
	memset(a,0,sizeof(int)*P);
	a[0]=1;
	int last=0;
	for (int i=0;i<K;i++)
	{
		int last2=min(last+n2[i]*v[i],P);//计算下一次的last
		memset(b,0,sizeof(int)*(last2+1));//只清空b[0..last2]
		for (int j=n1[i];j<=n2[i]&&j*v[i]<=last2;j++)//这里是last2
			for (int k=0;k<=last&&k+j*v[i]<=last2;k++)//这里一个是last，一个是last2
				b[k+j*v[i]]+=a[k];
		memcpy(a,b,sizeof(int)*(last2+1));//b赋值给a，只赋值0..last2
		last=last2;//更新last
	}
}


// 蔡勒公式
// 注意判断只有闰年才有2.29号
// 只适用于1582.10.4之后的计算
// Verified!
int caile(int year,int month,int day)//根据日期判断出星期几
{
    if(month==1||month==2)
    {
        month+=12;
        year--;
    }
    int c=year/100;
    int y=year%100;
    int m=month;
    int d=day;
    int W=c/4-2*c+y+y/4+26*(m+1)/10+d-1;
    W=(W%7+7)%7;
    if(!W) W=7;
    return W;
}



// 求解x^2=n(mod p)
struct T {
  LL p, d;
};

LL w;

T multi_er(T a, T b, LL m) {
  T ans;
  ans.p = (a.p * b.p % m + a.d * b.d % m * w % m) % m;
  ans.d = (a.p * b.d % m + a.d * b.p % m) % m;
  return ans;
}

T power(T a, LL b, LL m) {
  T ans;
  ans.p = 1;
  ans.d = 0;
  while(b) {
    if(b & 1) {
      ans = multi_er(ans, a, m);
      b--;
    }
    b >>= 1;
    a = multi_er(a, a, m);
  }
  return ans;
}

LL Legendre(LL a, LL p) {
  return quick_mod(a, (p-1)>>1, p);
}

LL mod(LL a, LL m) {
  a %= m;
  if(a < 0) a += m;
  return a;
}

LL Solve(LL n,LL p) {
  if(p == 2) return 1;
  if (Legendre(n, p) + 1 == p)
    return -1;
  LL a = -1, t;
  while(true) {
    a = rand() % p;
    t = a * a - n;
    w = mod(t, p);
    if(Legendre(w, p) + 1 == p) break;
  }
  T tmp;
  tmp.p = a;
  tmp.d = 1;
  T ans = power(tmp, (p + 1)>>1, p);
  return ans.p;
}



// 组合数性质
// C[n][m]=C[n-1][m]+C[n-1][m-1];
// C[n][m]=C[n][m-1]*(n-m+1)/m;
// m*C[n][m]=n*C[n-1][m]

// 费马小定理
// 假如p是质数，且gcd(a,p)=1，那么 a^(p-1)≡1（mod p）


// 威尔逊定理
// 当且仅当p为素数时：( p -1 )! ≡ -1 ( mod p )


// Burnside引理
// 对一个置换，若一个着色方案s经过置换后不变，称s为f的不动点。
// 将f的不动点个数记为C(f),则可证明等价类数目为所有C(f)的平均值。


// Polya定理
// 若置换f分解为m(f)个循环的乘积，那么每个循环内所有格子颜色必然相同。
// 假设涂k种颜色，则C(f)=k^m(f)。
// 等价类的个数等于所有置换f的k^m(f)的平均数。


// 莫比乌斯反演
// f(n)=sum(g(d)) (d|n) <==> g(n)=sum(u(n/d)f(d))  (d|n)
// f(n)=sum(g(d)) (n|d) <==> g(n)=sum(u(d/n)f(d))  (n|d)  （第二种形式）
// 若n可以被1以外完全平方数整除，u(n)=0
// 否则设n的质因数个数为k，u(n)=(-1)^k
// 考虑求长度为n的没有周期性的字符串个数，即相当于求周期是n的字符串的个数。
// f(d)表示周期是d的约数字符串的个数，g(d)表示周期恰好为d的字符串的个数。


// 威佐夫博弈。
// 威佐夫博弈简单的说，就是有两堆石子，双方每次可以取一堆石子中的任意个，不能不取
// 或者取两堆石子中的相同个。先取完者赢。
// 然后发现，对于（0，0）（1，2）（3，5）（4，7）（6，13）（8，13）…均为先手必输。
// 发现规律，差值递增，而每个数对的第一个数都是差值*1.618。


// 第一类Stirling数 s(p,k)
// s(p,k)的一个的组合学解释是：将p个物体排成k个非空循环排列的方法数。
// s(p,k)的递推公式： s(p,k)=(p-1)*s(p-1,k)+s(p-1,k-1) ,1<=k<=p-1
// 边界条件：s(p,0)=0 ,p>=1  s(p,p)=1  ,p>=0
 

// 第二类Stirling数 S(p,k)
// S(p,k)的一个组合学解释是：将p个物体划分成k个非空的不可辨别的（可以理解为盒子没有编号）集合的方法数。
// k!S(p,k)是把p个人分进k间有差别(如：被标有房号）的房间(无空房）的方法数。
// S(p,k)的递推公式是：S(p,k)=k*S(p-1,k)+S(p-1,k-1) ,1<= k<=p-1
// 边界条件：S(p,p)=1 ,p>=0    S(p,0)=0 ,p>=1
