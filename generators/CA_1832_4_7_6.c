#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#define K 7
#define V 6
#define NS 35
#define NT (NS*1296)
#define MAXN 2000
static volatile sig_atomic_t stop=0;
static void onterm(int s){(void)s;stop=1;}
static int sc[NS][4], csets[K][20], ncs[K];
static unsigned char A[MAXN][K],B[MAXN][K],S[MAXN][K],IN[MAXN][K];
static int N,BN,SN,INN;
static unsigned short cnt[NT];static int unc[NT],pos[NT],nunc;
static unsigned char vt[NT];static int vic=-1,nv1=0,W=2;
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int idx(const unsigned char*r,int s){return s*1296+r[sc[s][0]]*216+r[sc[s][1]]*36+r[sc[s][2]]*6+r[sc[s][3]];}
static inline void inc(int t){unsigned c=cnt[t]++;if(c==0){int p=pos[t];int l=unc[--nunc];unc[p]=l;pos[l]=p;pos[t]=-1;}else if(c==1&&vt[t])nv1--;}
static inline void dec(int t){unsigned c=--cnt[t];if(c==0){pos[t]=nunc;unc[nunc++]=t;}else if(c==1&&vt[t])nv1++;}
static void build(){memset(cnt,0,sizeof cnt);memset(vt,0,sizeof vt);vic=-1;nv1=0;nunc=0;
 for(int t=0;t<NT;t++){pos[t]=nunc;unc[nunc++]=t;}
 for(int r=0;r<N;r++)for(int s=0;s<NS;s++)inc(idx(A[r],s));}
static void clearvic(){if(vic>=0)for(int s=0;s<NS;s++)vt[idx(A[vic],s)]=0;vic=-1;nv1=0;}
static void setvic(int r){clearvic();vic=r;nv1=0;for(int s=0;s<NS;s++){int t=idx(A[r],s);vt[t]=1;if(cnt[t]==1)nv1++;}}
static int rowloss(int r){int l=0;for(int s=0;s<NS;s++)if(cnt[idx(A[r],s)]==1)l++;return l;}
static void removerow(int r){for(int s=0;s<NS;s++)dec(idx(A[r],s));memcpy(A[r],A[N-1],K);N--;}
static int pickvic(){static int L[MAXN];int mn=1<<30;for(int r=0;r<N;r++){L[r]=rowloss(r);if(L[r]<mn)mn=L[r];}
 int slack=(rnd()%4==0)?2:(int)(rnd()%2);int br=0,nb=0;
 for(int r=0;r<N;r++)if(L[r]<=mn+slack){nb++;if(rnd()%nb==0)br=r;}
 fprintf(stderr,"victim %d loss %d (min %d)\n",br,L[br],mn);return br;}
static int verify(unsigned char (*X)[K],int n){static unsigned char cv[NT];memset(cv,0,sizeof cv);
 for(int r=0;r<n;r++)for(int s=0;s<NS;s++)cv[idx(X[r],s)]=1;for(int t=0;t<NT;t++)if(!cv[t])return 0;return 1;}
static void setcell(int r,int c,int nv){for(int j=0;j<20;j++)dec(idx(A[r],csets[c][j]));A[r][c]=nv;
 for(int j=0;j<20;j++)inc(idx(A[r],csets[c][j]));}
static inline int delta(int r,int c,int nv){unsigned char row[K];memcpy(row,A[r],K);int d=0;
 for(int j=0;j<20;j++){int s=csets[c][j];int t=idx(row,s);if(cnt[t]==1)d+=W;else if(cnt[t]==2&&vt[t])d++;}
 row[c]=nv;for(int j=0;j<20;j++){int s=csets[c][j];int t=idx(row,s);if(cnt[t]==0)d-=W;else if(cnt[t]==1&&vt[t])d--;}
 return d;}
static long long walk(long long steps){long long acc=0;
 for(long long i=0;i<steps;i++){int r=rnd()%N,c=rnd()%K,nv=rnd()%V;if(nv==A[r][c])continue;
  int ok=1;for(int j=0;j<20;j++)if(cnt[idx(A[r],csets[c][j])]==1){ok=0;break;}
  if(ok){setcell(r,c,nv);acc++;}}
 return acc;}
static int dropfree(){int d=0;for(int r=0;r<N;){if(rowloss(r)==0){removerow(r);d++;}else r++;}return d;}
static void out(unsigned char (*X)[K],int n){for(int r=0;r<n;r++){for(int k=0;k<K;k++)printf(k?" %d":"%d",X[r][k]);printf("\n");}}
int main(){
 signal(SIGTERM,onterm);
 int s=0;for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){
  sc[s][0]=a;sc[s][1]=b;sc[s][2]=c;sc[s][3]=d;csets[a][ncs[a]++]=s;csets[b][ncs[b]++]=s;csets[c][ncs[c]++]=s;csets[d][ncs[d]++]=s;s++;}
 double T0=now(),TL=3420;rs^=(unsigned long long)time(NULL)*2654435761ULL;
 char*f=getenv("BEST_CANDIDATE");INN=0;
 if(f){FILE*F=fopen(f,"r");if(F){int x[K];while(INN<MAXN){int ok=1;for(int k=0;k<K;k++)if(fscanf(F,"%d",&x[k])!=1){ok=0;break;}if(!ok)break;for(int k=0;k<K;k++)IN[INN][k]=((x[k]%V)+V)%V;INN++;}fclose(F);}}
 if(INN==0){INN=1900;for(int r=0;r<INN;r++)for(int k=0;k<K;k++)IN[r][k]=rnd()%V;}
 int valid=verify(IN,INN);
 memcpy(A,IN,(size_t)INN*K);N=INN;build();memcpy(B,IN,(size_t)INN*K);BN=INN;memcpy(S,A,(size_t)N*K);SN=N;
 fprintf(stderr,"input N=%d valid=%d unc=%d\n",INN,valid,nunc);
 setvic(pickvic());
 enum{OFF=256};static double ex[OFF];
 const long long L=1LL<<21;const int ST=3;const double Thi=0.7,Tlo=0.1;
 long long it=0,cyc=0,lastImp=0;int bestE=1<<30,episodes=0;
 while(!stop){
  if((it&1023)==0){
   double tn=now();if(tn-T0>TL)break;
   if(nunc==0){if(!valid){valid=1;memcpy(B,A,(size_t)N*K);BN=N;}memcpy(S,A,(size_t)N*K);SN=N;}
   if(it>0&&it%L==0){cyc++;
    if(cyc-lastImp>=ST){
     if(nunc!=0){clearvic();memcpy(A,S,(size_t)SN*K);N=SN;build();}else clearvic();
     {long long steps=2000000LL*(1+episodes%10);long long acc=walk(steps);int d=dropfree();
      fprintf(stderr,"walk %lld acc %lld drop %d unc %d\n",steps,acc,d,nunc);
      if(nunc==0){memcpy(S,A,(size_t)N*K);SN=N;
       if(N<BN){memcpy(B,A,(size_t)N*K);BN=N;valid=1;fprintf(stderr,"PROGRESS %d\n",-N);fflush(stderr);}}}
     W=2+(rnd()%3==0);setvic(pickvic());bestE=1<<30;lastImp=cyc;episodes++;
     fprintf(stderr,"episode %d N=%d best=%d t=%.0f W=%d\n",episodes,N,BN,tn-T0,W);}}
   double ph=(double)(it%L)/L;double T=Thi*pow(Tlo/Thi,ph);
   for(int x=0;x<OFF;x++)ex[x]=exp(-x/T);
  }
  it++;
  if(nunc==0&&nv1==0){int r=vic;clearvic();removerow(r);memcpy(B,A,(size_t)N*K);BN=N;valid=1;memcpy(S,A,(size_t)N*K);SN=N;
   fprintf(stderr,"found N=%d at %.0fs\n",N,now()-T0);fprintf(stderr,"PROGRESS %d\n",-N);fflush(stderr);
   setvic(pickvic());bestE=1<<30;lastImp=cyc;continue;}
  int E=W*nunc+nv1;if(E<bestE){bestE=E;lastImp=cyc;}
  int Ul[NS],nul=0;for(int q=0;q<NS;q++){int tt=idx(A[vic],q);if(cnt[tt]==1)Ul[nul++]=tt;}
  int tot=nunc+nul;int k=rnd()%tot;int t=k<nunc?unc[k]:Ul[k-nunc];
  int ss=t/1296,code=t%1296;int v[4]={code/216,(code/36)%6,(code/6)%6,code%6};
  int c0=sc[ss][0],c1=sc[ss][1],c2=sc[ss][2],c3=sc[ss][3];
  int cr[256],cc[256],cv[256],cd[256],nc=0,dmin=1<<30;
  for(int r=0;r<N&&nc<256;r++){if(r==vic)continue;const unsigned char*a=A[r];
   int m0=a[c0]!=v[0],m1=a[c1]!=v[1],m2=a[c2]!=v[2],m3=a[c3]!=v[3];
   if(m0+m1+m2+m3!=1)continue;int j=m0?0:m1?1:m2?2:3;
   cr[nc]=r;cc[nc]=sc[ss][j];cv[nc]=v[j];cd[nc]=delta(r,cc[nc],cv[nc]);if(cd[nc]<dmin)dmin=cd[nc];nc++;}
  if(nc>0){double wt[256],sum=0;for(int i=0;i<nc;i++){int x=cd[i]-dmin;wt[i]=x<OFF?ex[x]:0;sum+=wt[i];}
   double u=(double)(rnd()>>11)*(1.0/9007199254740992.0)*sum;int ch=nc-1;
   for(int i=0;i<nc;i++){u-=wt[i];if(u<=0){ch=i;break;}}
   setcell(cr[ch],cc[ch],cv[ch]);}
  else{int r=rnd()%N;if(r==vic)r=(r+1)%N;if(r==vic)continue;
   for(int j=0;j<4;j++){int c=sc[ss][j];if(A[r][c]!=v[j])setcell(r,c,v[j]);}}
 }
 fprintf(stderr,"done BN=%d it=%lld episodes=%d\n",BN,it,episodes);
 if(valid&&verify(B,BN))out(B,BN);else out(IN,INN);
 return 0;}
