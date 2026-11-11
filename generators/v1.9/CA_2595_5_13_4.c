// DOMAIN: k=6..16
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#define T 5
#define K 13
#define V 4
#define MAXN 20000
#define TLIM 3480.0

static volatile sig_atomic_t stopflag=0;
static void onsig(int s){(void)s;stopflag=1;}
static double t0;
static double now(){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+ts.tv_nsec*1e-9;}
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}

int NS,VT,MS;
int (*sets)[T];
int *setsOf,*posIn; // [K][MS]
int pw[T];
uint8_t (*A)[K];
uint8_t (*Bst)[K]; int BN=0;
uint32_t *code; uint16_t *cnt; uint32_t *w; int *hpos,*holes; int nh=0; int64_t cost=0;
int N;
int *tabu;

static void addhole(int h){hpos[h]=nh;holes[nh++]=h;cost+=w[h];}
static void remhole(int h){int p=hpos[h];int l=holes[--nh];holes[p]=l;hpos[l]=p;hpos[h]=-1;cost-=w[h];}

static void rebuild(){
  memset(cnt,0,sizeof(uint16_t)*(size_t)NS*VT);
  for(int r=0;r<N;r++)for(int s=0;s<NS;s++){int c=0;for(int j=0;j<T;j++)c+=A[r][sets[s][j]]*pw[j];code[(size_t)r*NS+s]=c;cnt[(size_t)s*VT+c]++;}
  nh=0;cost=0;
  for(int h=0;h<NS*VT;h++){hpos[h]=-1;if(!cnt[h])addhole(h);}
}
static void setentry(int r,int c,int v){
  int old=A[r][c]; if(old==v)return;
  for(int i=0;i<MS;i++){int s=setsOf[c*MS+i],p=posIn[c*MS+i];
    uint32_t *cp=&code[(size_t)r*NS+s]; int o=s*VT+*cp; int nn=o+(v-old)*pw[p];
    if(--cnt[o]==0)addhole(o);
    if(cnt[nn]==0)remhole(nn);
    cnt[nn]++; *cp+= (v-old)*pw[p];}
  A[r][c]=v;
}
static int64_t delta(int r,int c,int v){
  int old=A[r][c];int64_t d=0;
  for(int i=0;i<MS;i++){int s=setsOf[c*MS+i],p=posIn[c*MS+i];
    int o=s*VT+code[(size_t)r*NS+s]; int nn=o+(v-old)*pw[p];
    if(cnt[o]==1)d+=w[o]; if(cnt[nn]==0)d-=w[nn];}
  return d;
}
static void delrow(int r){
  for(int s=0;s<NS;s++){int o=s*VT+code[(size_t)r*NS+s];if(--cnt[o]==0)addhole(o);}
  int l=N-1;
  if(r!=l){memcpy(A[r],A[l],K);memcpy(&code[(size_t)r*NS],&code[(size_t)l*NS],sizeof(uint32_t)*NS);
    for(int c=0;c<K;c++)tabu[r*K+c]=tabu[l*K+c];}
  N--;
}
static void addrow(const uint8_t *row){
  memcpy(A[N],row,K);
  for(int s=0;s<NS;s++){int c=0;for(int j=0;j<T;j++)c+=row[sets[s][j]]*pw[j];code[(size_t)N*NS+s]=c;int h=s*VT+c;if(cnt[h]==0)remhole(h);cnt[h]++;}
  for(int c=0;c<K;c++)tabu[N*K+c]=0;
  N++;
}
static void savebest(){memcpy(Bst,A,(size_t)N*K);BN=N;}
static void printbest(){
  for(int r=0;r<BN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",Bst[r][c]);printf("\n");}
  fflush(stdout);
}
static int timeup(){return stopflag||now()-t0>TLIM;}

static long long iter=0;
// repair: returns 1 if holes reach 0
static int repair(long long stall){
  int candr[MAXN>4096?4096:MAXN],candc[4096],candv[4096];
  int bestnh=nh; long long last=iter;
  while(nh>0){
    if((iter&1023)==0&&timeup())return 0;
    if(iter-last>stall)return 0;
    iter++;
    int h=holes[rnd()%nh]; int s=h/VT,pat=h%VT; int dg[T];
    for(int j=0;j<T;j++)dg[j]=(pat/pw[j])%V;
    int nc=0;
    for(int r=0;r<N&&nc<4096;r++){int mm=0,mc=-1;
      for(int j=0;j<T;j++)if(A[r][sets[s][j]]!=dg[j]){mm++;mc=j;if(mm>1)break;}
      if(mm==1){candr[nc]=r;candc[nc]=sets[s][mc];candv[nc]=dg[mc];nc++;}}
    if(nc==0){int r=rnd()%N;for(int j=0;j<T;j++)setentry(r,sets[s][j],dg[j]);continue;}
    int64_t bd=INT64_MAX;int bi=-1,ties=0;
    for(int i=0;i<nc;i++){int64_t d=delta(candr[i],candc[i],candv[i]);
      int tb=tabu[candr[i]*K+candc[i]]>iter;
      if(tb&&!(d<0&&cost+d==0))continue;
      if(d<bd){bd=d;bi=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)bi=i;}}
    if(bi<0){bi=rnd()%nc;bd=delta(candr[bi],candc[bi],candv[bi]);}
    int r=candr[bi],c=candc[bi],oldv=A[r][c];
    setentry(r,c,candv[bi]);
    tabu[r*K+c]=iter+5+rnd()%10; (void)oldv;
    if(bd>=0){for(int i=0;i<nh;i++)w[holes[i]]++;cost+=nh;}
    if(nh<bestnh){bestnh=nh;last=iter;}
  }
  return 1;
}

// Metropolis repair on unweighted hole count (w must be all 1)
static int sarepair(double tbud){
  static int candr[4096],candc[4096],candv[4096];
  double ts=now(),tend=ts+tbud; long long it=0; int bestnh=nh; double Tm=0.4;
  while(nh>0){
    if((it&255)==0){double tt=now();if(timeup()||tt>tend)return 0;
      Tm=0.25+0.35*(0.5+0.5*sin((tt-ts)*0.2));
      if(nh>bestnh+60)return 0;}
    it++;iter++;
    int h=holes[rnd()%nh]; int s=h/VT,pat=h%VT; int dg[T];
    for(int j=0;j<T;j++)dg[j]=(pat/pw[j])%V;
    int nc=0;
    for(int r=0;r<N&&nc<4096;r++){int mm=0,mc=-1;
      for(int j=0;j<T;j++)if(A[r][sets[s][j]]!=dg[j]){mm++;mc=j;if(mm>1)break;}
      if(mm==1){candr[nc]=r;candc[nc]=sets[s][mc];candv[nc]=dg[mc];nc++;}}
    if(nc==0)continue;
    int bi=-1,ties=0; int64_t bd=INT64_MAX;
    if(rnd()%10==0){bi=rnd()%nc;bd=delta(candr[bi],candc[bi],candv[bi]);}
    else for(int i=0;i<nc;i++){int64_t d=delta(candr[i],candc[i],candv[i]);
      if(tabu[candr[i]*K+candc[i]]>iter&&!(d<0&&nh+d==0))continue;
      if(d<bd){bd=d;bi=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)bi=i;}}
    if(bi<0)continue;
    if(bd>0&&exp(-(double)bd/Tm)*4294967296.0<(double)(rnd()&0xffffffffULL))continue;
    int r=candr[bi],c=candc[bi];
    setentry(r,c,candv[bi]);
    tabu[r*K+c]=iter+2+rnd()%4;
    if(nh<bestnh)bestnh=nh;
  }
  return 1;
}

static int *ulist=NULL;
#define GINV 1000000
static int dgain(int r,int r2,int c,int v){
  int old=A[r2][c];int g=0;
  for(int i=0;i<MS;i++){int s=setsOf[c*MS+i],p=posIn[c*MS+i];
    int cr2=code[(size_t)r2*NS+s]; int o=s*VT+cr2; int nn=o+(v-old)*pw[p];
    if(cnt[o]==1)return GINV;
    int cr=code[(size_t)r*NS+s];
    if(cnt[o]==2&&cr==cr2)g++;
    if(cnt[nn]==1&&cr==nn-s*VT)g--;}
  return g;
}
static int uqrow(int r){int u=0;for(int s=0;s<NS;s++)if(cnt[(size_t)s*VT+code[(size_t)r*NS+s]]==1)u++;return u;}
// plateau walk on valid arrays reducing unique coverage of row r
static int drain(int r,double tbud){
  if(!ulist)ulist=malloc(sizeof(int)*NS);
  double tend=now()+tbud; long long it=0;
  int uqr=uqrow(r);
  while(uqr>0){
    if((it&255)==0){if(timeup()||now()>tend)break;}
    it++;
    if(rnd()&1){
      int nu=0;for(int s=0;s<NS;s++)if(cnt[(size_t)s*VT+code[(size_t)r*NS+s]]==1)ulist[nu++]=s;
      uqr=nu;if(nu==0)break;
      int s=ulist[rnd()%nu];int pat=code[(size_t)r*NS+s];int dg[T];
      for(int j=0;j<T;j++)dg[j]=(pat/pw[j])%V;
      int bg=GINV,br=-1,bc=-1,bv=0,ties=0;
      for(int r2=0;r2<N;r2++){if(r2==r)continue;int mm=0,mc=-1;
        for(int j=0;j<T;j++)if(A[r2][sets[s][j]]!=dg[j]){mm++;mc=j;if(mm>1)break;}
        if(mm!=1)continue;
        int c=sets[s][mc];int g=dgain(r,r2,c,dg[mc]);
        if(g>=GINV)continue;
        if(g<bg){bg=g;br=r2;bc=c;bv=dg[mc];ties=1;}else if(g==bg){ties++;if(rnd()%ties==0){br=r2;bc=c;bv=dg[mc];}}}
      if(br<0)continue;
      if(bg<=0||(bg==1&&rnd()%16==0)){setentry(br,bc,bv);uqr+=bg;}
    } else {
      int r2=rnd()%N;if(r2==r)continue;int c=rnd()%K;int v=rnd()%(V-1);if(v>=A[r2][c])v++;
      int g=dgain(r,r2,c,v);
      if(g<=0){setentry(r2,c,v);uqr+=g;}
    }
  }
  return uqrow(r);
}

static int readfile(const char*fn){
  FILE*f=fopen(fn,"r");if(!f)return 0;
  char line[1<<16];int n=0;
  while(fgets(line,sizeof line,f)){int vals[64];int m=0;char*p=line;
    while(*p){while(*p==' '||*p=='\t'||*p==','||*p=='\r'||*p=='\n')p++;if(!*p)break;
      char*e;long x=strtol(p,&e,10);if(e==p){p++;continue;}if(m<64)vals[m]=x;m++;p=e;}
    if(m==0)continue; if(m!=K){fclose(f);return 0;}
    if(n>=MAXN){fclose(f);return 0;}
    for(int c=0;c<K;c++){if(vals[c]<0||vals[c]>=V)vals[c]=rnd()%V;A[n][c]=vals[c];}
    n++;}
  fclose(f);N=n;return n>0;
}

int main(){
  signal(SIGTERM,onsig);signal(SIGINT,onsig);
  t0=now(); rs^=(uint64_t)time(NULL);
  VT=1;for(int j=0;j<T;j++){pw[j]=VT;VT*=V;}
  // enumerate sets
  {int idx[T];NS=0;for(int j=0;j<T;j++)idx[j]=j;
    int cap=1;{double x=1;for(int i=0;i<T;i++)x=x*(K-i)/(i+1);cap=(int)(x+0.5);}
    sets=malloc(sizeof(int[T])*cap);
    while(1){for(int j=0;j<T;j++)sets[NS][j]=idx[j];NS++;
      int j=T-1;while(j>=0&&idx[j]==K-T+j)j--; if(j<0)break; idx[j]++;for(int l=j+1;l<T;l++)idx[l]=idx[l-1]+1;}}
  MS=NS*T/K;
  setsOf=malloc(sizeof(int)*K*MS);posIn=malloc(sizeof(int)*K*MS);
  {int *fill=calloc(K,sizeof(int));for(int s=0;s<NS;s++)for(int j=0;j<T;j++){int c=sets[s][j];setsOf[c*MS+fill[c]]=s;posIn[c*MS+fill[c]]=j;fill[c]++;}free(fill);}
  A=malloc((size_t)MAXN*K);Bst=malloc((size_t)MAXN*K);
  code=malloc(sizeof(uint32_t)*(size_t)MAXN*NS);
  cnt=malloc(sizeof(uint16_t)*(size_t)NS*VT);w=malloc(sizeof(uint32_t)*(size_t)NS*VT);
  hpos=malloc(sizeof(int)*(size_t)NS*VT);holes=malloc(sizeof(int)*(size_t)NS*VT);
  tabu=calloc((size_t)MAXN*K,sizeof(int));
  for(int i=0;i<NS*VT;i++)w[i]=1;
  const char*fn=getenv("BEST_CANDIDATE");
  int ok=fn&&readfile(fn);
  if(ok){rebuild();
    if(nh>0){ // try to repair the input if invalid
      repair(5000000);}
  }
  if(!ok||nh>0){ // greedy fallback
    N=0;nh=0;cost=0;memset(cnt,0,sizeof(uint16_t)*(size_t)NS*VT);
    for(int h=0;h<NS*VT;h++){hpos[h]=-1;addhole(h);}
    uint8_t row[K],br[K];
    while(nh>0&&N<MAXN-1){int bg=-1;
      for(int tr=0;tr<30;tr++){int h=holes[rnd()%nh];int s=h/VT,pat=h%VT;
        for(int c=0;c<K;c++)row[c]=rnd()%V;
        for(int j=0;j<T;j++)row[sets[s][j]]=(pat/pw[j])%V;
        int g=0;for(int s2=0;s2<NS;s2++){int cc=0;for(int j=0;j<T;j++)cc+=row[sets[s2][j]]*pw[j];if(!cnt[s2*VT+cc])g++;}
        if(g>bg){bg=g;memcpy(br,row,K);}}
      addrow(br);}
  }
  for(int i=0;i<NS*VT;i++)w[i]=1;
  savebest();
  fprintf(stderr,"start N=%d holes=%d t=%.1f\n",N,nh,now()-t0);
  fprintf(stderr,"PROGRESS %d\n",-BN);
  int attempt=0;
  int *uq=malloc(sizeof(int)*MAXN),*ord=malloc(sizeof(int)*MAXN);
  while(!timeup()){
    // unique counts
    for(int r=0;r<N;r++){int u=0;for(int s=0;s<NS;s++)if(cnt[(size_t)s*VT+code[(size_t)r*NS+s]]==1)u++;uq[r]=u;ord[r]=r;}
    // partial selection: pick attempt-th smallest (with jitter)
    int kk=attempt%((N<50)?N:50);
    for(int i=0;i<=kk;i++){int m=i;for(int j=i+1;j<N;j++)if(uq[ord[j]]<uq[ord[m]])m=j;int tmp=ord[i];ord[i]=ord[m];ord[m]=tmp;}
    int r=ord[kk];
    fprintf(stderr,"N=%d delete row uq=%d (rank %d) t=%.1f\n",N,uq[r],kk,now()-t0);
    int u2=drain(r,30.0);
    if(nh==0&&N==BN)savebest();
    fprintf(stderr,"drained uq %d -> %d t=%.1f\n",uq[r],u2,now()-t0);
    if(nh!=0){memcpy(A,Bst,(size_t)BN*K);N=BN;rebuild();attempt++;continue;}
    if(u2>3){attempt++;continue;}
    delrow(r);
    int okr=(nh==0)?1:sarepair(200.0);
    if(okr){savebest();attempt=0;
      fprintf(stderr,"VALID N=%d t=%.1f\n",N,now()-t0);fprintf(stderr,"PROGRESS %d\n",-BN);
    } else {
      if(timeup())break;
      fprintf(stderr,"fail nh=%d, restore\n",nh);
      memcpy(A,Bst,(size_t)BN*K);N=BN;rebuild();attempt++;
      memset(tabu,0,sizeof(int)*(size_t)MAXN*K);
    }
  }
  printbest();
  return 0;
}
