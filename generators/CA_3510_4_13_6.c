#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#define T 4
#define K 13
#define V 6
_Static_assert(T==4,"code is written for strength 4");
#define W2 (V*V)
#define W3 (V*V*V)
#define NT (V*V*V*V)
#define NS (K*(K-1)*(K-2)*(K-3)/24)
#define PC ((K-1)*(K-2)*(K-3)/6)
#define NG (NS*NT)
#define MAXN 8000
#define MAXM (2*K)
#define TL 3450.0
#define PG 0.2
typedef struct{int gb;unsigned char o0,o1,o2,pad;unsigned short w0,w1,w2,wc;}E;
typedef struct{int n,ngid,nfree;unsigned char A[MAXN][K],B[MAXN][K];int grp[MAXN],rot[MAXN],msz[MAXN],fr[MAXN];int mem[MAXN][MAXM];}St;
static St S,CK;
static int scol[NS][4];static E ce[K][PC];
static int cnt[NG],hl[NG],hp[NG],nh;
static unsigned char Best[MAXN][K];static int bestN=0;
static double lnt[4096],t0,hard;
static volatile sig_atomic_t stopf=0;
static void onterm(int s){(void)s;stopf=1;}
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double rnd01(void){return (rnd()>>11)*(1.0/9007199254740992.0);}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int tidx(const unsigned char*x,int s){int*c=scol[s];return x[c[0]]*W3+x[c[1]]*W2+x[c[2]]*V+x[c[3]];}
static inline void addh(int t){hp[t]=nh;hl[nh++]=t;}
static inline void delh(int t){int p=hp[t],l=hl[--nh];hl[p]=l;hp[l]=p;hp[t]=-1;}
static inline void setcell(int j,int c,int b){
  unsigned char*x=S.A[j];int a=x[c];if(a==b)return;
  E*ec=ce[c];
  for(int k=0;k<PC;k++){E*e=&ec[k];int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    int ia=base+a*e->wc,ib=base+b*e->wc;
    if(--cnt[ia]==0)addh(ia);
    if(cnt[ib]++==0)delh(ib);}
  x[c]=(unsigned char)b;
}
static inline int delta1(int j,int c,int b){
  unsigned char*x=S.A[j];int a=x[c];if(a==b)return 0;
  E*ec=ce[c];int d=0;
  for(int k=0;k<PC;k++){E*e=&ec[k];int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    d+=(cnt[base+a*e->wc]==1)-(cnt[base+b*e->wc]==0);}
  return d;
}
static void rebuild(void){
  memset(cnt,0,sizeof cnt);
  for(int j=0;j<S.n;j++)for(int s=0;s<NS;s++)cnt[s*NT+tidx(S.A[j],s)]++;
  nh=0;for(int t=0;t<NG;t++){if(cnt[t]==0){hp[t]=nh;hl[nh++]=t;}else hp[t]=-1;}
}
static int newg(void){int g=S.nfree?S.fr[--S.nfree]:S.ngid++;S.msz[g]=0;return g;}
static void freeg(int g){S.fr[S.nfree++]=g;}
static void memrm(int g,int j){
  for(int k=0;k<S.msz[g];k++)if(S.mem[g][k]==j){S.mem[g][k]=S.mem[g][--S.msz[g]];break;}
  if(S.msz[g]==0)freeg(g);
}
static void detach(int j){
  int g=S.grp[j];if(S.msz[g]<=1)return;
  memrm(g,j);int h=newg();memcpy(S.B[h],S.A[j],K);S.mem[h][0]=j;S.msz[h]=1;S.grp[j]=h;S.rot[j]=0;
}
static void gmove(int g,int c0,int b){
  for(int k=0;k<S.msz[g];k++){int j=S.mem[g][k];setcell(j,(c0+S.rot[j])%K,b);}
  S.B[g][c0]=(unsigned char)b;
}
static void delrow(int j){
  unsigned char*x=S.A[j];
  for(int s=0;s<NS;s++){int t=s*NT+tidx(x,s);if(--cnt[t]==0)addh(t);}
  memrm(S.grp[j],j);
  int L=S.n-1;
  if(j!=L){
    memcpy(S.A[j],S.A[L],K);int g=S.grp[L];S.grp[j]=g;S.rot[j]=S.rot[L];
    for(int k=0;k<S.msz[g];k++)if(S.mem[g][k]==L){S.mem[g][k]=j;break;}
  }
  S.n--;
}
static unsigned long long gkey[MAXN];
static int cmpk(const void*a,const void*b){unsigned long long x=gkey[*(const int*)a],y=gkey[*(const int*)b];return x<y?-1:x>y;}
static void group_init(void){
  static int ord[MAXN];
  S.ngid=0;S.nfree=0;
  for(int j=0;j<S.n;j++){
    unsigned long long bk=~0ULL;int br=0;
    for(int r=0;r<K;r++){unsigned long long code=0;for(int c=0;c<K;c++)code=code*V+S.A[j][(c+r)%K];if(code<bk){bk=code;br=r;}}
    gkey[j]=bk;S.rot[j]=br;ord[j]=j;
  }
  qsort(ord,(size_t)S.n,sizeof(int),cmpk);
  int g=-1;unsigned long long last=~0ULL;
  for(int i=0;i<S.n;i++){
    int j=ord[i];
    if(g<0||gkey[j]!=last||S.msz[g]>=MAXM){
      g=newg();for(int c=0;c<K;c++)S.B[g][c]=S.A[j][(c+S.rot[j])%K];last=gkey[j];
    }
    S.mem[g][S.msz[g]++]=j;S.grp[j]=g;
  }
  int hist[MAXM+1];memset(hist,0,sizeof hist);int ngr=0;
  for(int h=0;h<S.ngid;h++){hist[S.msz[h]]++;ngr++;}
  fprintf(stderr,"groups=%d:",ngr);for(int m=1;m<=MAXM;m++)if(hist[m])fprintf(stderr," %dx%d",hist[m],m);fprintf(stderr,"\n");
}
static int verify(unsigned char (*A)[K],int n){
  static unsigned char seen[NT];
  for(int s=0;s<NS;s++){
    memset(seen,0,sizeof seen);int d=0;
    for(int r=0;r<n;r++){int id=tidx(A[r],s);if(!seen[id]){seen[id]=1;d++;}}
    if(d!=NT)return 0;
  }
  return 1;
}
static int uniq(int j){int u=0;for(int s=0;s<NS;s++)if(cnt[s*NT+tidx(S.A[j],s)]==1)u++;return u;}
static void prune0(void){for(int j=S.n-1;j>=0;j--)if(uniq(j)==0){delrow(j);fprintf(stderr,"pruned redundant row, n=%d\n",S.n);}}
static int pick_weak(int rank){
  static double sc[MAXN];static char used[MAXN];
  for(int j=0;j<S.n;j++){sc[j]=uniq(j)+0.9*rnd01();used[j]=0;}
  int pick=-1;
  for(int k=0;k<=rank;k++){
    int b=-1;for(int j=0;j<S.n;j++)if(!used[j]&&(b<0||sc[j]<sc[b]))b=j;
    if(b<0)break;used[b]=1;pick=b;
  }
  return pick;
}
static int wt[NG];
static long long tabu[MAXN][K];
static inline void apply1(int j,int c,int b){
  detach(j);setcell(j,c,b);int g=S.grp[j];S.B[g][(c-S.rot[j]+K)%K]=(unsigned char)b;
}
static int repair(double dur){
  double ts=now(),te=ts+dur;if(te>hard)te=hard;
  for(int t=0;t<NG;t++)wt[t]=1;
  memset(tabu,0,sizeof tabu);
  int best=nh,h0s=nh;long long it=0,bumps=0,fb=0;
  while(nh>0){
    if((it&63)==0){double tn=now();if(tn>=te||stopf)break;}
    it++;
    int t=hl[rnd()%(unsigned long long)nh];int s=t/NT,u=t%NT;
    int dg[4]={u/W3,(u/W2)%V,(u/V)%V,u%V};
    const int*c4=scol[s];
    int bw=1<<30,bj=-1,bc=0,bb=0,nt=0;
    for(int j=0;j<S.n;j++){
      const unsigned char*x=S.A[j];int mis=0,mk=0;
      for(int k=0;k<4;k++)if(x[c4[k]]!=dg[k]){if(++mis>1)break;mk=k;}
      if(mis!=1)continue;
      int c=c4[mk],b=dg[mk],a=x[c];
      E*ec=ce[c];int d=0,w=0;
      for(int k=0;k<PC;k++){E*e=&ec[k];int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
        int ia=base+a*e->wc,ib=base+b*e->wc;
        if(cnt[ia]==1){d++;w+=wt[ia];}
        if(cnt[ib]==0){d--;w-=wt[ib];}}
      if(tabu[j][c]>it&&nh+d>=best)continue;
      if(w<bw){bw=w;bj=j;bc=c;bb=b;nt=1;}
      else if(w==bw&&rnd()%(unsigned long long)(++nt)==0){bj=j;bc=c;bb=b;}
    }
    if(bj<0){
      for(int q=0;q<20;q++){
        int j=(int)(rnd()%(unsigned long long)S.n);
        int mm[4],nm=0;for(int k=0;k<4;k++)if(S.A[j][c4[k]]!=dg[k])mm[nm++]=k;
        if(!nm)continue;
        int k=mm[rnd()%(unsigned long long)nm];
        if(tabu[j][c4[k]]>it)continue;
        apply1(j,c4[k],dg[k]);tabu[j][c4[k]]=it+3+(long long)(rnd()%8);fb++;break;
      }
    }else{
      apply1(bj,bc,bb);
      tabu[bj][bc]=it+3+(long long)(rnd()%8);
      if(bw>=0){for(int i=0;i<nh;i++)wt[hl[i]]++;bumps++;}
    }
    if(nh<best)best=nh;
  }
  fprintf(stderr,"wrepair n=%d holes %d->%d (best %d) it=%lld bumps=%lld fb=%lld t=%.0f\n",S.n,h0s,nh,best,it,bumps,fb,now()-t0);
  return nh==0;
}
static int repair_sa(double dur){
  double ts=now(),te=ts+dur;if(te>hard)te=hard;
  int best=nh,h0s=nh;double cyc=12.0,cs=ts,T0=0.7,T1=0.05,Tc=T0;
  long long it=0,ga=0,sa=0;
  while(nh>0){
    if((it&255)==0){double tn=now();if(tn>=te||stopf)break;double f=(tn-cs)/cyc;if(f>=1){cs=tn;f=0;}Tc=T0*pow(T1/T0,f);}
    it++;
    int t=hl[rnd()%(unsigned long long)nh];int s=t/NT,u=t%NT;
    int dg[4]={u/W3,(u/W2)%V,(u/V)%V,u%V};
    if(rnd01()<PG){
      for(int tr=0;tr<6;tr++){
        int j=(int)(rnd()%(unsigned long long)S.n),g=S.grp[j];
        if(S.msz[g]<2)continue;
        int mm[4],nm=0;for(int k=0;k<4;k++)if(S.A[j][scol[s][k]]!=dg[k])mm[nm++]=k;
        if(!nm)continue;
        int k=mm[rnd()%(unsigned long long)nm];int c=scol[s][k],b=dg[k];
        int c0=(c-S.rot[j]+K)%K,a=S.B[g][c0];
        int h0=nh;gmove(g,c0,b);int d=nh-h0;
        if(d>0&&d>=Tc*(1.0+S.msz[g]/4.0)*lnt[rnd()&4095])gmove(g,c0,a);else ga++;
        break;
      }
    }else{
      int bd=1<<30,bj=-1,bc=0,bb=0;
      for(int q=0;q<8;q++){
        int j=(int)(rnd()%(unsigned long long)S.n);
        int mm[4],nm=0;for(int k=0;k<4;k++)if(S.A[j][scol[s][k]]!=dg[k])mm[nm++]=k;
        if(!nm)continue;
        int k=mm[rnd()%(unsigned long long)nm];int c=scol[s][k],b=dg[k];
        int d=delta1(j,c,b);
        if(d<bd||(d==bd&&(rnd()&1))){bd=d;bj=j;bc=c;bb=b;}
      }
      if(bj>=0&&(bd<=0||bd<Tc*lnt[rnd()&4095])){
        detach(bj);setcell(bj,bc,bb);
        int g=S.grp[bj];S.B[g][(bc-S.rot[bj]+K)%K]=(unsigned char)bb;
        sa++;
      }
    }
    if(nh<best)best=nh;
  }
  fprintf(stderr,"repair n=%d holes %d->%d (best %d) it=%lld ga=%lld sa=%lld t=%.0f\n",S.n,h0s,nh,best,it,ga,sa,now()-t0);
  return nh==0;
}
static void record(void){
  if(verify(S.A,S.n)&&(bestN==0||S.n<bestN)){
    bestN=S.n;memcpy(Best,S.A,(size_t)S.n*K);
    fprintf(stderr,"NEW BEST N=%d t=%.0f\n",bestN,now()-t0);
    fprintf(stderr,"PROGRESS %d\n",-bestN);fflush(stderr);
  }
}
int main(void){
  t0=now();hard=t0+TL;
  signal(SIGTERM,onterm);
  {struct timespec ts;clock_gettime(CLOCK_REALTIME,&ts);rs^=(unsigned long long)ts.tv_nsec*0x2545F4914F6CDD1DULL^(unsigned long long)ts.tv_sec;if(!rs)rs=1;for(int i=0;i<16;i++)rnd();}
  for(int i=0;i<4096;i++)lnt[i]=-log((i+0.5)/4096.0);
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={W3,W2,V,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){
    int c=scol[s][i],j=cc[c]++;E*e=&ce[c][j];
    e->gb=s*NT;e->wc=(unsigned short)wts[i];e->pad=0;
    int oc[3],ow[3],q=0;
    for(int k=0;k<4;k++)if(k!=i){oc[q]=scol[s][k];ow[q]=wts[k];q++;}
    e->o0=(unsigned char)oc[0];e->o1=(unsigned char)oc[1];e->o2=(unsigned char)oc[2];
    e->w0=(unsigned short)ow[0];e->w1=(unsigned short)ow[1];e->w2=(unsigned short)ow[2];
  }
  const char*fn=getenv("BEST_CANDIDATE");
  int n0=0;
  if(fn){
    FILE*f=fopen(fn,"r");
    if(f){int v,c=0;while(n0<MAXN&&fscanf(f,"%d",&v)==1){S.A[n0][c++]=(unsigned char)(v<0?0:v>=V?V-1:v);if(c==K){c=0;n0++;}}fclose(f);}
  }
  S.n=n0;
  if(n0==0){fprintf(stderr,"no input\n");return 1;}
  static unsigned char In[MAXN][K];memcpy(In,S.A,(size_t)n0*K);
  group_init();
  rebuild();
  fprintf(stderr,"input N=%d holes=%d\n",n0,nh);
  if(nh>0)repair(600);
  if(nh==0){prune0();record();}
  double dur=60;int rank=0;
  while(bestN>0&&!stopf&&now()<hard-5){
    CK=S;
    int j=pick_weak(rank);if(j<0)break;
    delrow(j);
    double d=dur;if(d>hard-now())d=hard-now();
    if(repair(d)){prune0();record();rank=0;continue;}
    S=CK;rebuild();
    rank++;if(rank>=4){rank=0;dur*=1.5;if(dur>600)dur=600;}
  }
  if(bestN==0){fprintf(stderr,"nothing valid, printing input\n");bestN=n0;memcpy(Best,In,(size_t)n0*K);}
  fprintf(stderr,"final N=%d verified=%d t=%.0f\n",bestN,verify(Best,bestN),now()-t0);
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",Best[r][c]);printf("\n");}
  fflush(stdout);
  return 0;
}
