#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#define T 3
#define K 15
#define V 12
#if T!=3
#error "this generator is strength-3 specific"
#endif
#define VT (V*V*V)
#define PC ((K-1)*(K-2)/2)
#define NTRI (K*(K-1)*(K-2)/6)
#define NTUP (NTRI*VT)
#define MAXN 6000
#define TLIM 3480.0
#define ALIM 420.0
#define TPER 40.0
#define TMAX 1.5
#define TMIN 0.3

static volatile sig_atomic_t stopf=0;
static void hnd(int s){(void)s;stopf=1;}
static int tri[NTRI][3];
static int ct_tid[K][PC],ct_o1[K][PC],ct_o2[K][PC],ct_mc[K][PC],ct_m1[K][PC],ct_m2[K][PC];
static int cnt[NTUP],w[NTUP],hpos[NTUP],holes[NTUP],nh;
static unsigned char A[MAXN][K],B[MAXN][K];
static int N,BN;
static long long tabu[MAXN][K];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addh(int i){hpos[i]=nh;holes[nh++]=i;}
static inline void remh(int i){int p=hpos[i];int l=holes[--nh];holes[p]=l;hpos[l]=p;}
static inline int tidx(int t,unsigned char*x){return t*VT+x[tri[t][0]]*V*V+x[tri[t][1]]*V+x[tri[t][2]];}
static void rebuild(void){
  memset(cnt,0,sizeof(cnt));nh=0;
  for(int r=0;r<N;r++)for(int t=0;t<NTRI;t++)cnt[tidx(t,A[r])]++;
  for(int i=0;i<NTUP;i++)if(!cnt[i])addh(i);
}
static inline long long delta(int r,int c,int nv){
  int ov=A[r][c];long long g=0;
  for(int j=0;j<PC;j++){
    int base=ct_tid[c][j]*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov,in=base+ct_mc[c][j]*nv;
    if(cnt[io]==1)g-=w[io];
    if(cnt[in]==0)g+=w[in];
  }
  return g;
}
static inline void apply(int r,int c,int nv){
  int ov=A[r][c];if(ov==nv)return;
  for(int j=0;j<PC;j++){
    int base=ct_tid[c][j]*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov,in=base+ct_mc[c][j]*nv;
    if(--cnt[io]==0)addh(io);
    if(cnt[in]++==0)remh(in);
  }
  A[r][c]=nv;
}
static void delrow(int r){
  for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(--cnt[i]==0)addh(i);}
  memcpy(A[r],A[N-1],K);N--;
}
static void outbest(void){
  for(int r=0;r<BN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",B[r][c]);printf("\n");}
  fflush(stdout);
}
#define CHD 6
#define NODECAP 4000
#define CHMAXH 6
#define CCAP 512
static int stk_r[CHD+2],stk_c[CHD+2],sd;static long nodes,chainok=0;
static inline int lossinfo(int r,int c,int nv,int*lost){
  int ov=A[r][c];int L=0;
  for(int j=0;j<PC;j++){
    int base=ct_tid[c][j]*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov;
    if(cnt[io]==1){L++;*lost=io;if(L>1)return L;}
  }
  return L;
}
static int chain(int u,int depth){
  if(++nodes>NODECAP)return 0;
  if(cnt[u]>0)return 1;
  int t=u/VT,v=u%VT;
  int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
  int qr[CCAP],qc[CCAP],qv[CCAP],ql[CCAP],qL[CCAP],nq=0;
  for(int r=0;r<N&&nq<CCAP;r++){
    unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
    if(m1+m2+m3!=2)continue;
    int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
    int skip=0;for(int s=0;s<sd;s++)if(stk_r[s]==r&&stk_c[s]==c){skip=1;break;}
    if(skip)continue;
    int lost=-1;int L=lossinfo(r,c,nv,&lost);
    if(L==0){apply(r,c,nv);return 1;}
    if(L==1){qr[nq]=r;qc[nq]=c;qv[nq]=nv;ql[nq]=lost;qL[nq]=L;nq++;}
  }
  if(depth==0)return 0;
  for(int i=0;i<nq;i++){int j=i+rnd()%(nq-i);
    int a;a=qr[i];qr[i]=qr[j];qr[j]=a;a=qc[i];qc[i]=qc[j];qc[j]=a;a=qv[i];qv[i]=qv[j];qv[j]=a;a=ql[i];ql[i]=ql[j];ql[j]=a;a=qL[i];qL[i]=qL[j];qL[j]=a;}
  for(int i=0;i<nq;i++){
    if(nodes>NODECAP)return 0;
    int r=qr[i],c=qc[i],ov=A[r][c];
    apply(r,c,qv[i]);stk_r[sd]=r;stk_c[sd]=c;sd++;
    if(cnt[ql[i]]>0){return 1;}
    if(chain(ql[i],depth-1))return 1;
    sd--;apply(r,c,ov);
  }
  return 0;
}
static void trychains(void){
  int hs[CHMAXH+2],m=nh<CHMAXH?nh:CHMAXH;
  for(int i=0;i<m;i++)hs[i]=holes[i];
  for(int i=0;i<m;i++){
    if(cnt[hs[i]]>0)continue;
    nodes=0;sd=0;
    if(chain(hs[i],CHD))chainok++;
  }
}
static int pickrow(int spread){
  static int u[MAXN];static int ord[MAXN];
  for(int r=0;r<N;r++){int s=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1)s++;u[r]=s;ord[r]=r;}
  for(int i=0;i<N;i++){int j=i+rnd()%(N-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  int m=spread<N?spread:N;
  for(int i=0;i<m;i++){int b=i;for(int j=i+1;j<N;j++)if(u[ord[j]]<u[ord[b]])b=j;int x=ord[i];ord[i]=ord[b];ord[b]=x;}
  int k=rnd()%m;
  fprintf(stderr,"delete row uniq=%d (min %d)\n",u[ord[k]],u[ord[0]]);
  return ord[k];
}
int main(void){
  signal(SIGTERM,hnd);signal(SIGINT,hnd);
  double t0=now();rs^=(unsigned long long)time(NULL)*2654435761ULL;
  int n=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int d=b+1;d<K;d++){tri[n][0]=a;tri[n][1]=b;tri[n][2]=d;n++;}
  int pcn[K];memset(pcn,0,sizeof(pcn));
  int mul[3]={V*V,V,1};
  for(int t=0;t<NTRI;t++)for(int p=0;p<3;p++){
    int c=tri[t][p],j=pcn[c]++;int q1=(p+1)%3,q2=(p+2)%3;
    ct_tid[c][j]=t;ct_mc[c][j]=mul[p];ct_o1[c][j]=tri[t][q1];ct_m1[c][j]=mul[q1];ct_o2[c][j]=tri[t][q2];ct_m2[c][j]=mul[q2];
  }
  const char*fn=getenv("BEST_CANDIDATE");N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int x,c=0;while(N<MAXN&&fscanf(f,"%d",&x)==1){A[N][c++]=(unsigned char)(((x%V)+V)%V);if(c==K){c=0;N++;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 0;}
  rebuild();
  BN=N;memcpy(B,A,(size_t)N*K);
  int bestvalid=(nh==0);
  fprintf(stderr,"loaded N=%d holes=%d\n",N,nh);
  if(!bestvalid)BN=0;
  for(int i=0;i<NTUP;i++)w[i]=1;
  memset(tabu,0,sizeof(tabu));
  long long it=0;double ast=now();int spread=1;int minh=1<<30;double Temp=TMAX;
  static int cr[MAXN],cc[MAXN],cv[MAXN];static long long cdl[MAXN];static double cw[MAXN];
  while(!stopf){
    if(nh==0){
      if(!bestvalid||N<BN){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;fprintf(stderr,"PROGRESS %d\n",-N);fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);}
      int r=pickrow(spread);delrow(r);
      for(int i=0;i<NTUP;i++)w[i]=1;
      memset(tabu,0,sizeof(tabu));ast=now();minh=nh;spread=1;
      fprintf(stderr,"N=%d holes=%d\n",N,nh);
      if(nh<=CHMAXH){trychains();fprintf(stderr,"after chains holes=%d (chainok=%ld)\n",nh,chainok);minh=nh;}
      if(nh==0)continue;
    }
    it++;
    if((it&1023)==0){
      double tn=now();
      if(tn-t0>TLIM)break;
      {double ph=fmod(tn-ast,TPER)/TPER;Temp=TMAX*pow(TMIN/TMAX,ph);}
      if((it&((1<<20)-1))==0)fprintf(stderr,"N=%d nh=%d minh=%d T=%.2f t=%.0f\n",N,nh,minh,Temp,tn-t0);
      if(tn-ast>ALIM&&bestvalid){
        fprintf(stderr,"stall N=%d minh=%d, restore\n",N,minh);
        N=BN;memcpy(A,B,(size_t)N*K);rebuild();spread=20;
        int r=pickrow(spread);delrow(r);
        for(int i=0;i<NTUP;i++)w[i]=1;memset(tabu,0,sizeof(tabu));ast=now();minh=nh;spread=1;
        continue;
      }
    }
    if((it&4095)==0&&nh<=CHMAXH){
      trychains();if(nh<minh)minh=nh;
      if(nh==0)continue;
    }
    int h=holes[rnd()%nh];int t=h/VT,v=h%VT;
    int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
    int nc=0;long long bd=-(1LL<<60);
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
      if(m1+m2+m3!=2)continue;
      int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
      long long d=delta(r,c,nv);
      cr[nc]=r;cc[nc]=c;cv[nc]=nv;cdl[nc]=d;nc++;
      if(d>bd)bd=d;
    }
    if(nc==0){
      int r=rnd()%N;apply(r,c1,v1);apply(r,c2,v2);apply(r,c3,v3);
    }else{
      double tot=0;
      for(int i=0;i<nc;i++){double e=(double)(cdl[i]-bd)/Temp;cw[i]=e<-40?0:exp(e);tot+=cw[i];}
      double u=(double)(rnd()>>11)*(1.0/9007199254740992.0)*tot;int k=0;
      for(;k<nc-1;k++){u-=cw[k];if(u<=0)break;}
      apply(cr[k],cc[k],cv[k]);
    }
    if(nh<minh)minh=nh;
  }
  if(nh==0&&(!bestvalid||N<BN)){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;}
  if(!bestvalid){BN=N;memcpy(B,A,(size_t)N*K);}
  fprintf(stderr,"final N=%d\n",BN);
  outbest();
  return 0;
}
