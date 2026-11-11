// DOMAIN: k=12..24
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#define T 3
#define K 19
#define V 15
#if T!=3
#error "this generator is strength-3 specific"
#endif
#if K>64
#error "merge masks need K<=64"
#endif
#define VT (V*V*V)
#define PC ((K-1)*(K-2)/2)
#define NTRI (K*(K-1)*(K-2)/6)
#define NTUP (NTRI*VT)
#define MAXN 12000
#define TLIM 3470.0
#define ALIM 900.0
#define SLIM 180.0
#define FRCAP 64
#define DRIFT 12
#define KICKT 1e9
#define KICKM 2
#define TMAX 1.5
#define TMIN 0.3
#define WMAX 250
#define PLAT_STALL 30.0
#define PLAT_LEVEL 1.5
#define PLAT_INIT 40.0
#define DRAIN_T 8.0
#define DRAIN_STALL 20.0
#define DRAIN_IDLE 1.5
#define CHMAXH 6
typedef unsigned long long u64;

static volatile sig_atomic_t stopf=0;
static void hnd(int s){(void)s;stopf=1;}
static int tri[NTRI][3];
static int ct_tid[K][PC],ct_o1[K][PC],ct_o2[K][PC],ct_mc[K][PC],ct_m1[K][PC],ct_m2[K][PC];
static unsigned short cnt[NTUP];static unsigned char w[NTUP];
static int hpos[NTUP],holes[NTUP],nh;
static unsigned char A[MAXN][K],B[MAXN][K];
static int N,BN;
static double t0;
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static long long Wsum=0;
static inline void addh(int i){hpos[i]=nh;holes[nh++]=i;Wsum+=w[i];}
static inline void remh(int i){int p=hpos[i];int l=holes[--nh];holes[p]=l;hpos[l]=p;Wsum-=w[i];}
static inline int tidx(int t,unsigned char*x){return t*VT+x[tri[t][0]]*V*V+x[tri[t][1]]*V+x[tri[t][2]];}
static void rebuild(void){
  memset(cnt,0,sizeof(cnt));nh=0;Wsum=0;
  for(int r=0;r<N;r++)for(int t=0;t<NTRI;t++)cnt[tidx(t,A[r])]++;
  for(int i=0;i<NTUP;i++)if(!cnt[i])addh(i);
}
static inline long long delta2(int r,int c,int nv,int*dn){
  int ov=A[r][c];long long g=0;int d=0;
  for(int j=0;j<PC;j++){
    int base=ct_tid[c][j]*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov,in=base+ct_mc[c][j]*nv;
    if(cnt[io]==1){g-=w[io];d++;}
    if(cnt[in]==0){g+=w[in];d--;}
  }
  *dn=d;return g;
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
#define CCAP 512
static int stk_r[CHD+4],stk_c[CHD+4],sd;static long nodes,chainok=0,compok=0;
static int lastu,minu;
static inline int lossinfo(int r,int c,int nv,int*lost){
  int ov=A[r][c];int L=0;(void)nv;
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
  int qr[CCAP],qc[CCAP],qv[CCAP],ql[CCAP],nq=0;
  for(int r=0;r<N&&nq<CCAP;r++){
    unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
    if(m1+m2+m3!=2)continue;
    int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
    int skip=0;for(int s=0;s<sd;s++)if(stk_r[s]==r&&stk_c[s]==c){skip=1;break;}
    if(skip)continue;
    int lost=-1;int L=lossinfo(r,c,nv,&lost);
    if(L==0){apply(r,c,nv);return 1;}
    if(L==1){qr[nq]=r;qc[nq]=c;qv[nq]=nv;ql[nq]=lost;nq++;}
  }
  if(depth==0)return 0;
  for(int i=0;i<nq;i++){int j=i+rnd()%(nq-i);
    int a;a=qr[i];qr[i]=qr[j];qr[j]=a;a=qc[i];qc[i]=qc[j];qc[j]=a;a=qv[i];qv[i]=qv[j];qv[j]=a;a=ql[i];ql[i]=ql[j];ql[j]=a;}
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
#define CPNODE 600
static long cbud;
static int compound(int u){
  if(cnt[u]>0)return 1;
  int t=u/VT,v=u%VT;
  int cs[3]={tri[t][0],tri[t][1],tri[t][2]};int vs[3]={v/(V*V),(v/V)%V,v%V};
  static int cand[MAXN],cmi[MAXN];int n=0;
  for(int r=0;r<N;r++){
    int m=0,mi=-1;for(int d=0;d<3;d++)if(A[r][cs[d]]==vs[d]){m++;mi=d;}
    if(m==1){cand[n]=r;cmi[n]=mi;n++;}
  }
  for(int i=0;i<n;i++){int j=i+rnd()%(n-i);int x=cand[i];cand[i]=cand[j];cand[j]=x;x=cmi[i];cmi[i]=cmi[j];cmi[j]=x;}
  int orig[CHMAXH+2],no=nh<CHMAXH?nh:CHMAXH;
  for(int i=0;i<no;i++)orig[i]=holes[i];
  for(int i=0;i<n&&!stopf;i++){
    int r=cand[i],mi=cmi[i];int a=(mi+1)%3,b=(mi+2)%3;
    int oa=A[r][cs[a]],ob=A[r][cs[b]];int h0=nh;
    apply(r,cs[a],vs[a]);apply(r,cs[b],vs[b]);
    if(nh<h0)return 1;
    if(nh==h0&&cbud>0&&nh<=CHMAXH+1){
      int nw=-1,nn=0;
      for(int k=0;k<nh;k++){int hh=holes[k],f=0;for(int q=0;q<no;q++)if(orig[q]==hh){f=1;break;}if(!f){nw=hh;nn++;}}
      if(nn==1){
        nodes=NODECAP-CPNODE;cbud-=CPNODE;
        sd=2;stk_r[0]=r;stk_c[0]=cs[a];stk_r[1]=r;stk_c[1]=cs[b];
        if(chain(nw,CHD-2))return 1;
      }
    }
    apply(r,cs[b],ob);apply(r,cs[a],oa);
  }
  return 0;
}
static void trychains(int comp,int maxh){
  static int hs[4096];int m=nh<maxh?nh:maxh;
  for(int i=0;i<m;i++)hs[i]=holes[i];
  for(int i=0;i<m;i++){
    if(cnt[hs[i]]>0)continue;
    nodes=0;sd=0;
    if(chain(hs[i],CHD)){chainok++;continue;}
  }
  if(comp&&nh>0&&nh<=3){
    cbud=12000;
    int m2=nh;for(int i=0;i<m2;i++)hs[i]=holes[i];
    for(int i=0;i<m2;i++){if(cnt[hs[i]]>0)continue;if(compound(hs[i]))compok++;}
  }
}
#define RRCAP 4096
#define RRTRY 24
#define LOWN 256
static long rrok=0,rrside=0;
static int lowr[LOWN],nlow=0;static double lowt=-1e9;
static int rowuniq(int r);
static void buildlow(void){
  static int u[MAXN];static int hist[NTRI+2];
  memset(hist,0,sizeof(hist));
  for(int r=0;r<N;r++){u[r]=rowuniq(r);hist[u[r]]++;}
  int th=0,acc=0;while(th<=NTRI&&acc+hist[th]<LOWN){acc+=hist[th];th++;}
  nlow=0;
  for(int r=0;r<N&&nlow<LOWN;r++)if(u[r]<th)lowr[nlow++]=r;
  for(int r=0;r<N&&nlow<LOWN;r++)if(u[r]==th)lowr[nlow++]=r;
  lowt=now();
}
static void rowsub(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(--cnt[i]==0)addh(i);}}
static void rowadd(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]++==0)remh(i);}}
static int rowresolve(int tries){
  static int H[RRCAP];static int fx[K];static unsigned char best[K],x[K];
  int imp=0;
  if(nh>0&&now()-lowt>3.0)buildlow();
  for(int a=0;a<tries&&nh>0&&!stopf;a++){
    int nh0=nh;if(nh0>RRCAP/2)return imp;
    int r=(nlow>0&&(rnd()&3))?lowr[rnd()%nlow]:(int)(rnd()%N);if(r>=N)r=rnd()%N;int no=nh0;
    rowsub(r);
    int nH=nh;if(nH>RRCAP){rowadd(r);continue;}
    memcpy(H,holes,(size_t)nH*sizeof(int));
    int bm=-1;
    for(int q=0;q<RRTRY;q++){
      for(int c=0;c<K;c++)fx[c]=-1;
      if(q%4==3){for(int i=0;i<nH;i++){int j=i+rnd()%(nH-i);int y=H[i];H[i]=H[j];H[j]=y;}}
      else{
        for(int i=0;i<no;i++){int j=i+rnd()%(no-i);int y=H[i];H[i]=H[j];H[j]=y;}
        for(int i=no;i<nH;i++){int j=i+rnd()%(nH-i);int y=H[i];H[i]=H[j];H[j]=y;}
      }
      for(int i=0;i<nH;i++){int h=H[i];int t=h/VT,v=h%VT;
        int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
        if((fx[c1]<0||fx[c1]==v1)&&(fx[c2]<0||fx[c2]==v2)&&(fx[c3]<0||fx[c3]==v3)){fx[c1]=v1;fx[c2]=v2;fx[c3]=v3;}
      }
      for(int c=0;c<K;c++)x[c]=fx[c]<0?A[r][c]:(unsigned char)fx[c];
      int m=0;
      for(int i=0;i<nH;i++){int h=H[i];int t=h/VT,v=h%VT;
        if(x[tri[t][0]]==v/(V*V)&&x[tri[t][1]]==(v/V)%V&&x[tri[t][2]]==v%V)m++;}
      if(m>bm){bm=m;memcpy(best,x,K);}
    }
    int nn=nH-bm;
    if(nn<nh0||(nn==nh0&&(rnd()&1))){memcpy(A[r],best,K);if(nn<nh0){imp++;rrok++;}else rrside++;}
    rowadd(r);
  }
  return imp;
}
static unsigned char Ck[MAXN][K];static int ckh;static double ckt,cksv,tbase;static long kicks=0;
static void cksave(void){memcpy(Ck,A,(size_t)N*K);cksv=now();}
static void ckreset(void){ckh=nh;cksave();ckt=tbase=cksv;}
static void kick(void){
  memcpy(A,Ck,(size_t)N*K);rebuild();
  for(int k=0;k<nh;k++){int hh=holes[k];if(w[hh]<WMAX)w[hh]++;}
  for(int m=0;m<KICKM;m++){int r=rnd()%N,c=rnd()%K,nv=rnd()%V;apply(r,c,nv);}
  kicks++;ckt=tbase=now();
}
static inline int lossdu(int r,int c,int nv,int*du){
  int ov=A[r][c];int d=0;
  for(int j=0;j<PC;j++){
    int base=ct_tid[c][j]*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov,in=base+ct_mc[c][j]*nv;
    if(cnt[io]==1)return 1;
    if(cnt[io]==2)d++;
    if(cnt[in]==1)d--;else if(cnt[in]==0)d++;
  }
  *du=d;return 0;
}
static inline int guided(int r,int c,int nv){
  int du;if(lossdu(r,c,nv,&du))return 0;
  if(du<=0||(du==1&&(rnd()&15)==0)){apply(r,c,nv);return 1;}
  return 0;
}
static long countU(void){long u=0;for(int i=0;i<NTUP;i++)if(cnt[i]==1)u++;return u;}
static void plateau(double secs){
  double ts=now();long acc=0,it=0;
  double rem=TLIM-(ts-t0);if(secs>rem)secs=rem;
  while(!stopf&&now()-ts<secs){
    for(int k=0;k<4096;k++){
      it++;int r=rnd()%N,c=rnd()%K,nv=rnd()%V;if(nv==A[r][c])continue;
      if(guided(r,c,nv))acc++;
    }
  }
  fprintf(stderr,"plateau %.0fs it=%ld acc=%ld nh=%d U=%ld\n",secs,it,acc,nh,countU());
}
static int rowuniq(int r){int s=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1)s++;return s;}
static void deepchains(void){
  for(int rd=0;rd<10&&nh>0&&nh<=64&&!stopf;rd++){
    int h0=nh;trychains(0,64);
    if(nh>=h0)break;
  }
  fprintf(stderr,"after deepchains holes=%d (chainok=%ld)\n",nh,chainok);
}
static int drain(int r,double secs){
  double ts=now();double rem=TLIM-(ts-t0);if(secs>rem)secs=rem;
  static int lst[NTRI];int u0=rowuniq(r),bestu=u0;long mv=0;double lastimp=ts;
  while(!stopf&&now()-ts<secs){
    int n=0;for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]==1)lst[n++]=i;}
    if(n<bestu){bestu=n;lastimp=now();}
    if(n==0)break;
    if(now()-lastimp>DRAIN_IDLE)break;
    for(int q=0;q<n;q++){int u=lst[q];if(cnt[u]!=1)continue;
      int t=u/VT,v=u%VT;
      int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
      int st=rnd()%N;
      for(int k=0;k<N;k++){int s=(st+k)%N;if(s==r)continue;
        unsigned char*x=A[s];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
        if(m1+m2+m3!=2)continue;
        int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
        int lost;if(lossinfo(s,c,nv,&lost)==0){apply(s,c,nv);mv++;break;}
      }
    }
    for(int k=0;k<30000;k++){
      int s=rnd()%N;if(s==r)continue;int c=rnd()%K,nv=rnd()%V;if(nv==A[s][c])continue;
      guided(s,c,nv);
    }
  }
  int uf=rowuniq(r);
  fprintf(stderr,"drain row uniq %d -> %d (min %d) moves %ld\n",u0,uf,bestu,mv);
  return uf;
}
static unsigned char FR[FRCAP][K],lastdel[K];static int nfr=0;
static int isfailed(int r){int m=nfr<FRCAP?nfr:FRCAP;for(int f=0;f<m;f++){int d=0;for(int c=0;c<K&&d<=3;c++)d+=A[r][c]!=FR[f][c];if(d<=3)return 1;}return 0;}
static int pickrow(int spread){
  static int u[MAXN];static int ord[MAXN];
  for(int r=0;r<N;r++){u[r]=rowuniq(r);if(nfr&&isfailed(r))u[r]+=100000;ord[r]=r;}
  for(int i=0;i<N;i++){int j=i+rnd()%(N-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  int m=spread<N?spread:N;
  for(int i=0;i<m;i++){int b=i;for(int j=i+1;j<N;j++)if(u[ord[j]]<u[ord[b]])b=j;int x=ord[i];ord[i]=ord[b];ord[b]=x;}
  int k=rnd()%m;
  lastu=u[ord[k]];minu=u[ord[0]];
  fprintf(stderr,"pick row uniq=%d (min %d)\n",u[ord[k]],u[ord[0]]);
  return ord[k];
}
#define UCAP 20
#define MLIST 400
#define PMAX 12
#define MTOP 64
static int ul[MAXN][UCAP],uc[MAXN];
static int mg_s[MTOP],mg_l[MTOP],mg_r[MTOP],nmg;static u64 mg_m[MTOP];
static int cmpuc(const void*a,const void*b){return uc[*(const int*)a]-uc[*(const int*)b];}
static int mergecand(int lim){
  static int L[MAXN];int nl=0;
  for(int r=0;r<N;r++){int k=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1){if(k<UCAP)ul[r][k]=t;k++;}uc[r]=k;}
  for(int r=0;r<N;r++)if(uc[r]<=lim&&uc[r]<=UCAP&&!(nfr&&isfailed(r)))L[nl++]=r;
  for(int i=0;i<nl;i++){int j=i+rnd()%(nl-i);int x=L[i];L[i]=L[j];L[j]=x;}
  qsort(L,nl,sizeof(int),cmpuc);
  if(nl>MLIST)nl=MLIST;
  nmg=0;
  static u64 um[1<<PMAX];
  for(int i=0;i<nl;i++)for(int j=i+1;j<nl;j++){
    int s=L[i],l=L[j];if(uc[s]>uc[l]){int x=s;s=l;l=x;}
    int p=uc[s],q=uc[l];if(p>PMAX)continue;
    u64 cs[UCAP],cl[UCAP];
    for(int k=0;k<p;k++){int t=ul[s][k];u64 m=0;for(int d=0;d<3;d++){int c=tri[t][d];if(A[s][c]!=A[l][c])m|=1ULL<<c;}cs[k]=m;}
    for(int k=0;k<q;k++){int t=ul[l][k];u64 m=0;for(int d=0;d<3;d++){int c=tri[t][d];if(A[s][c]!=A[l][c])m|=1ULL<<c;}cl[k]=m;}
    int bs=-1;u64 bm=0;um[0]=0;
    for(int S=0;S<(1<<p);S++){
      if(S)um[S]=um[S&(S-1)]|cs[__builtin_ctz(S)];
      u64 uu=um[S];int sat=__builtin_popcount(S);
      for(int k=0;k<q;k++)if(!(cl[k]&uu))sat++;
      if(sat>bs){bs=sat;bm=uu;}
    }
    int res=p+q-bs;
    if(nmg<MTOP||res<mg_r[nmg-1]){
      int k=nmg<MTOP?nmg++:MTOP-1;
      while(k>0&&mg_r[k-1]>res){mg_s[k]=mg_s[k-1];mg_l[k]=mg_l[k-1];mg_m[k]=mg_m[k-1];mg_r[k]=mg_r[k-1];k--;}
      mg_s[k]=s;mg_l[k]=l;mg_m[k]=bm;mg_r[k]=res;
    }
  }
  return nmg;
}
#define NPROBE 6
#define PROBE_T 25.0
static unsigned char PB[MAXN][K],PS[MAXN][K];
static void probe_reduce(int spread){
  static int u[MAXN],ord[MAXN];
  for(int r=0;r<N;r++){u[r]=rowuniq(r);if(nfr&&isfailed(r))u[r]+=100000;ord[r]=r;}
  for(int i=0;i<N;i++){int j=i+rnd()%(N-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  int m=NPROBE*spread;if(m>N)m=N;
  for(int i=0;i<m;i++){int b=i;for(int j=i+1;j<N;j++)if(u[ord[j]]<u[ord[b]])b=j;int x=ord[i];ord[i]=ord[b];ord[b]=x;}
  for(int i=0;i<m;i++){int j=i+rnd()%(m-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  int np=NPROBE<m?NPROBE:m;
  static int cand[NPROBE],cu[NPROBE];
  for(int i=0;i<np;i++){cand[i]=ord[i];cu[i]=u[ord[i]];}
  int N0=N;memcpy(PS,A,(size_t)N0*K);
  int bh=1<<30,bN=0,bp=-1;unsigned char bdel[K],tdel[K];double ts=now();
  for(int p=0;p<np&&!stopf;p++){
    if(p>0){memcpy(A,PS,(size_t)N0*K);N=N0;rebuild();}
    int r=cand[p];memcpy(tdel,A[r],K);
    delrow(r);lowt=-1e9;
    rowresolve(1500);
    if(nh>0)deepchains();
    fprintf(stderr,"probe %d uniq=%d -> holes %d\n",p,cu[p],nh);
    if(nh<bh){bh=nh;bN=N;bp=p;memcpy(PB,A,(size_t)N*K);memcpy(bdel,tdel,K);}
    if(nh==0)break;
    if(now()-ts>PROBE_T||TLIM-(now()-t0)<60)break;
  }
  if(bp<0){memcpy(A,PS,(size_t)N0*K);N=N0;rebuild();int r=cand[0];memcpy(lastdel,A[r],K);delrow(r);return;}
  memcpy(A,PB,(size_t)bN*K);N=bN;rebuild();memcpy(lastdel,bdel,K);lowt=-1e9;
  lastu=cu[bp];
  fprintf(stderr,"probe chose %d (uniq %d) holes %d\n",bp,cu[bp],nh);
}
static void reduce(int spread,double dt){
  if(dt>0){(void)dt;probe_reduce(spread);return;}
  int r=pickrow(spread);
  lastu=drain(r,dt);
  if(lastu<minu)minu=lastu;
  int n=mergecand(minu+8);
  if(n>0){
    int k=rnd()%(spread<n?spread:n);
    if(mg_r[k]<lastu||(mg_r[k]==lastu&&(rnd()&1))){
      int s=mg_s[k],l=mg_l[k];u64 m=mg_m[k];
      int us=uc[s],ul2=uc[l];
      for(int c=0;c<K;c++){int nv=((m>>c)&1ULL)?A[s][c]:A[l][c];apply(l,c,nv);}
      memcpy(lastdel,A[s],K);
      delrow(s);
      fprintf(stderr,"merge uniq %d+%d -> holes %d (pred %d, best %d)\n",us,ul2,nh,mg_r[k],mg_r[0]);
      return;
    }
  }
  memcpy(lastdel,A[r],K);
  delrow(r);
}
#define SWAPH 10
static int coverrow(void){
  int fx[K];for(int c=0;c<K;c++)fx[c]=-1;
  if(N>=MAXN)return 0;
  for(int k=0;k<nh;k++){int h=holes[k];int t=h/VT,v=h%VT;
    int cs[3]={tri[t][0],tri[t][1],tri[t][2]},vs[3]={v/(V*V),(v/V)%V,v%V};
    for(int d=0;d<3;d++){if(fx[cs[d]]>=0&&fx[cs[d]]!=vs[d])return 0;}
    for(int d=0;d<3;d++)fx[cs[d]]=vs[d];
  }
  unsigned char x[K];int br=rnd()%N;
  for(int c=0;c<K;c++)x[c]=fx[c]>=0?(unsigned char)fx[c]:A[br][c];
  for(int rd=0;rd<4;rd++)for(int c=0;c<K;c++){if(fx[c]>=0)continue;
    int bv=x[c],bs=-1;
    for(int v=0;v<V;v++){int s=0;
      for(int j=0;j<PC;j++){int idx=ct_tid[c][j]*VT+ct_m1[c][j]*x[ct_o1[c][j]]+ct_m2[c][j]*x[ct_o2[c][j]]+ct_mc[c][j]*v;if(cnt[idx]==1)s++;}
      if(s>bs||(s==bs&&(rnd()&1))){bs=s;bv=v;}
    }
    x[c]=(unsigned char)bv;
  }
  memcpy(A[N],x,K);N++;rowadd(N-1);
  memcpy(FR[nfr%FRCAP],x,K);nfr++;
  return nh==0;
}
#define ORB_T 1500.0
#define OSTALL 150.0
static unsigned char OB[MAXN][K];static int OBN=0,OBm=0,OBf=0,om=0,of=0;
static inline void omove(int r,int c,int nv){
  if(r<om*K){int s=r/K,i=r%K,sc=(c-i+K)%K;for(int j=0;j<K;j++)apply(s*K+j,(sc+j)%K,nv);}
  else apply(r,c,nv);
}
static void wreset(void){for(int i=0;i<NTUP;i++)w[i]=1;Wsum=nh;}
static void oreduce(int spread){
  if(N<=1)return;
  static int u[MAXN],L[MAXN];
  int lo=of>0?om*K:0,mn=1<<30,nl=0;
  for(int r=lo;r<N;r++){u[r]=rowuniq(r);if(u[r]<mn)mn=u[r];}
  for(int r=lo;r<N;r++)if(u[r]<=mn+spread-1)L[nl++]=r;
  int r=L[rnd()%nl];int ur=u[r];
  if(of>0){delrow(r);of--;}
  else{
    int s=r/K,i0=r%K;static unsigned char tmp[K][K];
    memcpy(tmp,A[s*K],K*K);
    if(s!=om-1)memcpy(A[s*K],A[(om-1)*K],K*K);
    om--;int p=om*K;for(int j=0;j<K;j++)if(j!=i0)memcpy(A[p++],tmp[j],K);
    of=K-1;N=p;rebuild();
  }
  wreset();
  fprintf(stderr,"orbit reduce uniq=%d -> N=%d (m=%d f=%d) holes=%d\n",ur,N,om,of,nh);
}
static void orbit_phase(double budget,int start,int*bv){
  double ts=now();
  if(ts-t0+budget>TLIM)budget=TLIM-(ts-t0);
  om=(start+K-1)/K;of=0;N=om*K;if(N>MAXN||om<1)return;
  for(int s=0;s<om;s++)for(int c=0;c<K;c++){int x=rnd()%V;for(int j=0;j<K;j++)A[s*K+j][(c+j)%K]=x;}
  rebuild();wreset();
  static int cr[MAXN],cc[MAXN],cv[MAXN];
  double lastv=now();long long it=0;
  while(!stopf){
    if(nh==0){
      if(!OBN||N<OBN){OBN=N;OBm=om;OBf=of;memcpy(OB,A,(size_t)N*K);}
      fprintf(stderr,"orbit valid N=%d m=%d f=%d t=%.0f\n",N,om,of,now()-t0);
      if(!*bv||N<BN){BN=N;memcpy(B,A,(size_t)N*K);*bv=1;fprintf(stderr,"PROGRESS %d\n",-N);}
      if(N<=1)break;
      oreduce(1);lastv=now();continue;
    }
    it++;
    if((it&255)==0){
      double tn=now();if(tn-ts>budget)break;
      if((it&((1<<18)-1))==0){
        for(int i=0;i<NTUP;i++)if(w[i]>1)w[i]=1+((w[i]-1)>>1);
        Wsum=0;for(int k=0;k<nh;k++)Wsum+=w[holes[k]];
        fprintf(stderr,"orbit N=%d m=%d f=%d nh=%d W=%lld t=%.0f\n",N,om,of,nh,Wsum,tn-t0);
      }
      if(tn-lastv>OSTALL&&OBN){
        om=OBm;of=OBf;N=OBN;memcpy(A,OB,(size_t)N*K);rebuild();wreset();
        oreduce(3);lastv=now();continue;
      }
    }
    int h=holes[rnd()%nh];int t=h/VT,v=h%VT;
    int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
    int nc=0;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
      if(m1+m2+m3!=2)continue;
      int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
      cr[nc]=r;cc[nc]=c;cv[nc]=nv;nc++;
    }
    if(nc==0){int r=rnd()%N;omove(r,c1,v1);omove(r,c2,v2);omove(r,c3,v3);continue;}
    long long bd=(1LL<<62);int bk=0,ties=0;
    for(int k=0;k<nc;k++){
      int r=cr[k],c=cc[k],ov=A[r][c];long long w0=Wsum;
      omove(r,c,cv[k]);long long d=Wsum-w0;omove(r,c,ov);
      if(d<bd){bd=d;bk=k;ties=1;}else if(d==bd&&rnd()%(++ties)==0)bk=k;
    }
    if(bd<0)omove(cr[bk],cc[bk],cv[bk]);
    else{
      if(w[h]<250){w[h]++;Wsum++;}
      if(bd==0||rnd()%100<3)omove(cr[bk],cc[bk],cv[bk]);
    }
  }
  fprintf(stderr,"orbit phase done best N=%d\n",OBN);
}
int main(void){
  signal(SIGTERM,hnd);signal(SIGINT,hnd);
  t0=now();rs^=(unsigned long long)time(NULL)*2654435761ULL;
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
  if(N==0){
    N=(int)(VT*3);if(N>MAXN)N=MAXN;
    for(int r=0;r<N;r++)for(int c=0;c<K;c++)A[r][c]=rnd()%V;
  }
  rebuild();
  BN=N;memcpy(B,A,(size_t)N*K);
  int bestvalid=(nh==0);
  fprintf(stderr,"loaded N=%d holes=%d\n",N,nh);
  if(!bestvalid)BN=0;
  for(int i=0;i<NTUP;i++)w[i]=1;
  if(!bestvalid){
    orbit_phase(ORB_T,VT*3,&bestvalid);
    if(bestvalid){N=BN;memcpy(A,B,(size_t)N*K);}
    rebuild();
    for(int i=0;i<NTUP;i++)w[i]=1;
  }
  if(stopf){if(!bestvalid){BN=N;memcpy(B,A,(size_t)N*K);}outbest();return 0;}
  if(bestvalid){fprintf(stderr,"initial U=%ld\n",countU());plateau(PLAT_INIT);
    if(nh==0){memcpy(B,A,(size_t)N*K);}}
  ckreset();
  long long it=0;double ast=now();int spread=1;int minh=1<<30;double Temp=TMAX;
  double limp=ast;int lastminh=1<<30;
  static int cr[MAXN],cc[MAXN],cv[MAXN];static long long cdl[MAXN];
  while(!stopf){
    if(nh==0){
      if(!bestvalid||N<BN){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);}
      plateau(PLAT_LEVEL);
      if(stopf)break;
      reduce(spread,DRAIN_T);
      if(stopf)break;
      for(int i=0;i<NTUP;i++)w[i]=1;
      ckreset();
      ast=now();minh=nh;spread=1;limp=ast;lastminh=1<<30;
      fprintf(stderr,"N=%d holes=%d\n",N,nh);
      rowresolve(3000);
      fprintf(stderr,"after rowresolve holes=%d (rrok=%ld side=%ld)\n",nh,rrok,rrside);
      deepchains();
      if(nh>0){rowresolve(1000);deepchains();}
      minh=nh;
      if(nh<=CHMAXH){trychains(1,CHMAXH);fprintf(stderr,"after chains holes=%d (chainok=%ld compok=%ld)\n",nh,chainok,compok);minh=nh;}
      if(nh==0)continue;
    }
    it++;
    if((it&1023)==0){
      double tn=now();
      if(tn-t0>TLIM)break;
      {double ph=(tn-tbase)/KICKT;if(ph>1.0)ph=1.0;Temp=TMAX*pow(TMIN/TMAX,ph);(void)Temp;}
      if((it&((1<<20)-1))==0)fprintf(stderr,"N=%d nh=%d minh=%d ckh=%d t=%.0f compok=%ld kicks=%ld\n",N,nh,minh,ckh,tn-t0,compok,kicks);
      if((it&((1<<17)-1))==0){for(int i=0;i<NTUP;i++)if(w[i]>1)w[i]=(unsigned char)(1+((w[i]-1)*3)/4);
        Wsum=0;for(int k=0;k<nh;k++)Wsum+=w[holes[k]];}
      if(tn-ckt>KICKT&&nh>0){kick();if(nh==0)continue;}
      if(nh>ckh+3*DRIFT){memcpy(A,Ck,(size_t)N*K);rebuild();if(nh==0)continue;}
      if(minh<lastminh){lastminh=minh;limp=tn;}
      if((tn-ast>ALIM||tn-limp>SLIM)&&bestvalid){
        fprintf(stderr,"stall N=%d minh=%d after %.0fs, restore+plateau (failed rows %d)\n",N,minh,tn-ast,nfr+1);
        memcpy(FR[nfr%FRCAP],lastdel,K);nfr++;
        int swapped=0;
        if(N==BN-1&&ckh<=SWAPH){
          memcpy(A,Ck,(size_t)N*K);rebuild();
          if(nh==0){swapped=1;fprintf(stderr,"valid N=%d (checkpoint) t=%.0f\n",N,now()-t0);}
          else if(nh<=SWAPH){int h0=nh;swapped=coverrow();fprintf(stderr,"swap cover %d holes -> %s N=%d\n",h0,swapped?"ok":"fail",N);}
        }
        if(swapped&&nh==0){BN=N;memcpy(B,A,(size_t)N*K);}
        else{N=BN;memcpy(A,B,(size_t)N*K);rebuild();}
        plateau(PLAT_STALL);
        if(stopf)break;
        if(nh==0){memcpy(B,A,(size_t)N*K);}
        spread=4;
        reduce(spread,DRAIN_STALL);
        if(stopf)break;
        for(int i=0;i<NTUP;i++)w[i]=1;ckreset();ast=now();minh=nh;spread=1;limp=ast;lastminh=1<<30;
        rowresolve(3000);deepchains();minh=nh;
        if(nh<=CHMAXH){trychains(1,CHMAXH);minh=nh;}
        continue;
      }
    }
    if(nh==0)continue;
    if((it&2047)==0&&nh<=40){
      rowresolve(64);
      if(nh<minh)minh=nh;
      if(nh<ckh){ckh=nh;cksave();ckt=cksv;}
      if(nh==0)continue;
    }
    if((it&4095)==0&&nh<=CHMAXH){
      trychains((it&16383)==0,CHMAXH);if(nh<minh)minh=nh;
      if(nh==0)continue;
      for(int k=0;k<nh;k++){int hh=holes[k];if(w[hh]<WMAX)w[hh]++;}
    }
    int h=holes[rnd()%nh];int t=h/VT,v=h%VT;
    int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
    int nc=0,na=0;long long bd=-(1LL<<60);int cap=ckh+DRIFT;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
      if(m1+m2+m3!=2)continue;
      int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
      int dn;long long d=delta2(r,c,nv,&dn);
      if(nh+dn>cap){nc++;continue;}
      cr[na]=r;cc[na]=c;cv[na]=nv;cdl[na]=d;na++;nc++;
      if(d>bd)bd=d;
    }
    if(nc==0){
      int r=rnd()%N;apply(r,c1,v1);apply(r,c2,v2);apply(r,c3,v3);
    }else if(na==0){
      if(w[h]<WMAX){w[h]++;Wsum++;}
    }else{
      int bk=0,ties=0;
      for(int i=0;i<na;i++)if(cdl[i]==bd){if(rnd()%(++ties)==0)bk=i;}
      if(bd>0)apply(cr[bk],cc[bk],cv[bk]);
      else if(bd==0&&(rnd()%100)<30)apply(cr[bk],cc[bk],cv[bk]);
      else{
        if(w[h]<WMAX){w[h]++;Wsum++;}
        if((rnd()%1000)<5){int k=rnd()%na;apply(cr[k],cc[k],cv[k]);}
      }
    }
    if(nh<minh)minh=nh;
    if(nh<ckh){ckh=nh;cksave();ckt=cksv;}
    else if(nh==ckh&&(it&255)==0&&now()-cksv>1.0)cksave();
  }
  if(nh==0&&(!bestvalid||N<BN)){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;}
  if(!bestvalid){BN=N;memcpy(B,A,(size_t)N*K);}
  fprintf(stderr,"final N=%d\n",BN);
  outbest();
  return 0;
}
