// DOMAIN: k=40..64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#define T 3
#define K 58
#define V 7
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
#define MAXN 4000
#define TLIM 3470.0
#define ALIM 900.0
#define DRIFT 400
#define KICKT 1e9
#define KICKM 2
#define TMAX 1.5
#define TMIN 0.3
#define WMAX 250
#define PLAT_STALL 2.0
#define PLAT_LEVEL 1.0
#define PLAT_INIT 3.0
#define DRAIN_T 60.0
#define DRAIN_STALL 150.0
#define DRAIN_IDLE 6.0
#define CHMAXH 6
typedef unsigned long long u64;

static volatile sig_atomic_t stopf=0;
static void hnd(int s){(void)s;stopf=1;}
static int tri[NTRI][3];
static int ct_tid[K][PC],ct_o1[K][PC],ct_o2[K][PC],ct_mc[K][PC],ct_m1[K][PC],ct_m2[K][PC];
static unsigned char cnt[NTUP];static unsigned char w[NTUP];
typedef struct{int b;unsigned char o1,o2,p,x;}CT;
static CT ct[K][PC];
static const int MC3[3]={V*V,V,1},M13[3]={V,1,V*V},M23[3]={1,V*V,V};
static unsigned char st[(NTUP+3)/4];
static inline int gs(int i){return (st[i>>2]>>((i&3)<<1))&3;}
static inline void ss(int i){int c=cnt[i];if(c>3)c=3;int sh=(i&3)<<1;st[i>>2]=(unsigned char)((st[i>>2]&~(3<<sh))|(c<<sh));}
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
  for(int i=0;i<NTUP;i++){if(!cnt[i])addh(i);ss(i);}
}
static inline long long delta2(int r,int c,int nv,int*dn){
  int ov=A[r][c];long long g=0;int d=0;const unsigned char*x=A[r];const CT*e=ct[c];
  for(int j=0;j<PC;j++){
    int p=e[j].p;int base=e[j].b+M13[p]*x[e[j].o1]+M23[p]*x[e[j].o2];
    int io=base+MC3[p]*ov,in=base+MC3[p]*nv;
    if(gs(io)==1){g-=w[io];d++;}
    if(gs(in)==0){g+=w[in];d--;}
  }
  *dn=d;return g;
}
static inline void apply(int r,int c,int nv){
  int ov=A[r][c];if(ov==nv)return;const unsigned char*x=A[r];const CT*e=ct[c];
  for(int j=0;j<PC;j++){
    int p=e[j].p;int base=e[j].b+M13[p]*x[e[j].o1]+M23[p]*x[e[j].o2];
    int io=base+MC3[p]*ov,in=base+MC3[p]*nv;
    if(--cnt[io]==0)addh(io);
    if(cnt[io]<=3)ss(io);
    if(cnt[in]++==0)remh(in);
    if(cnt[in]<=3)ss(in);
  }
  A[r][c]=nv;
}
static void delrow(int r){
  for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(--cnt[i]==0)addh(i);if(cnt[i]<=3)ss(i);}
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
static int lastu,minu;static int exclr=-1;
static inline int lossinfo(int r,int c,int nv,int*lost){
  int ov=A[r][c];int L=0;(void)nv;const unsigned char*x=A[r];const CT*e=ct[c];
  for(int j=0;j<PC;j++){
    int p=e[j].p;int io=e[j].b+M13[p]*x[e[j].o1]+M23[p]*x[e[j].o2]+MC3[p]*ov;
    if(gs(io)==1){L++;*lost=io;if(L>1)return L;}
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
    if(r==exclr)continue;
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
static void trychains(int comp){
  int hs[CHMAXH+2],m=nh<CHMAXH?nh:CHMAXH;
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
#define RRCAP 32768
#define RRTRY 24
static long rrok=0,rrside=0;
static void rowsub(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(--cnt[i]==0)addh(i);if(cnt[i]<=3)ss(i);}}
static void rowadd(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]++==0)remh(i);if(cnt[i]<=3)ss(i);}}
static int rowresolve(int tries){
  static int H[RRCAP];static int fx[K];static unsigned char best[K],x[K];
  int imp=0;
  for(int a=0;a<tries&&nh>0&&!stopf;a++){
    int nh0=nh;if(nh0>RRCAP/2)return imp;
    int r=rnd()%N;int no=nh0;
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
  int ov=A[r][c];int d=0;const unsigned char*x=A[r];const CT*e=ct[c];
  for(int j=0;j<PC;j++){
    int p=e[j].p;int base=e[j].b+M13[p]*x[e[j].o1]+M23[p]*x[e[j].o2];
    int io=base+MC3[p]*ov,in=base+MC3[p]*nv;
    int so=gs(io);
    if(so==1)return 1;
    if(so==2)d++;
    int sn=gs(in);
    if(sn==1)d--;else if(sn==0)d++;
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
    for(int k=0;k<1024;k++){
      it++;int r=rnd()%N,c=rnd()%K,nv=rnd()%V;if(nv==A[r][c])continue;
      if(guided(r,c,nv))acc++;
    }
  }
  fprintf(stderr,"plateau %.0fs it=%ld acc=%ld nh=%d U=%ld\n",secs,it,acc,nh,countU());
}
static int rowuniq(int r){int s=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1)s++;return s;}
#define DCN 1200
#define DCD 3
static long dchok=0;
/* move unique tuple u (covered only by row r) onto another row, allowing one loss repaired by a chain */
static int drainchain(int u,int r){
  int t=u/VT,v=u%VT;
  int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
  int qr[64],qc[64],qv[64],ql[64],nq=0;
  for(int s=0;s<N&&nq<64;s++){if(s==r)continue;
    unsigned char*x=A[s];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
    if(m1+m2+m3!=2)continue;
    int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
    int lost=-1;int L=lossinfo(s,c,nv,&lost);
    if(L==0){apply(s,c,nv);return 1;}
    if(L==1){qr[nq]=s;qc[nq]=c;qv[nq]=nv;ql[nq]=lost;nq++;}
  }
  for(int i=0;i<nq;i++){int j=i+rnd()%(nq-i);
    int a;a=qr[i];qr[i]=qr[j];qr[j]=a;a=qc[i];qc[i]=qc[j];qc[j]=a;a=qv[i];qv[i]=qv[j];qv[j]=a;a=ql[i];ql[i]=ql[j];ql[j]=a;}
  exclr=r;nodes=NODECAP-DCN;
  for(int i=0;i<nq&&nodes<=NODECAP;i++){
    int s=qr[i],c=qc[i],ov=A[s][c];
    apply(s,c,qv[i]);sd=1;stk_r[0]=s;stk_c[0]=c;
    if(cnt[ql[i]]>0){exclr=-1;return 1;}
    if(chain(ql[i],DCD)){exclr=-1;return 1;}
    sd=0;apply(s,c,ov);
  }
  exclr=-1;return 0;
}
static int drain(int r,double secs){
  double ts=now();double rem=TLIM-(ts-t0);if(secs>rem)secs=rem;
  static int lst[NTRI];int u0=rowuniq(r),bestu=u0;long mv=0;double lastimp=ts;
  while(!stopf&&now()-ts<secs){
    int n=0;for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]==1)lst[n++]=i;}
    if(n<bestu){bestu=n;lastimp=now();}
    if(n==0)break;
    if(now()-lastimp>DRAIN_IDLE)break;
    for(int q=0;q<n&&!stopf;q++){int u=lst[q];if(cnt[u]!=1)continue;
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
    for(int q=0;q<n&&!stopf;q++){int u=lst[q];if(cnt[u]!=1)continue;
      if(drainchain(u,r)){mv++;dchok++;}
    }
    if(nh>0){/* safety: chains must never leave holes */ }
    for(int k=0;k<20000;k++){
      int s=rnd()%N;if(s==r)continue;int c=rnd()%K,nv=rnd()%V;if(nv==A[s][c])continue;
      guided(s,c,nv);
    }
  }
  int uf=rowuniq(r);
  fprintf(stderr,"drain row uniq %d -> %d (min %d) moves %ld chains %ld nh=%d\n",u0,uf,bestu,mv,dchok,nh);
  return uf;
}
static int pickrow(int spread){
  static int u[MAXN];static int ord[MAXN];
  for(int r=0;r<N;r++){u[r]=rowuniq(r);ord[r]=r;}
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
  for(int r=0;r<N;r++)if(uc[r]<=lim&&uc[r]<=UCAP)L[nl++]=r;
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
static void reduce(int spread,double dt){
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
      delrow(s);
      fprintf(stderr,"merge uniq %d+%d -> holes %d (pred %d, best %d)\n",us,ul2,nh,mg_r[k],mg_r[0]);
      return;
    }
  }
  delrow(r);
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
  for(int c=0;c<K;c++)for(int j=0;j<PC;j++){
    ct[c][j].b=ct_tid[c][j]*VT;ct[c][j].o1=(unsigned char)ct_o1[c][j];ct[c][j].o2=(unsigned char)ct_o2[c][j];
    ct[c][j].p=(unsigned char)(ct_mc[c][j]==V*V?0:(ct_mc[c][j]==V?1:2));ct[c][j].x=0;
  }
  for(int i=0;i<NTUP;i++)w[i]=1;
  const char*fn=getenv("BEST_CANDIDATE");N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int x,c=0;while(N<MAXN&&fscanf(f,"%d",&x)==1){A[N][c++]=(unsigned char)(((x%V)+V)%V);if(c==K){c=0;N++;}}fclose(f);}}
  if(N==0){
    /* fallback: random array, greedily extended until valid */
    N=VT*3;if(N>MAXN)N=MAXN;
    for(int r=0;r<N;r++)for(int c=0;c<K;c++)A[r][c]=rnd()%V;
    rebuild();
    while(nh>0&&N<MAXN){
      for(int c=0;c<K;c++)A[N][c]=rnd()%V;
      int h=holes[rnd()%nh];int t=h/VT,v=h%VT;
      A[N][tri[t][0]]=v/(V*V);A[N][tri[t][1]]=(v/V)%V;A[N][tri[t][2]]=v%V;
      N++;rowadd(N-1);
    }
  }
  rebuild();
  BN=N;memcpy(B,A,(size_t)N*K);
  int bestvalid=(nh==0);
  fprintf(stderr,"loaded N=%d holes=%d\n",N,nh);
  if(stopf){outbest();return 0;}
  if(bestvalid){fprintf(stderr,"initial U=%ld\n",countU());plateau(PLAT_INIT);
    if(nh==0){memcpy(B,A,(size_t)N*K);}}
  ckreset();
  long long it=0;double ast=now();int spread=1;int minh=1<<30;double Temp=TMAX;
  static int cr[MAXN],cc[MAXN],cv[MAXN];static long long cdl[MAXN];
  while(!stopf){
    if(nh==0){
      if(!bestvalid||N<BN){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);fprintf(stderr,"PROGRESS %d\n",-N);}
      plateau(PLAT_LEVEL);
      if(stopf)break;
      reduce(spread,DRAIN_T);
      if(stopf)break;
      for(int i=0;i<NTUP;i++)w[i]=1;Wsum=nh;
      ckreset();
      ast=now();minh=nh;spread=1;
      fprintf(stderr,"N=%d holes=%d\n",N,nh);
      rowresolve(3000);
      fprintf(stderr,"after rowresolve holes=%d (rrok=%ld side=%ld)\n",nh,rrok,rrside);
      minh=nh;
      if(nh<=CHMAXH){trychains(1);fprintf(stderr,"after chains holes=%d (chainok=%ld compok=%ld)\n",nh,chainok,compok);minh=nh;}
      if(nh==0)continue;
    }
    it++;
    if((it&255)==0){
      double tn=now();
      if(tn-t0>TLIM)break;
      {double ph=(tn-tbase)/KICKT;if(ph>1.0)ph=1.0;Temp=TMAX*pow(TMIN/TMAX,ph);(void)Temp;}
      if((it&((1<<18)-1))==0)fprintf(stderr,"N=%d nh=%d minh=%d ckh=%d t=%.0f compok=%ld kicks=%ld\n",N,nh,minh,ckh,tn-t0,compok,kicks);
      if((it&((1<<16)-1))==0){for(int i=0;i<NTUP;i++)if(w[i]>1)w[i]=(unsigned char)(1+((w[i]-1)*3)/4);
        Wsum=0;for(int k=0;k<nh;k++)Wsum+=w[holes[k]];}
      if(tn-ckt>KICKT&&nh>0){kick();if(nh==0)continue;}
      if(nh>ckh+3*DRIFT){memcpy(A,Ck,(size_t)N*K);rebuild();if(nh==0)continue;}
      if(tn-ast>ALIM&&bestvalid){
        fprintf(stderr,"stall N=%d minh=%d, restore+plateau\n",N,minh);
        N=BN;memcpy(A,B,(size_t)N*K);rebuild();
        plateau(PLAT_STALL);
        if(stopf)break;
        if(nh==0){memcpy(B,A,(size_t)N*K);}
        spread=4;
        reduce(spread,DRAIN_STALL);
        if(stopf)break;
        for(int i=0;i<NTUP;i++)w[i]=1;Wsum=nh;ckreset();ast=now();minh=nh;spread=1;
        rowresolve(3000);minh=nh;
        if(nh<=CHMAXH){trychains(1);minh=nh;}
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
      trychains((it&16383)==0);if(nh<minh)minh=nh;
      if(nh==0)continue;
      for(int k=0;k<nh;k++){int hh=holes[k];if(w[hh]<WMAX){w[hh]++;Wsum++;}}
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
    if((rnd()&3)==0){
      /* 2-entry moves on rows matching exactly one of the hole's columns */
      int cs3[3]={c1,c2,c3},vs3[3]={v1,v2,v3};
      long long b2=-(1LL<<60);int br=-1,bmi=-1,ties=0,st=rnd()%N,ev=0;
      for(int k=0;k<N;k++){int r=(st+k)%N;unsigned char*x=A[r];
        int m=0,mi=-1;for(int d=0;d<3;d++)if(x[cs3[d]]==vs3[d]){m++;mi=d;}
        if(m!=1)continue;
        if(++ev>48)break;
        int a=(mi+1)%3,b=(mi+2)%3;int oa=x[cs3[a]];long long W0=Wsum;int h0=nh;
        apply(r,cs3[a],vs3[a]);int dn;long long d=delta2(r,cs3[b],vs3[b],&dn);
        long long g=(W0-Wsum)+d;int nn=nh+dn;
        apply(r,cs3[a],oa);(void)h0;
        if(nn>cap)continue;
        if(g>b2){b2=g;br=r;bmi=mi;ties=1;}else if(g==b2&&rnd()%(++ties)==0){br=r;bmi=mi;}
      }
      if(br>=0&&(b2>0||(na>0&&b2>bd)||(b2>-8&&(rnd()%1000)<(unsigned)(20+2*b2>0?20+2*b2:1)))){
        int a=(bmi+1)%3,b=(bmi+2)%3;apply(br,cs3[a],vs3[a]);apply(br,cs3[b],vs3[b]);
        if(nh<minh)minh=nh;
        if(nh<ckh){ckh=nh;cksave();ckt=cksv;}
        continue;
      }
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
        if(bd>-6&&(rnd()%1000)<(unsigned)(10+3*bd>0?10+3*bd:1))apply(cr[bk],cc[bk],cv[bk]);
        else if((rnd()%1000)<3){int k=rnd()%na;apply(cr[k],cc[k],cv[k]);}
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
