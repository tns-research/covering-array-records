// DOMAIN: k=8..24
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
#define ALIM 600.0
#define KICKT 25.0
#define KICKM 2
#define TPER 40.0
#define TMAX 1.5
#define TMIN 0.3
#define WMAX 8
#define PLAT_STALL 30.0
#define PLAT_LEVEL 2.0
#define UC2 64
#define SCAP 48
#define HCAP (2*UC2+NTRI+4)
#define PNODE 400
#define PAIRGAP 6.0
#define PAIRBUD 0.8

static volatile sig_atomic_t stopf=0;
static void hnd(int s){(void)s;stopf=1;}
static int tri[NTRI][3];
static int ct_tid[K][PC],ct_o1[K][PC],ct_o2[K][PC],ct_mc[K][PC],ct_m1[K][PC],ct_m2[K][PC];
static int cnt[NTUP],w[NTUP],hpos[NTUP],holes[NTUP],nh;
static unsigned char A[MAXN][K],B[MAXN][K];
static int N,BN;
static double t0;
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
static int stk_r[CHD+4],stk_c[CHD+4],sd;static long nodes,chainok=0,compok=0;
static int lastu,minu;
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
/* ---- 2->2 exact pair re-solve ---- */
static int ub[MAXN][UC2],ubn[MAXN];
static void buildU(void){
  for(int r=0;r<N;r++){int k=0;for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]==1){if(k<UC2)ub[r][k]=i;k++;}}ubn[r]=k;}
}
static int Hc[HCAP][3],Hv[HCAP][3],nH;
static int X[2][K];static long pnodes,pairok=0,pairtry=0;
static inline void addH(int i){int t=i/VT,v=i%VT;
  Hc[nH][0]=tri[t][0];Hc[nH][1]=tri[t][1];Hc[nH][2]=tri[t][2];
  Hv[nH][0]=v/(V*V);Hv[nH][1]=(v/V)%V;Hv[nH][2]=v%V;nH++;}
static int csp(void){
  if(++pnodes>PNODE)return 0;
  int best=-1,bo=3,bok[2]={0,0};
  for(int i=0;i<nH;i++){
    int ok[2]={0,0},s=0;
    for(int z=0;z<2&&!s;z++){int sat=1,can=1;
      for(int d=0;d<3;d++){int xv=X[z][Hc[i][d]];if(xv!=Hv[i][d]){sat=0;if(xv>=0){can=0;break;}}}
      if(sat)s=1;ok[z]=can;}
    if(s)continue;
    int no=ok[0]+ok[1];
    if(no==0)return 0;
    if(no<bo){bo=no;best=i;bok[0]=ok[0];bok[1]=ok[1];if(no==1)break;}
  }
  if(best<0)return 1;
  int z0=rnd()&1;
  for(int zz=0;zz<2;zz++){int z=zz^z0;if(!bok[z])continue;
    int sv[3];for(int d=0;d<3;d++){sv[d]=X[z][Hc[best][d]];X[z][Hc[best][d]]=Hv[best][d];}
    if(csp())return 1;
    for(int d=2;d>=0;d--)X[z][Hc[best][d]]=sv[d];
    if(pnodes>PNODE)return 0;
  }
  return 0;
}
static int pairsolve(int u,double budget){
  if(cnt[u]>0)return 1;
  double ts=now();
  buildU();
  int t=u/VT,v=u%VT;int cs[3]={tri[t][0],tri[t][1],tri[t][2]},vs[3]={v/(V*V),(v/V)%V,v%V};
  static int R[MAXN],S[MAXN];int nr=0,ns=0;
  for(int m=2;m>=0;m--){int st=nr;
    for(int r=0;r<N;r++){if(ubn[r]>UC2)continue;int mm=0;for(int d=0;d<3;d++)mm+=A[r][cs[d]]==vs[d];if(mm==m)R[nr++]=r;}
    for(int i=st;i<nr;i++){int j=i+rnd()%(nr-i);int x=R[i];R[i]=R[j];R[j]=x;}
  }
  for(int r=0;r<N;r++)if(ubn[r]<=SCAP)S[ns++]=r;
  for(int i=0;i<ns;i++){int j=i+rnd()%(ns-i);int x=S[i];S[i]=S[j];S[j]=x;}
  for(int a=0;a<nr&&!stopf;a++){
    int r=R[a];if((a&3)==0&&now()-ts>budget)break;
    for(int b=0;b<ns;b++){int s=S[b];if(s==r)continue;
      nH=0;
      for(int k=0;k<ubn[r];k++)addH(ub[r][k]);
      for(int k=0;k<ubn[s];k++)addH(ub[s][k]);
      int msk=0;for(int c=0;c<K;c++)if(A[r][c]==A[s][c])msk|=1<<c;
      if(__builtin_popcount(msk)>=3)for(int tt=0;tt<NTRI;tt++){
        if(((msk>>tri[tt][0])&1)&&((msk>>tri[tt][1])&1)&&((msk>>tri[tt][2])&1)){int i=tidx(tt,A[r]);if(cnt[i]==2)addH(i);}
      }
      for(int c=0;c<K;c++){X[0][c]=-1;X[1][c]=-1;}
      for(int d=0;d<3;d++)X[0][cs[d]]=vs[d];
      pnodes=0;pairtry++;
      if(csp()){
        int nx[K],ny[K];
        for(int c=0;c<K;c++){nx[c]=X[0][c]>=0?X[0][c]:A[r][c];ny[c]=X[1][c]>=0?X[1][c]:A[s][c];}
        int h0=nh;
        for(int c=0;c<K;c++)apply(r,c,nx[c]);
        for(int c=0;c<K;c++)apply(s,c,ny[c]);
        pairok++;
        fprintf(stderr,"pairsolve ok: holes %d -> %d (tries %ld, t=%.2fs)\n",h0,nh,pairtry,now()-ts);
        return cnt[u]>0;
      }
    }
  }
  return 0;
}
static void pairall(double bud){
  if(bud<1e9)return; /* disabled: 0 successes in 579k tries */
  int hs[CHMAXH+2],m=nh<CHMAXH?nh:CHMAXH;
  for(int i=0;i<m;i++)hs[i]=holes[i];
  for(int i=0;i<m&&!stopf;i++){if(cnt[hs[i]]>0)continue;if(now()-t0>TLIM)break;pairsolve(hs[i],bud);}
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
static unsigned char Ck[MAXN][K];static int ckh;static double ckt,cksv,tbase;static long kicks=0;
static void cksave(void){memcpy(Ck,A,(size_t)N*K);cksv=now();}
static void ckreset(void){ckh=nh;cksave();ckt=tbase=cksv;}
static void kick(void){
  memcpy(A,Ck,(size_t)N*K);rebuild();
  for(int k=0;k<nh;k++){int hh=holes[k];if(w[hh]<WMAX)w[hh]++;}
  for(int m=0;m<KICKM;m++){int r=rnd()%N,c=rnd()%K,nv=rnd()%V;apply(r,c,nv);}
  kicks++;ckt=tbase=now();
}
static void plateau(double secs){
  double ts=now();long acc=0,it=0;
  double rem=TLIM-(ts-t0);if(secs>rem)secs=rem;
  while(!stopf&&now()-ts<secs){
    for(int k=0;k<4096;k++){
      it++;int r=rnd()%N,c=rnd()%K,nv=rnd()%V;if(nv==A[r][c])continue;
      int lost;if(lossinfo(r,c,nv,&lost)==0){apply(r,c,nv);acc++;}
    }
  }
  fprintf(stderr,"plateau %.0fs it=%ld acc=%ld nh=%d\n",secs,it,acc,nh);
}
static int pickrow(int spread){
  static int u[MAXN];static int ord[MAXN];
  for(int r=0;r<N;r++){int s=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1)s++;u[r]=s;ord[r]=r;}
  for(int i=0;i<N;i++){int j=i+rnd()%(N-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  int m=spread<N?spread:N;
  for(int i=0;i<m;i++){int b=i;for(int j=i+1;j<N;j++)if(u[ord[j]]<u[ord[b]])b=j;int x=ord[i];ord[i]=ord[b];ord[b]=x;}
  int k=rnd()%m;
  lastu=u[ord[k]];minu=u[ord[0]];
  fprintf(stderr,"delete row uniq=%d (min %d)\n",u[ord[k]],u[ord[0]]);
  return ord[k];
}
#define UCAP 20
#define MLIST 400
#define PMAX 12
#define MTOP 64
static int ul[MAXN][UCAP],uc[MAXN];
static int mg_s[MTOP],mg_l[MTOP],mg_m[MTOP],mg_r[MTOP],nmg;
static int cmpuc(const void*a,const void*b){return uc[*(const int*)a]-uc[*(const int*)b];}
static int mergecand(int lim){
  static int L[MAXN];int nl=0;
  for(int r=0;r<N;r++){int k=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1){if(k<UCAP)ul[r][k]=t;k++;}uc[r]=k;}
  for(int r=0;r<N;r++)if(uc[r]<=lim&&uc[r]<=UCAP)L[nl++]=r;
  for(int i=0;i<nl;i++){int j=i+rnd()%(nl-i);int x=L[i];L[i]=L[j];L[j]=x;}
  qsort(L,nl,sizeof(int),cmpuc);
  if(nl>MLIST)nl=MLIST;
  nmg=0;
  static int um[1<<PMAX];
  for(int i=0;i<nl;i++)for(int j=i+1;j<nl;j++){
    int s=L[i],l=L[j];if(uc[s]>uc[l]){int x=s;s=l;l=x;}
    int p=uc[s],q=uc[l];if(p>PMAX)continue;
    int cs[UCAP],cl[UCAP];
    for(int k=0;k<p;k++){int t=ul[s][k],m=0;for(int d=0;d<3;d++){int c=tri[t][d];if(A[s][c]!=A[l][c])m|=1<<c;}cs[k]=m;}
    for(int k=0;k<q;k++){int t=ul[l][k],m=0;for(int d=0;d<3;d++){int c=tri[t][d];if(A[s][c]!=A[l][c])m|=1<<c;}cl[k]=m;}
    int bs=-1,bm=0;um[0]=0;
    for(int S=0;S<(1<<p);S++){
      if(S)um[S]=um[S&(S-1)]|cs[__builtin_ctz(S)];
      int uu=um[S],sat=__builtin_popcount(S);
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
static void reduce(int spread){
  int r=pickrow(spread);
  int n=mergecand(minu+8);
  if(n>0){
    int k=rnd()%(spread<n?spread:n);
    if(mg_r[k]<lastu||(mg_r[k]==lastu&&(rnd()&1))){
      int s=mg_s[k],l=mg_l[k],m=mg_m[k];
      int us=uc[s],ul2=uc[l];
      for(int c=0;c<K;c++){int nv=((m>>c)&1)?A[s][c]:A[l][c];apply(l,c,nv);}
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
  const char*fn=getenv("BEST_CANDIDATE");N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int x,c=0;while(N<MAXN&&fscanf(f,"%d",&x)==1){A[N][c++]=(unsigned char)(((x%V)+V)%V);if(c==K){c=0;N++;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 0;}
  rebuild();
  BN=N;memcpy(B,A,(size_t)N*K);
  int bestvalid=(nh==0);
  fprintf(stderr,"loaded N=%d holes=%d\n",N,nh);
  if(!bestvalid)BN=0;
  for(int i=0;i<NTUP;i++)w[i]=1;
  ckreset();
  long long it=0;double ast=now();int spread=1;int minh=1<<30;double Temp=TMAX;double lastpair=0;
  static int cr[MAXN],cc[MAXN],cv[MAXN];static long long cdl[MAXN];static double cw[MAXN];
  while(!stopf){
    if(nh==0){
      if(!bestvalid||N<BN){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;fprintf(stderr,"valid N=%d t=%.0f pairok=%ld\n",N,now()-t0,pairok);}
      plateau(PLAT_LEVEL);
      if(stopf)break;
      reduce(spread);
      for(int i=0;i<NTUP;i++)w[i]=1;
      ckreset();
      ast=now();minh=nh;spread=1;
      fprintf(stderr,"N=%d holes=%d\n",N,nh);
      if(nh<=CHMAXH){trychains(1);fprintf(stderr,"after chains holes=%d (chainok=%ld compok=%ld)\n",nh,chainok,compok);
        if(nh>0){pairall(PAIRBUD);lastpair=now();fprintf(stderr,"after pairs holes=%d\n",nh);}
        minh=nh;}
      if(nh==0)continue;
    }
    it++;
    if((it&1023)==0){
      double tn=now();
      if(tn-t0>TLIM)break;
      {double ph=(tn-tbase)/KICKT;if(ph>1.0)ph=1.0;Temp=TMAX*pow(TMIN/TMAX,ph);}
      if((it&((1<<20)-1))==0)fprintf(stderr,"N=%d nh=%d minh=%d ckh=%d T=%.2f t=%.0f compok=%ld kicks=%ld\n",N,nh,minh,ckh,Temp,tn-t0,compok,kicks);
      if(tn-ckt>KICKT&&nh>0){kick();if(nh==0)continue;}
      if(tn-ast>ALIM&&bestvalid){
        fprintf(stderr,"stall N=%d minh=%d, restore+plateau\n",N,minh);
        N=BN;memcpy(A,B,(size_t)N*K);rebuild();
        plateau(PLAT_STALL);
        if(stopf)break;
        if(nh==0){memcpy(B,A,(size_t)N*K);}
        spread=20;
        reduce(spread);
        for(int i=0;i<NTUP;i++)w[i]=1;ckreset();ast=now();minh=nh;spread=1;
        if(nh<=CHMAXH){trychains(1);if(nh>0){pairall(PAIRBUD);lastpair=now();}minh=nh;}
        continue;
      }
    }
    if((it&4095)==0&&nh<=CHMAXH){
      trychains((it&16383)==0);if(nh<minh)minh=nh;
      if(nh==0)continue;
      if(nh<=3){double tn=now();if(tn-lastpair>PAIRGAP){pairall(PAIRBUD);lastpair=now();if(nh<minh)minh=nh;if(nh==0)continue;}}
      for(int k=0;k<nh;k++){int hh=holes[k];if(w[hh]<WMAX)w[hh]++;}
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
    if(nh<ckh){ckh=nh;cksave();ckt=cksv;}
    else if(nh==ckh&&(it&255)==0&&now()-cksv>1.0)cksave();
  }
  if(nh==0&&(!bestvalid||N<BN)){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;}
  if(!bestvalid){BN=N;memcpy(B,A,(size_t)N*K);}
  fprintf(stderr,"final N=%d pairok=%ld pairtry=%ld\n",BN,pairok,pairtry);
  outbest();
  return 0;
}
