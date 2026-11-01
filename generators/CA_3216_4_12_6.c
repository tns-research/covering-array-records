#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <signal.h>
#define T 4
#define K 12
#define V 6
#define MUL 3 /* column multiplier on Z_(K-1); 3 has order 5 mod 11 */
#define NT (T==2?V*V:T==3?V*V*V:T==4?V*V*V*V:T==5?V*V*V*V*V:V*V*V*V*V*V)
#define NS (T==2?K*(K-1)/2:T==3?K*(K-1)*(K-2)/6:T==4?K*(K-1)*(K-2)*(K-3)/24:T==5?K*(K-1)*(K-2)*(K-3)*(K-4)/120:K*(K-1)*(K-2)*(K-3)*(K-4)*(K-5)/720)
#define CK1 (NS*T/K)
#define NG (NS*NT)
#define CAP 5000
#define HOLE 1000
#define TAILS 10
#define TAILT 25.0
#define GMAX 4096
#define BW (CAP/64+1)
typedef unsigned long long ull;
static int m,nO,N0,bestN=INT_MAX,rawmode=0;
static unsigned char A[CAP][K],B[CAP][K],I0[CAP][K];
static unsigned short cnt[NG],xr[NG],xs[NG];
static int orb[NG],uf[NG],memb[NG],mst[NG+1],mis[NG],pos[NG],M=0;
static int vc[NG];static unsigned vst[NG],ep=0;
static long long tb[CAP][K],fit=0;
static short brk[CAP][K];
static ull bset[K][V][BW];
static int scol[NS][T],sid[1<<K],wts[T];
static int csg[K][CK1],cot[K][CK1][T-1],cw[K][CK1],cwo[K][CK1][T-1],csi[K][CK1];
static int remc[65536],addc[65536];
static int gp[GMAX][K],gs[GMAX][V],ng;
static int genp[6][K],gens[6][V],ngen;
static int F=0; /* rows [0,F) have short orbits: frozen, never moved or deleted */
static const int SYA[6]={1,2,3,4,0,5},SYB[6]={0,2,4,1,3,5},SYBB[6]={0,4,3,2,1,5},SYC[6]={5,4,2,3,1,0};
static double expt[4096],slk=0,Tm=60,Th[2]={450,450},t0g,TL=3480;
static ull rs=0xD1B54A32D192ED03ULL;
static volatile sig_atomic_t stopf=0;
static void onterm(int x){(void)x;stopf=1;}
static inline ull rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(const unsigned char*x,int s){int t=0;for(int i=0;i<T;i++)t+=x[scol[s][i]]*wts[i];return t;}
static inline void addb(int r,int s,int d){short*b=brk[r];int*sc=scol[s];for(int i=0;i<T;i++)b[sc[i]]+=d;}
static void setsoft(int s){
  for(int i=0;i<8;i++)remc[i]=addc[i]=0;
  slk=5.0*s;
  remc[1]=addc[0]=HOLE;
  remc[2]=addc[1]=s;
  remc[3]=addc[2]=s/3;
  remc[4]=addc[3]=s/9;
}
static inline int rowbase(const unsigned char*x,int c,int j){
  int g=csg[c][j];
  for(int i=0;i<T-1;i++)g+=x[cot[c][j][i]]*cwo[c][j][i];
  return g;
}
static inline int delta_raw(int r,int c,int a,int b){
  unsigned char*x=A[r];int d=0;
  for(int j=0;j<CK1;j++){
    int g=rowbase(x,c,j);
    d+=remc[cnt[g+a*cw[c][j]]]-addc[cnt[g+b*cw[c][j]]];
  }
  return d;
}
/* exact under orbit collisions: several tuples of one move may share an orbit */
static int delta_x(int r,int c,int a,int b){
  unsigned char*x=A[r];int d=0;
  ep++;
  for(int j=0;j<CK1;j++){
    int g=rowbase(x,c,j);
    int o=orb[g+a*cw[c][j]],n=orb[g+b*cw[c][j]];
    int q=(vst[o]==ep)?vc[o]:cnt[o];
    d+=remc[q];vst[o]=ep;vc[o]=q-1;
    q=(vst[n]==ep)?vc[n]:cnt[n];
    d-=addc[q];vst[n]=ep;vc[n]=q+1;
  }
  return d;
}
static inline int delta(int r,int c,int a,int b){return rawmode?delta_raw(r,c,a,b):delta_x(r,c,a,b);}
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  unsigned char*x=A[r];
  for(int j=0;j<CK1;j++){
    int g=rowbase(x,c,j);
    int o=orb[g+a*cw[c][j]],n=orb[g+b*cw[c][j]],s=csi[c][j];
    xr[o]^=(unsigned short)r;xs[o]^=(unsigned short)s;
    int q=--cnt[o];
    if(q==0){addm(o);addb(r,s,-1);}
    else if(q==1)addb(xr[o],xs[o],1);
    q=cnt[n]++;
    if(q==0){delm(n);addb(r,s,1);}
    else if(q==1)addb(xr[n],xs[n],-1);
    xr[n]^=(unsigned short)r;xs[n]^=(unsigned short)s;
  }
  bset[c][a][r>>6]&=~(1ULL<<(r&63));
  bset[c][b][r>>6]|=1ULL<<(r&63);
  A[r][c]=(unsigned char)b;
}
static void rebuild(void){
  memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);memset(xs,0,sizeof xs);
  memset(brk,0,sizeof brk);memset(bset,0,sizeof bset);M=0;
  for(int r=0;r<m;r++){
    for(int s=0;s<NS;s++){int o=orb[s*NT+tidx(A[r],s)];cnt[o]++;xr[o]^=(unsigned short)r;xs[o]^=(unsigned short)s;}
    for(int c=0;c<K;c++)bset[c][A[r][c]][r>>6]|=1ULL<<(r&63);
  }
  for(int o=0;o<nO;o++){
    pos[o]=-1;
    if(cnt[o]==0)addm(o);
    else if(cnt[o]==1)addb(xr[o],xs[o],1);
  }
}
static void delweak(void){
  int br=-1,bu=1<<30,ties=0;
  for(int r=F;r<m;r++){
    int u=0;
    for(int s=0;s<NS;s++)if(cnt[orb[s*NT+tidx(A[r],s)]]==1)u++;
    if(u<bu){bu=u;br=r;ties=1;}
    else if(u==bu){ties++;if(rnd()%ties==0)br=r;}
  }
  m--;
  if(br!=m){memcpy(A[br],A[m],K);memcpy(tb[br],tb[m],sizeof tb[0]);}
  rebuild();
  fprintf(stderr,"deleted orbit uq=%d, m=%d rows<=%d M=%d t=%.0f\n",bu,m,m*ng,M,now()-t0g);
}
/* distance-1 rows for every member of a random hole orbit, noisy Boltzmann choice, short tabu */
static void focus(void){
  int O=mis[rnd()%M],br=-1,bc=-1,bv=-1;
  double bs=0;
  int W=(m+63)/64;
  fit++;
  for(int mi=mst[O];mi<mst[O+1];mi++){
    int g=memb[mi],s=g/NT,id=g%NT,tv[T];
    for(int i=T-1;i>=0;i--){tv[i]=id%V;id/=V;}
    int*sc=scol[s];
    for(int w=0;w<W;w++){
      ull ok=~0ULL,one=0;
      for(int i=0;i<T;i++){ull p=bset[sc[i]][tv[i]][w];one=(one&p)|(ok&~p);ok&=p;}
      while(one){
        int r=w*64+__builtin_ctzll(one);
        one&=one-1;
        if(r<F)continue;
        unsigned char*x=A[r];
        int mc=0;
        while(x[sc[mc]]==tv[mc])mc++;
        int c=sc[mc];
        double th=Tm*expt[rnd()>>52];
        if((brk[r][c]-1)*HOLE-slk-th>=bs)continue;
        int d=delta(r,c,x[c],tv[mc]);
        if(tb[r][c]>fit&&d>-HOLE/2)continue;
        double q=d-th;
        if(q<bs){bs=q;br=r;bc=c;bv=tv[mc];}
      }
    }
  }
  if(br>=0){
    change(br,bc,bv);
    tb[br][bc]=fit+10+(long long)(rnd()%20);
  }
}
static int find(int x){while(uf[x]!=x){uf[x]=uf[uf[x]];x=uf[x];}return x;}
static void addsym(const int*p){
  for(int c=0;c<K;c++)genp[ngen][c]=c;
  for(int v=0;v<V;v++)gens[ngen][v]=p[v];
  ngen++;
}
/* sym: 0 none, 1/2 additive shift, 3 PGL(2,5)=S5, 4 PSL(2,5)=A5, 5 D5, 6 Z5 (3..6 need V==6) */
static int build_group(int rot,int mul,int sym){
  ngen=0;
  if(rot){
    for(int c=0;c<K-1;c++)genp[ngen][c]=(c+1)%(K-1);
    genp[ngen][K-1]=K-1;
    for(int v=0;v<V;v++)gens[ngen][v]=v;
    ngen++;
  }
  if(mul){
    int seen=0;
    for(int c=0;c<K-1;c++)seen|=1<<(c*MUL%(K-1));
    if(seen!=(1<<(K-1))-1)return 0;
    for(int c=0;c<K-1;c++)genp[ngen][c]=c*MUL%(K-1);
    genp[ngen][K-1]=K-1;
    for(int v=0;v<V;v++)gens[ngen][v]=v;
    ngen++;
  }
  if(sym==1||sym==2){int p[V];for(int v=0;v<V;v++)p[v]=(v+sym)%V;addsym(p);}
  if(sym==3){addsym(SYA);addsym(SYB);addsym(SYC);}
  if(sym==4){addsym(SYA);addsym(SYBB);addsym(SYC);}
  if(sym==5){addsym(SYA);addsym(SYBB);}
  if(sym==6){addsym(SYA);}
  ng=1;for(int c=0;c<K;c++)gp[0][c]=c;
  for(int v=0;v<V;v++)gs[0][v]=v;
  for(int i=0;i<ng;i++)for(int h=0;h<ngen;h++){
    int np[K],ns[V],f=0;
    for(int c=0;c<K;c++)np[c]=genp[h][gp[i][c]];
    for(int v=0;v<V;v++)ns[v]=gens[h][gs[i][v]];
    for(int j=0;j<ng&&!f;j++)if(memcmp(gs[j],ns,sizeof ns)==0&&memcmp(gp[j],np,sizeof np)==0)f=1;
    if(!f){
      if(ng==GMAX)return 0;
      memcpy(gp[ng],np,sizeof np);memcpy(gs[ng],ns,sizeof ns);ng++;
    }
  }
  for(int g=0;g<NG;g++)uf[g]=g;
  for(int g=0;g<NG;g++){
    int s=g/NT,id=g%NT,tv[T];
    for(int i=T-1;i>=0;i--){tv[i]=id%V;id/=V;}
    for(int h=0;h<ngen;h++){
      int tmp[K],mask=0;
      for(int i=0;i<T;i++){int c=genp[h][scol[s][i]];tmp[c]=gens[h][tv[i]];mask|=1<<c;}
      int s2=sid[mask],idx=0;
      for(int i=0;i<T;i++)idx+=tmp[scol[s2][i]]*wts[i];
      int a=find(g),b=find(s2*NT+idx);
      if(a!=b)uf[a]=b;
    }
  }
  nO=0;
  static int lab[NG];
  for(int g=0;g<NG;g++)lab[g]=-1;
  for(int g=0;g<NG;g++){int r=find(g);if(lab[r]<0)lab[r]=nO++;orb[g]=lab[r];}
  memset(mst,0,sizeof mst);
  for(int g=0;g<NG;g++)mst[orb[g]+1]++;
  for(int o=0;o<nO;o++)mst[o+1]+=mst[o];
  static int fill[NG];
  for(int o=0;o<nO;o++)fill[o]=mst[o];
  for(int g=0;g<NG;g++)memb[fill[orb[g]]++]=g;
  rawmode=(ng==1);
  return 1;
}
static inline void imgrow(int e,const unsigned char*x,unsigned char*y){for(int c=0;c<K;c++)y[gp[e][c]]=(unsigned char)gs[e][x[c]];}
static inline ull rowcode(const unsigned char*x){ull v=0;for(int c=0;c<K;c++)v=v*V+x[c];return v;}
static int cmpull(const void*a,const void*b){ull x=*(const ull*)a,y=*(const ull*)b;return x<y?-1:x>y;}
static int verify(unsigned char (*R)[K],int n){
  static unsigned char seen[NG];
  memset(seen,0,sizeof seen);
  for(int r=0;r<n;r++)for(int s=0;s<NS;s++)seen[s*NT+tidx(R[r],s)]=1;
  for(int i=0;i<NG;i++)if(!seen[i])return 0;
  return 1;
}
/* expand the orbit representatives, dedupe, verify; returns 1 on a new best */
static int record(void){
  long long tot=(long long)m*ng;
  ull*E=malloc((size_t)tot*sizeof(ull));
  long long n=0;
  for(int r=0;r<m;r++)for(int e=0;e<ng;e++){unsigned char y[K];imgrow(e,A[r],y);E[n++]=rowcode(y);}
  qsort(E,(size_t)n,sizeof(ull),cmpull);
  long long u=0;
  for(long long i=0;i<n;i++)if(i==0||E[i]!=E[i-1])E[u++]=E[i];
  int ret=0;
  if(u<bestN&&u<=CAP){
    static unsigned char R[CAP][K];
    for(long long i=0;i<u;i++){ull v=E[i];for(int c=K-1;c>=0;c--){R[i][c]=(unsigned char)(v%V);v/=V;}}
    if(verify(R,(int)u)){
      bestN=(int)u;memcpy(B,R,(size_t)u*K);ret=1;
      fprintf(stderr,"valid N=%d (order %d, m=%d) t=%.0f\n",bestN,ng,m,now()-t0g);
    } else fprintf(stderr,"BUG: expansion failed verification\n");
  }
  free(E);
  return ret;
}
static void load_reps(unsigned char (*src)[K],int n){
  m=0;
  for(int r=0;r<n;r++){
    ull c0=rowcode(src[r]);
    int mn=1;
    for(int e=1;e<ng&&mn;e++){unsigned char y[K];imgrow(e,src[r],y);if(rowcode(y)<c0)mn=0;}
    if(mn)memcpy(A[m++],src[r],K);
  }
}
/* number of distinct rows after expanding the current representatives */
static long long expsize(void){
  long long tot=(long long)m*ng,n=0,u=0;
  ull*E=malloc((size_t)tot*sizeof(ull));
  for(int r=0;r<m;r++)for(int e=0;e<ng;e++){unsigned char y[K];imgrow(e,A[r],y);E[n++]=rowcode(y);}
  qsort(E,(size_t)n,sizeof(ull),cmpull);
  for(long long i=0;i<n;i++)if(i==0||E[i]!=E[i-1])u++;
  free(E);
  return u;
}
/* move short-orbit representatives to the front and freeze them */
static void setfrozen(void){
  F=0;
  if(ng==1)return;
  ull*E=malloc((size_t)ng*sizeof(ull));
  for(int r=0;r<m;r++){
    for(int e=0;e<ng;e++){unsigned char y[K];imgrow(e,A[r],y);E[e]=rowcode(y);}
    qsort(E,(size_t)ng,sizeof(ull),cmpull);
    int u=0;
    for(int i=0;i<ng;i++)if(i==0||E[i]!=E[i-1])u++;
    if(u<ng){
      unsigned char t[K];
      memcpy(t,A[r],K);memcpy(A[r],A[F],K);memcpy(A[F],t,K);
      F++;
    }
  }
  free(E);
}
/* constant row + non-constant Paley-Hadamard rows (two symbols) + random full rows */
static void seed_s5(void){
  memset(A[0],0,K);
  m=1;
  int p=K-1,prime=(p>2&&p%4==3);
  for(int d=2;d*d<=p;d++)if(p%d==0)prime=0;
  if(prime){
    char qr[K];
    memset(qr,0,sizeof qr);
    for(int x=1;x<p;x++)qr[x*x%p]=1;
    for(int i=0;i<p;i++){
      for(int j=0;j<p;j++)A[m][j]=(unsigned char)((j==i||qr[(j-i+p)%p])?0:1);
      A[m][K-1]=1;
      m++;
    }
  }
  int base=V+(m-1)*V*(V-1);
  int add=(N0-base)/ng+3;
  if(add<1)add=1;
  for(int r=0;r<add;r++,m++)for(int c=0;c<K;c++)A[m][c]=(unsigned char)(rnd()%V);
}
typedef struct{int rot,mul,sym;double fr;}Stage;
static const Stage CHA[5]={{0,0,3,0.18},{0,0,4,0.20},{0,0,5,0.40},{0,0,6,0.50},{0,0,0,1.0}};
static const Stage CHB[5]={{1,1,1,0.10},{1,0,1,0.35},{1,0,2,0.35},{1,0,0,0.60},{0,0,0,1.0}};
/* first stage of chain B whose group leaves the input array invariant */
static int firstinv(void){
  for(int i=0;i<4;i++){
    if(!build_group(CHB[i].rot,CHB[i].mul,CHB[i].sym))continue;
    load_reps(I0,N0);
    if(expsize()==N0)return i;
  }
  return 0;
}
static void run_stage(double tend,double pat){
  double scale=m*K/47000.0,D0=40*scale;
  if(D0<0.25)D0=0.25;
  double Dmax=900*scale;if(Dmax<120)Dmax=120;
  double D=D0,Dc=D,Tc0=Th[0],Tc1=0.65*Th[0],Mt=3,sumM=0;
  long long nM=0,it=0,accU=0;
  int phase=0,cycle=-1,focE=16,bestMs=M,succ=0,si=0;
  double tn=now(),phStart=tn,tlog=tn,lastImp=tn;
  setsoft(TAILS);Tm=TAILT;
  while(1){
    if((it&4095)==0){
      tn=now();
      if(stopf||tn-t0g>TL||tn>tend||tn-lastImp>pat)break;
      double el=tn-phStart;
      if(tn-tlog>60){
        tlog=tn;
        fprintf(stderr,"t=%.0f order=%d cyc=%d ph=%d m=%d M=%d bestMs=%d T=%.1f Th=%.0f/%.0f D=%.1f Mt=%.0f bestN=%d accU=%lld it=%lld\n",
          tn-t0g,ng,cycle,phase,m,M,bestMs,Tm,Th[0],Th[1],D,Mt,bestN,accU,it);
      }
      if(phase==0){
        double dur=D*0.2*(M<=2?3:1);
        if(el>dur){
          cycle++;
          si=(cycle%4==3);
          setsoft(si?15:30);
          Mt=(si?1.0e-4:1.5e-4)*nO;if(Mt<3)Mt=3;
          Tm=Th[si];focE=128;phase=1;phStart=tn;sumM=0;nM=0;
        }
      } else if(phase==1){
        if(el>D/4){sumM+=M;nM++;}
        if(el>D/2){
          double avg=nM?sumM/nM:M;
          double f=avg>2*Mt?0.9:avg>Mt?0.97:avg<Mt/2?1.1:1.03;
          Th[si]*=f;
          if(Th[si]<60)Th[si]=60;
          if(Th[si]>3000)Th[si]=3000;
          Tm=Tc0=Th[si];Tc1=0.65*Tc0;Dc=D;
          phase=2;phStart=tn;focE=64;succ=0;
        }
      } else {
        double f=el/Dc;
        if(f>=1){
          if(!succ){D*=1.2;if(D>Dmax)D=Dmax;}
          phase=0;phStart=tn;focE=16;
          setsoft(TAILS);Tm=TAILT;
        } else Tm=Tc0*pow(Tc1/Tc0,f);
      }
    }
    it++;
    if(M==0){
      if(record())lastImp=now();
      if(m<=F+1)break;
      delweak();
      bestMs=M;
      tn=now();
      if(phase==2){
        /* chain: keep descending on the next orbit instead of restarting the cycle */
        succ=1;
        D*=0.9;if(D<0.4*D0)D=0.4*D0;
        Tc0=Tm*1.02;if(Tc0>Th[si])Tc0=Th[si];
        if(Tc0<Tc1*1.05)Tc0=Tc1*1.05;
        Dc=D*log(Tc0/Tc1)/log(1/0.65);if(Dc<0.05)Dc=0.05;
        phStart=tn;
      } else if(phase==0)phStart=tn;
      continue;
    }
    if(M<bestMs)bestMs=M;
    if((it&(focE-1))==0){focus();continue;}
    {
      ull u=rnd();
      int r=F+(int)(((u&0xFFFFFULL)*(ull)(m-F))>>20);
      int c=(int)((((u>>20)&0xFFFFULL)*K)>>16);
      double th=Tm*expt[u>>52];
      if(brk[r][c]*HOLE>th+slk)continue;
      int a=A[r][c],b=(a+1+(int)((((u>>36)&0xFFFFULL)*(V-1))>>16))%V;
      int d=delta(r,c,a,b);
      if((double)d<=th){change(r,c,b);accU++;}
    }
  }
  if(M==0)record();
  fprintf(stderr,"stage order=%d end t=%.0f m=%d M=%d bestMs=%d bestN=%d it=%lld\n",ng,now()-t0g,m,M,bestMs,bestN,it);
}
int main(void){
  t0g=now();
  signal(SIGTERM,onterm);
  rs^=(ull)time(NULL)*0x9E3779B97F4A7C15ULL;if(!rs)rs=0xD1B54A32D192ED03ULL;
  for(int i=0;i<20;i++)rnd();
  const char*fn=getenv("BEST_CANDIDATE");
  N0=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){I0[N0][c++]=(unsigned char)v;if(c==K){c=0;N0++;if(N0>=CAP)break;}}fclose(f);}}
  if(N0==0){fprintf(stderr,"no input\n");return 1;}
  {int w=1;for(int i=T-1;i>=0;i--){wts[i]=w;w*=V;}}
  int s=0;
  for(int mask=0;mask<(1<<K);mask++)if(__builtin_popcount(mask)==T){
    sid[mask]=s;int q=0;
    for(int c=0;c<K;c++)if(mask>>c&1)scol[s][q++]=c;
    s++;
  }
  int cc[K];memset(cc,0,sizeof cc);
  for(s=0;s<NS;s++)for(int i=0;i<T;i++){
    int c=scol[s][i],j=cc[c]++,q=0;
    csg[c][j]=s*NT;cw[c][j]=wts[i];csi[c][j]=s;
    for(int k=0;k<T;k++)if(k!=i){cot[c][j][q]=scol[s][k];cwo[c][j][q]=wts[k];q++;}
  }
  for(int i=0;i<4096;i++)expt[i]=-log((i+0.5)/4096.0);
  /* chain A: S5 > A5 > D5 > Z5 > 1 on symbols; chain B: the parent's Z_11 x Z_6 chain.
     Every stage's array stays invariant under the next stage's group. */
  int chain=(V==6)?0:1;
  int st=chain?firstinv():0;
  while(st<5){
    double tn=now(),rem=TL-(tn-t0g);
    if(stopf||rem<=0)break;
    Stage S=chain?CHB[st]:CHA[st];
    if(!build_group(S.rot,S.mul,S.sym)){fprintf(stderr,"stage %d skipped\n",st);st++;continue;}
    if(chain==0&&st==0){
      seed_s5();
    } else if(ng==1){
      if(bestN<N0)load_reps(B,bestN);
      else load_reps(I0,N0);
    } else if(bestN<INT_MAX){
      load_reps(B,bestN);
    } else {
      load_reps(I0,N0);
      if(expsize()!=N0){
        m=N0/ng+1;
        for(int r=0;r<m;r++)for(int c=0;c<K;c++)A[r][c]=(unsigned char)(rnd()%V);
      }
    }
    setfrozen();
    memset(tb,0,sizeof tb);fit=0;
    rebuild();
    double bud=rem*S.fr;
    fprintf(stderr,"chain %d stage %d order=%d orbits=%d m=%d F=%d M=%d budget=%.0f t=%.0f\n",chain,st,ng,nO,m,F,M,bud,tn-t0g);
    run_stage(tn+bud,ng==1?1e18:bud*0.5);
    st++;
    if(chain==0&&((st==1&&(bestN==INT_MAX||bestN>N0*103/100))||(st==2&&bestN>=N0))){
      fprintf(stderr,"chain A gave bestN=%d against input %d: switching to chain B t=%.0f\n",bestN,N0,now()-t0g);
      chain=1;
      bestN=INT_MAX;
      st=firstinv();
    }
  }
  int useB=(bestN<=N0);
  if(useB&&!verify(B,bestN))useB=0;
  fprintf(stderr,"final N=%d\n",useB?bestN:N0);
  unsigned char (*R)[K]=useB?B:I0;
  int n=useB?bestN:N0;
  for(int r=0;r<n;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",R[r][c]);printf("\n");}
  return 0;
}
