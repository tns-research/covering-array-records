#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#define T 3
#define K 13
#define V 10
#define VT (V*V*V)
#define PC ((K-1)*(K-2)/2)
#define NTRI (K*(K-1)*(K-2)/6)
#define NTUP (NTRI*VT)
#define MAXN 6000
#define TLIM 3500.0
#define ALIM 450.0
#define DRT 150.0
#define DTH 6
#define DSTALL 300000

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
static void rowrem(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(--cnt[i]==0)addh(i);}}
static void rowadd(int r){for(int t=0;t<NTRI;t++){int i=tidx(t,A[r]);if(cnt[i]++==0)remh(i);}}
static long long rowgain(unsigned char*x){long long g=0;for(int t=0;t<NTRI;t++){int i=tidx(t,x);if(!cnt[i])g+=w[i];}return g;}
static inline long long colgain(unsigned char*x,int c,int v){
  long long g=0;
  for(int j=0;j<PC;j++){int idx=ct_tid[c][j]*VT+ct_m1[c][j]*x[ct_o1[c][j]]+ct_m2[c][j]*x[ct_o2[c][j]]+ct_mc[c][j]*v;if(!cnt[idx])g+=w[idx];}
  return g;
}
/* evaluate a row re-solve of row r forced to cover the hole; returns weighted delta, new row in out */
static long long evalres(int r,int c1,int v1,int c2,int v2,int c3,int v3,unsigned char*out){
  rowrem(r);
  long long g0=rowgain(A[r]);
  unsigned char x[K];memcpy(x,A[r],K);x[c1]=v1;x[c2]=v2;x[c3]=v3;
  int ord[K];for(int c=0;c<K;c++)ord[c]=c;
  for(int sw=0;sw<3;sw++){
    for(int i=0;i<K;i++){int j=i+rnd()%(K-i);int y=ord[i];ord[i]=ord[j];ord[j]=y;}
    int ch=0;
    for(int i=0;i<K;i++){int c=ord[i];if(c==c1||c==c2||c==c3)continue;
      long long bg=-1;int bv=x[c],tie=0;
      for(int v=0;v<V;v++){long long g=colgain(x,c,v);
        if(g>bg){bg=g;bv=v;tie=1;}else if(g==bg){tie++;if(rnd()%tie==0)bv=v;}}
      if(bv!=x[c]){x[c]=(unsigned char)bv;ch=1;}
    }
    if(!ch)break;
  }
  long long g1=rowgain(x);
  rowadd(r);
  memcpy(out,x,K);
  return g1-g0;
}
/* loss-free move test: returns 0 if (r,c)->nv breaks a uniquely covered tuple; else net drop of R-unique count */
static inline int nmove(int R,int r,int c,int nv,int*net){
  int ov=A[r][c],n=0;
  for(int j=0;j<PC;j++){
    int t=ct_tid[c][j];
    int base=t*VT+ct_m1[c][j]*A[r][ct_o1[c][j]]+ct_m2[c][j]*A[r][ct_o2[c][j]];
    int io=base+ct_mc[c][j]*ov,in=base+ct_mc[c][j]*nv;
    if(cnt[io]==1)return 0;
    int rt=tidx(t,A[R]);
    if(cnt[io]==2&&rt==io)n--;
    if(cnt[in]==1&&rt==in)n++;
  }
  *net=n;return 1;
}
static void drain(int R,double tl){
  static int U[NTRI];double ts=now();
  int best=1<<30,k0=-1,k=0;long long s=0,li=0;
  while(!stopf){
    s++;
    if((s&255)==0&&now()-ts>tl)break;
    k=0;for(int t=0;t<NTRI;t++){int i=tidx(t,A[R]);if(cnt[i]==1)U[k++]=i;}
    if(k0<0)k0=k;
    if(k<best){best=k;li=s;}
    if(k<=DTH||s-li>DSTALL)break;
    int bn=0,br=-1,bc=0,bv=0,ties=0;
    if(rnd()%8){
      int h=U[rnd()%k];int t=h/VT,v=h%VT;
      int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
      for(int r=0;r<N;r++){
        if(r==R)continue;
        unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
        if(m1+m2+m3!=2)continue;
        int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
        int net;if(!nmove(R,r,c,nv,&net))continue;
        if(net>bn){bn=net;br=r;bc=c;bv=nv;ties=1;}
        else if(net==bn&&bn>0){ties++;if(rnd()%ties==0){br=r;bc=c;bv=nv;}}
      }
    }
    if(br>=0){apply(br,bc,bv);continue;}
    for(int a=0;a<40;a++){
      int r=rnd()%N;if(r==R)continue;int c=rnd()%K,v=rnd()%V;if(v==A[r][c])continue;
      int net;if(!nmove(R,r,c,v,&net))continue;
      if(net>=0||rnd()%20==0){apply(r,c,v);break;}
    }
  }
  fprintf(stderr,"drain uniq %d -> %d (best %d) steps=%lld holes=%d\n",k0,k,best,s,nh);
}
static void outbest(void){
  for(int r=0;r<BN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",B[r][c]);printf("\n");}
  fflush(stdout);
}
static int pickrow(int spread){
  static int u[MAXN];static int ord[MAXN];
  for(int r=0;r<N;r++){int s=0;for(int t=0;t<NTRI;t++)if(cnt[tidx(t,A[r])]==1)s++;u[r]=s;ord[r]=r;}
  for(int i=0;i<N;i++){int j=i+rnd()%(N-i);int x=ord[i];ord[i]=ord[j];ord[j]=x;}
  // partial selection of the spread smallest
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
  if(fn){FILE*f=fopen(fn,"r");if(f){int x,c=0;while(N<MAXN&&fscanf(f,"%d",&x)==1){A[N][c++]=(unsigned char)(x%V);if(c==K){c=0;N++;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 0;}
  rebuild();
  BN=N;memcpy(B,A,(size_t)N*K);
  int bestvalid=(nh==0);
  fprintf(stderr,"loaded N=%d holes=%d\n",N,nh);
  if(!bestvalid)BN=0;
  for(int i=0;i<NTUP;i++)w[i]=1;
  memset(tabu,0,sizeof(tabu));
  long long it=0;double ast=now();int spread=1;long long bumps=0;int minh=1<<30;
  while(!stopf){
    if(nh==0){
      if(!bestvalid||N<BN){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;fprintf(stderr,"PROGRESS %d\n",-N);fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);}
      int r=pickrow(spread);
      if(now()-t0<TLIM-ALIM){drain(r,DRT);if(nh==0&&N==BN)memcpy(B,A,(size_t)N*K);}
      if(stopf)break;
      delrow(r);
      for(int i=0;i<NTUP;i++)w[i]=1;
      memset(tabu,0,sizeof(tabu));ast=now();minh=nh;spread=1;
      fprintf(stderr,"N=%d holes=%d\n",N,nh);
    }
    it++;
    if((it&1023)==0){
      double tn=now();
      if(tn-t0>TLIM)break;
      if(tn-ast>ALIM&&bestvalid){
        fprintf(stderr,"stall N=%d minh=%d, restore\n",N,minh);
        N=BN;memcpy(A,B,(size_t)N*K);rebuild();spread=20;
        int r=pickrow(spread);
        if(tn-t0<TLIM-ALIM){drain(r,DRT);if(nh==0&&N==BN)memcpy(B,A,(size_t)N*K);}
        if(stopf)break;
        delrow(r);
        for(int i=0;i<NTUP;i++)w[i]=1;memset(tabu,0,sizeof(tabu));ast=now();minh=nh;spread=1;
        continue;
      }
    }
    int h=holes[rnd()%nh];int t=h/VT,v=h%VT;
    int c1=tri[t][0],c2=tri[t][1],c3=tri[t][2],v1=v/(V*V),v2=(v/V)%V,v3=v%V;
    long long bd=-(1LL<<60);int br=-1,bc=-1,bv=0,ties=0;
    int noise=(rnd()%100)<2;
    int cand[6],nc=0,seenc=0;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int m1=x[c1]==v1,m2=x[c2]==v2,m3=x[c3]==v3;
      int mm=m1+m2+m3;
      if(mm==1||mm==2){seenc++;if(nc<6)cand[nc++]=r;else{int q=rnd()%seenc;if(q<6)cand[q]=r;}}
      if(mm!=2)continue;
      int c,nv;if(!m1){c=c1;nv=v1;}else if(!m2){c=c2;nv=v2;}else{c=c3;nv=v3;}
      if(tabu[r][c]>it)continue;
      long long d=noise?(long long)(rnd()%1000):delta(r,c,nv);
      if(d>bd){bd=d;br=r;bc=c;bv=nv;ties=1;}
      else if(d==bd){ties++;if(rnd()%ties==0){br=r;bc=c;bv=nv;}}
    }
    if(!noise&&(br<0||bd<=0)&&nc>0&&(br<0||rnd()%3==0)){
      long long rb=-(1LL<<60);int rr=-1;unsigned char rx[K],tx[K];
      for(int i=0;i<nc;i++){
        int r=cand[i];
        if(tabu[r][c1]>it&&A[r][c1]!=v1)continue;
        if(tabu[r][c2]>it&&A[r][c2]!=v2)continue;
        if(tabu[r][c3]>it&&A[r][c3]!=v3)continue;
        long long d=evalres(r,c1,v1,c2,v2,c3,v3,tx);
        if(d>rb||(d==rb&&rnd()%2)){rb=d;rr=r;memcpy(rx,tx,K);}
      }
      if(rr>=0&&(br<0||rb>=bd)){
        if(rb<=0){
          for(int i=0;i<nh;i++)w[holes[i]]++;
          if(++bumps%3000==0){for(int i=0;i<NTUP;i++)if(w[i]>1)w[i]=1+(w[i]-1)*3/4;}
        }
        for(int c=0;c<K;c++)if(rx[c]!=A[rr][c])tabu[rr][c]=it+3+rnd()%8;
        rowrem(rr);memcpy(A[rr],rx,K);rowadd(rr);
        if(nh<minh)minh=nh;
        continue;
      }
    }
    if(br<0){
      int r=rnd()%N;apply(r,c1,v1);apply(r,c2,v2);apply(r,c3,v3);
      tabu[r][c1]=tabu[r][c2]=tabu[r][c3]=it+5+rnd()%6;
    }else{
      if(bd<=0&&!noise){
        for(int i=0;i<nh;i++)w[holes[i]]++;
        if(++bumps%3000==0){for(int i=0;i<NTUP;i++)if(w[i]>1)w[i]=1+(w[i]-1)*3/4;}
      }
      apply(br,bc,bv);tabu[br][bc]=it+3+rnd()%8;
    }
    if(nh<minh)minh=nh;
  }
  if(nh==0&&(!bestvalid||N<BN)){BN=N;memcpy(B,A,(size_t)N*K);bestvalid=1;}
  if(!bestvalid){BN=N;memcpy(B,A,(size_t)N*K);}
  fprintf(stderr,"final N=%d\n",BN);
  outbest();
  return 0;
}
