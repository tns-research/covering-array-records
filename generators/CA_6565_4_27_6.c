#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define K 27
#define NS 17550
#define NT 1296
#define NG (NS*NT)
#define CAP 7000
#define PC 2600
#define EVMAX 200
#define DUCAP 2
#define MARG 12
#define UPH 5
#define JCAP 100000
#define POLSEC 40.0
static int N;
static unsigned char A[CAP][K], B[CAP][K], I0[CAP][K];
static unsigned short cnt[NG], xr[NG];
static int w[NG];
static int mis[NG],pos[NG],M=0;
static int hv[NG],hpos[NG],H=0;
static int L[CAP][K],Lc[CAP][K];
static long long tb[CAP][K];
static int scol[NS][4];
static unsigned char mdc[NG][8];
typedef struct{int g,s;unsigned char o0,o1,o2,pad;unsigned short w0,w1,w2,wc;}E;
static E ce[K][PC];
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(int r,int s){unsigned char*x=A[r];return x[scol[s][0]]*216+x[scol[s][1]]*36+x[scol[s][2]]*6+x[scol[s][3]];}
static int U[CAP][NS],Uc[CAP],upos[NG];
static unsigned smask[NS];
static inline void addL(int r,int s,int dw,int dc,int g){int*sc=scol[s];
  if(dc>0){upos[g]=Uc[r];U[r][Uc[r]++]=g;}
  else if(dc<0){int p=upos[g];int l=U[r][--Uc[r]];U[r][p]=l;upos[l]=p;}
  L[r][sc[0]]+=dw;L[r][sc[1]]+=dw;L[r][sc[2]]+=dw;L[r][sc[3]]+=dw;
  Lc[r][sc[0]]+=dc;Lc[r][sc[1]]+=dc;Lc[r][sc[2]]+=dc;Lc[r][sc[3]]+=dc;}
static int oi[PC];
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  unsigned char*x=A[r]; int dd=b-a;
  E*ec=ce[c];
  for(int j=0;j<PC;j++){E*e=&ec[j];oi[j]=e->g+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2+a*e->wc;}
  for(int j=0;j<PC;j++){
    if(j+8<PC){int p=oi[j+8];int pn=p+dd*ec[j+8].wc;
      __builtin_prefetch(&cnt[p]);__builtin_prefetch(&xr[p]);__builtin_prefetch(&cnt[pn]);__builtin_prefetch(&xr[pn]);}
    int o=oi[j],n=o+dd*ec[j].wc,s=ec[j].s;
    unsigned short cO=--cnt[o]; xr[o]^=(unsigned short)r;
    if(cO==0){addm(o);addL(r,s,-w[o],-1,o);}
    else if(cO==1){addL(xr[o],s,w[o],1,o);}
    unsigned short cN=cnt[n]++;
    if(cN==0){delm(n);addL(r,s,w[n],1,n);}
    else if(cN==1){addL(xr[n],s,-w[n],-1,n);}
    xr[n]^=(unsigned short)r;
  }
  A[r][c]=b;
}
static void rebuild(){
  memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);
  for(int i=0;i<NG;i++)pos[i]=-1;
  M=0;
  for(int s=0;s<NS;s++){int*sc=scol[s];int g0=s*NT;
    for(int r=0;r<N;r++){unsigned char*x=A[r];int g=g0+x[sc[0]]*216+x[sc[1]]*36+x[sc[2]]*6+x[sc[3]];cnt[g]++;xr[g]^=(unsigned short)r;}}
  memset(L,0,sizeof L);memset(Lc,0,sizeof Lc);memset(Uc,0,sizeof Uc);
  for(int g=0;g<NG;g++){if(cnt[g]==0)addm(g);else if(cnt[g]==1)addL(xr[g],g/NT,w[g],1,g);}
}
static void resetw(){for(int i=0;i<H;i++){w[hv[i]]=1;hpos[hv[i]]=-1;}H=0;}
static void decayw(int num,int den){
  for(int i=0;i<H;){
    int g=hv[i];int ow=w[g];
    int nw=1+(int)((long long)(ow-1)*num/den); if(nw<1)nw=1;
    if(nw!=ow&&cnt[g]==1)addL(xr[g],g/NT,nw-ow,0,g);
    w[g]=nw;
    if(nw<=1){hpos[g]=-1;int l=hv[--H];hv[i]=l;hpos[l]=i;}
    else i++;
  }
}
static void delrow(int br){
  for(int s=0;s<NS;s++){int g=s*NT+tidx(br,s);unsigned short c=--cnt[g];xr[g]^=(unsigned short)br;
    if(c==0)addm(g);else if(c==1)addL(xr[g],s,w[g],1,g);}
  int last=N-1;
  if(br!=last){
    memcpy(A[br],A[last],K);memcpy(L[br],L[last],sizeof L[0]);memcpy(Lc[br],Lc[last],sizeof Lc[0]);
    memcpy(U[br],U[last],Uc[last]*sizeof(int));Uc[br]=Uc[last];
    unsigned short f=(unsigned short)(last^br);
    for(int s=0;s<NS;s++){int g=s*NT+tidx(br,s);xr[g]^=f;}
  }
  Uc[last]=0;
  N--;
}
static int pickrow(int rank){
  static int us[CAP],ord[CAP];
  for(int r=0;r<N;r++){int u=0;for(int c=0;c<K;c++)u+=Lc[r][c];us[r]=u/4;ord[r]=r;}
  for(int i=N-1;i>0;i--){int j=rnd()%(i+1);int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  if(rank>=N)rank=0;
  for(int i=0;i<=rank;i++){
    int bi=i;for(int j=i+1;j<N;j++)if(us[ord[j]]<us[ord[bi]])bi=j;
    int t=ord[i];ord[i]=ord[bi];ord[bi]=t;
  }
  int r=ord[rank];
  fprintf(stderr,"delete row %d unique=%d\n",r,us[r]);
  return r;
}
static inline int evald(int r,int c,int b,int*du){
  unsigned char*x=A[r];int gw=0,gc=0;
  for(int q=0;q<M;q++){unsigned char*u=mdc[q];
    if(u[0]!=c&&u[1]!=c&&u[2]!=c&&u[3]!=c)continue;
    int ok=1;
    for(int k=0;k<4;k++){int col=u[k],v=u[4+k];
      if(col==c){if(v!=b){ok=0;break;}}else if(x[col]!=v){ok=0;break;}}
    if(ok){gw+=w[mis[q]];gc++;}
  }
  *du=Lc[r][c]-gc;return L[r][c]-gw;
}
static inline int evalm(int r,int mask,int*sc,int*tv,int*du){
  unsigned char t[K];memcpy(t,A[r],K);
  unsigned cm=0;int n=0,c0=0;
  for(int k=0;k<4;k++)if(mask>>k&1){t[sc[k]]=(unsigned char)tv[k];cm|=1u<<sc[k];n++;c0=sc[k];}
  int loss=0,lc=0;
  if(n==1){loss=L[r][c0];lc=Lc[r][c0];}
  else{int*u=U[r];int m=Uc[r];for(int j=0;j<m;j++){int g=u[j];if(smask[g/NT]&cm){loss+=w[g];lc++;}}}
  int gw=0,gc=0;
  for(int q=0;q<M;q++){unsigned char*h=mdc[q];
    if(t[h[0]]==h[4]&&t[h[1]]==h[5]&&t[h[2]]==h[6]&&t[h[3]]==h[7]){gw+=w[mis[q]];gc++;}}
  *du=lc-gc;return loss-gw;
}
/* Called only with M==0: reassign free cells to drain uniques from weak rows, delete empty rows. */
static int pcell[CAP*K],pbase[PC];
static void polish(double t0,double TL,double budget){
  double ts=now();int stop=0,dels=0;long acc=0,evs=0;
  static int W[NS+2];
  W[0]=1<<20;for(int u=1;u<NS+2;u++)W[u]=(1<<20)/u;
  while(!stop){
    int F=0;
    for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(Lc[r][c]==0)pcell[F++]=r*K+c;
    for(int i=F-1;i>0;i--){int j=rnd()%(i+1);int t=pcell[i];pcell[i]=pcell[j];pcell[j]=t;}
    int improved=0;
    for(int i=0;i<F;i++){
      if((i&15)==0){double tn=now();if(tn-t0>TL||tn-ts>budget){stop=1;break;}}
      int r=pcell[i]/K,c=pcell[i]%K;
      if(Lc[r][c]!=0)continue;
      unsigned char*x=A[r];int a=x[c];E*ec=ce[c];
      long long gn[6]={0,0,0,0,0,0},lossv=0;
      for(int j=0;j<PC;j++){E*e=&ec[j];pbase[j]=e->g+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;}
      for(int j=0;j<PC;j++){
        if(j+8<PC){int p=pbase[j+8],wc8=ec[j+8].wc;for(int b=0;b<6;b++)__builtin_prefetch(&cnt[p+b*wc8]);}
        int p=pbase[j],wc=ec[j].wc;
        for(int b=0;b<6;b++){
          int idx=p+b*wc;int cv=cnt[idx];
          if(b==a){if(cv==2){int own=xr[idx]^r;lossv+=W[Uc[own]+1];}}
          else if(cv==1){int own=xr[idx];gn[b]+=W[Uc[own]];}
        }
      }
      evs++;
      int bb=-1;long long bg=0;
      for(int b=0;b<6;b++)if(b!=a&&gn[b]-lossv>bg){bg=gn[b]-lossv;bb=b;}
      if(bb>=0){change(r,c,bb);improved=1;acc++;}
    }
    for(int r=N-1;r>=0;r--)if(Uc[r]==0&&N>1){delrow(r);dels++;}
    if(!improved)break;
  }
  long long tu=0;int mn=1<<30,fr=0;
  for(int r=0;r<N;r++){tu+=Uc[r];if(Uc[r]<mn)mn=Uc[r];for(int c=0;c<K;c++)if(Lc[r][c]==0)fr++;}
  fprintf(stderr,"polish N=%d M=%d evals=%ld acc=%ld dels=%d uniq=%lld minU=%d free=%d t=%.0f\n",N,M,evs,acc,dels,tu,mn,fr,now()-t0);
}
static void output(unsigned char (*X)[K],int n){
  for(int r=0;r<n;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",X[r][c]);printf("\n");}
  fflush(stdout);
}
int main(){
  double t0=now(),TL=3380,STALL=240;
  const char*fn=getenv("BEST_CANDIDATE");
  N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;if(N>=CAP)break;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 1;}
  int N0=N; memcpy(I0,A,sizeof(A[0])*N);
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  for(s=0;s<NS;s++)smask[s]=(1u<<scol[s][0])|(1u<<scol[s][1])|(1u<<scol[s][2])|(1u<<scol[s][3]);
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={216,36,6,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){int c=scol[s][i];int j=cc[c]++;E*e=&ce[c][j];e->g=s*NT;e->s=s;e->wc=wts[i];
    int oc[3],ow[3],q=0;for(int k=0;k<4;k++)if(k!=i){oc[q]=scol[s][k];ow[q]=wts[k];q++;}
    e->o0=oc[0];e->o1=oc[1];e->o2=oc[2];e->w0=ow[0];e->w1=ow[1];e->w2=ow[2];}
  for(int i=0;i<NG;i++){w[i]=1;hpos[i]=-1;}
  rebuild();
  fprintf(stderr,"input N=%d missing=%d t=%.1f\n",N,M,now()-t0);
  int bestN=-1;
  if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
  static int cr[CAP],cmk[CAP],clb[CAP];
  static int Jr[JCAP],Jc[JCAP],Jv[JCAP]; int J=0;
  long long it=0,bumps=0; int bestMs=M; double tlast=now(),tlog=now(); int restarts=0;
  memset(tb,0,sizeof tb);
  while(1){
    if((it&63)==0){
      double tn=now();
      if(tn-t0>TL)break;
      if(bestN>0 && tn-tlast>STALL && M>0){
        restarts++;
        N=bestN; memcpy(A,B,sizeof(A[0])*N); resetw(); rebuild();
        memset(tb,0,sizeof tb);
        polish(t0,TL,POLSEC);
        bestN=N; memcpy(B,A,sizeof(A[0])*N);
        int br=pickrow((int)(rnd()%(3+restarts)));
        delrow(br); bestMs=M; tlast=now(); J=0;
        fprintf(stderr,"restart %d N=%d missing=%d t=%.0f\n",restarts,N,M,tn-t0);
        continue;
      }
      if(tn-tlog>30){tlog=tn;fprintf(stderr,"it=%lld N=%d M=%d bestMs=%d H=%d best=%d t=%.0f\n",it,N,M,bestMs,H,bestN,tn-t0);}
    }
    it++;
    if(M==0){
      if(bestN<0||N<bestN){bestN=N;memcpy(B,A,sizeof(A[0])*N);fprintf(stderr,"valid N=%d t=%.0f it=%lld\n",N,now()-t0,it);}
      polish(t0,TL,POLSEC);
      if(N<bestN)fprintf(stderr,"valid N=%d after polish t=%.0f\n",N,now()-t0);
      bestN=N; memcpy(B,A,sizeof(A[0])*N);
      int br=pickrow(0);
      delrow(br);
      memset(tb,0,sizeof tb);
      decayw(1,2);
      bestMs=M; tlast=now(); J=0;
      continue;
    }
    int g=mis[rnd()%M]; s=g/NT; int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int *sc=scol[s]; int nc=0; int minub=1<<30;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r]; int*Lr=L[r];
      int mask=0,lb=0,ub=0;
      for(int i=0;i<4;i++)if(x[sc[i]]!=tv[i]){int l=Lr[sc[i]];ub+=l;if(l>lb)lb=l;mask|=1<<i;}
      cr[nc]=r;cmk[nc]=mask;clb[nc]=lb;nc++;
      if(ub<minub)minub=ub;
    }
    for(int q=0;q<M;q++){int h=mis[q];int hs=h/NT,hid=h%NT;unsigned char*u=mdc[q];int*hc=scol[hs];
      u[0]=hc[0];u[1]=hc[1];u[2]=hc[2];u[3]=hc[3];u[4]=hid/216;u[5]=(hid/36)%6;u[6]=(hid/6)%6;u[7]=hid%6;}
    int pick=-1,bd=1<<30,ties=0,ev=0;
    int off=rnd()%nc;
    for(int ii=0;ii<nc&&ev<EVMAX;ii++){
      int i=ii+off; if(i>=nc)i-=nc;
      if(clb[i]>minub)continue;
      int r=cr[i],mask=cmk[i];
      int tabu=0;
      for(int k=0;k<4;k++)if((mask>>k&1)&&tb[r][sc[k]]>it)tabu=1;
      ev++;
      int du;int d=evalm(r,mask,sc,tv,&du);
      if(du>DUCAP)continue;
      if(tabu && !(M+du<bestMs))continue;
      if(d<bd){bd=d;pick=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)pick=i;}
    }
    int apply=0;
    if(pick>=0&&bd<=0)apply=1;
    else{
      for(int i=0;i<M;i++){int q=mis[i];if(w[q]==1){hpos[q]=H;hv[H++]=q;}w[q]++;}
      bumps++;
      if(bumps%2000==0)decayw(4,5);
      if(pick>=0&&rnd()%100<UPH)apply=1;
    }
    if(apply){
      int r=cr[pick],mask=cmk[pick];
      for(int k=0;k<4;k++)if(mask>>k&1){
        int c=sc[k];
        if(J<JCAP){Jr[J]=r;Jc[J]=c;Jv[J]=A[r][c];J++;}
        change(r,c,tv[k]);
        tb[r][c]=it+4+(long long)(rnd()%10);
      }
      if(M<=bestMs){if(M<bestMs)tlast=now();bestMs=M;J=0;}
      else if(M>bestMs+MARG||J+4>=JCAP){while(J>0){J--;change(Jr[J],Jc[J],Jv[J]);}}
    }
  }
  if(bestN<0){
    if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
    else {bestN=N0;memcpy(B,I0,sizeof(A[0])*N0);}
  }
  /* full verification */
  N=bestN; memcpy(A,B,sizeof(A[0])*N); resetw(); rebuild();
  fprintf(stderr,"final N=%d missing=%d iters=%lld\n",bestN,M,it);
  if(M==0) output(B,bestN);
  else { fprintf(stderr,"verification failed, printing input\n"); output(I0,N0); }
  return 0;
}
