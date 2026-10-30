#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <signal.h>
#define K 8
#define V 6
#define NT 1296
#define NS 70
#define MAXN 2100
static volatile sig_atomic_t stopf=0;
static void onsig(int s){(void)s;stopf=1;}
static int N,bestN=0,snapN=0,inN=0;
static unsigned char G[MAXN][K],BEST[MAXN][K],SNAP[MAXN][K],INP[MAXN][K];
static int cs[NS][4],cj[K][35],cp[K][35];static const int pw[4]={1,6,36,216};
static uint16_t idx[MAXN][NS],cnt[NS*NT],xr[NS*NT];static int uq[MAXN];
static int ulist[NS*NT],upos[NS*NT],unc;
static long tabu[MAXN][K];static long step=0;
static double t0g,TL=3440;
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addhole(int id){upos[id]=unc;ulist[unc++]=id;}
static inline void remhole(int id){int p=upos[id];int last=ulist[--unc];ulist[p]=last;upos[last]=p;}
static inline int cidx(const unsigned char*g,int j){return g[cs[j][0]]+6*g[cs[j][1]]+36*g[cs[j][2]]+216*g[cs[j][3]];}
static void build(){
  memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);
  for(int r=0;r<N;r++)for(int j=0;j<NS;j++){int t=cidx(G[r],j);idx[r][j]=t;cnt[j*NT+t]++;xr[j*NT+t]^=r;}
  unc=0;memset(uq,0,sizeof uq);
  for(int id=0;id<NS*NT;id++){if(cnt[id]==0)addhole(id);else if(cnt[id]==1)uq[xr[id]]++;}
}
static void setcell(int r,int c,int b){
  int a=G[r][c];if(a==b)return;
  for(int t=0;t<35;t++){int j=cj[c][t],p=cp[c][t];
    int o=idx[r][j],n=o+(b-a)*pw[p];int ido=j*NT+o,idn=j*NT+n;
    int k=cnt[ido]--;xr[ido]^=r;
    if(k==1){uq[r]--;addhole(ido);}else if(k==2){uq[xr[ido]]++;}
    k=cnt[idn]++;
    if(k==0){uq[r]++;remhole(idn);}else if(k==1){uq[xr[idn]]--;}
    xr[idn]^=r;idx[r][j]=n;}
  G[r][c]=b;
}
static void delrow(int r){memcpy(G[r],G[N-1],K);memcpy(tabu[r],tabu[N-1],sizeof(long)*K);N--;build();}
static void saveBest(){
  if(unc==0&&(bestN==0||N<bestN)){bestN=N;memcpy(BEST,G,(size_t)K*N);memcpy(SNAP,G,(size_t)K*N);snapN=N;
    fprintf(stderr,"FOUND N=%d at %.1fs\n",N,now()-t0g);fprintf(stderr,"PROGRESS %d\n",-N);}
}
static void pickdelete(){
  int mn=1<<30;for(int r=0;r<N;r++)if(uq[r]<mn)mn=uq[r];
  int lim=mn+(rnd()%3==0?2:0);
  static int c[MAXN];int nc=0;for(int r=0;r<N;r++)if(uq[r]<=lim)c[nc++]=r;
  int r=c[rnd()%nc];fprintf(stderr,"delete uq=%d (min %d) t=%.0f\n",uq[r],mn,now()-t0g);delrow(r);
}
static double phi[160];
static int ds[80],dv[80],nd;
static inline void addd(int s,int v){for(int i=0;i<nd;i++)if(ds[i]==s){dv[i]+=v;return;}ds[nd]=s;dv[nd]=v;nd++;}
static void freeDeletes(){
  int again=1;while(again&&unc==0){again=0;for(int r=0;r<N;r++)if(uq[r]==0){delrow(r);saveBest();again=1;break;}}
}
static void runA(double dur,double T0,double T1){
  double ts=now(),lastlog=ts;long it=0;double T=T0;
  freeDeletes();
  while(!stopf){
    if((++it&4095)==0){double tn=now();if(tn-t0g>TL||tn-ts>dur)break;
      T=T0*pow(T1/T0,(tn-ts)/dur);
      if(tn-lastlog>20){lastlog=tn;int mn=1<<30;double S=0;for(int r=0;r<N;r++){if(uq[r]<mn)mn=uq[r];S+=phi[uq[r]];}
        fprintf(stderr,"A t=%.0f N=%d T=%.2f minuq=%d obj=%.1f\n",tn-t0g,N,T,mn,S);}}
    int r=rnd()%N,c=rnd()%K,b=rnd()%(V-1);if(b>=G[r][c])b++;
    int a=G[r][c];int ok=1;
    for(int t=0;t<35;t++){int j=cj[c][t];if(cnt[j*NT+idx[r][j]]==1){ok=0;break;}}
    if(!ok)continue;
    nd=0;
    for(int t=0;t<35;t++){int j=cj[c][t],p=cp[c][t];int o=idx[r][j],n=o+(b-a)*pw[p];
      int ido=j*NT+o,idn=j*NT+n;
      if(cnt[ido]==2)addd(xr[ido]^r,1);
      if(cnt[idn]==1)addd(xr[idn],-1);}
    double dc=0;for(int i=0;i<nd;i++){int u=uq[ds[i]];dc+=phi[u+dv[i]]-phi[u];}
    if(dc>0&&urand()>=exp(-dc/T))continue;
    setcell(r,c,b);
    if(unc==0){int z=-1;for(int i=0;i<nd;i++)if(uq[ds[i]]==0){z=ds[i];break;}
      if(z>=0){delrow(z);saveBest();freeDeletes();}}
  }
}
static int cr[MAXN],cbest[MAXN];
static void runB(double dur){
  double last=now();
  memset(tabu,0,sizeof tabu);
  pickdelete();int bu=unc;double T=0.6;
  while(!stopf){
    step++;
    if((step&255)==0){double tn=now();if(tn-t0g>TL||tn-last>dur)return;}
    if(unc==0){saveBest();pickdelete();bu=unc;last=now();continue;}
    int id=ulist[rnd()%unc];int j=id/NT,t=id%NT;
    int v[4];for(int i=0;i<4;i++)v[i]=(t/pw[i])%V;
    int nc=0;
    for(int r=0;r<N;r++){int d=0;for(int i=0;i<4;i++)d+=(G[r][cs[j][i]]!=v[i]);if(d<=2)cr[nc++]=r;}
    int bd=1<<30,nb=0;
    for(int q=0;q<nc;q++){int r=cr[q];int istabu=0;unsigned char sv[4];
      for(int i=0;i<4;i++){int c=cs[j][i];sv[i]=G[r][c];if(sv[i]!=v[i]&&tabu[r][c]>step)istabu=1;}
      int u0=unc;
      for(int i=0;i<4;i++)setcell(r,cs[j][i],v[i]);
      int d=unc-u0;
      for(int i=3;i>=0;i--)setcell(r,cs[j][i],sv[i]);
      if(istabu&&unc+d>=bu)continue;
      if(d<bd){bd=d;nb=0;}
      if(d==bd)cbest[nb++]=r;
    }
    int r;
    if(nb==0||urand()<0.02){if(nc==0)continue;r=cr[rnd()%nc];}
    else{if(bd>0&&urand()>exp(-bd/T))continue;r=cbest[rnd()%nb];}
    for(int i=0;i<4;i++){int c=cs[j][i];if(G[r][c]!=v[i]){setcell(r,c,v[i]);tabu[r][c]=step+8+rnd()%8;}}
    if(unc<bu){bu=unc;}
  }
}
static int verify(unsigned char a[][K],int n){
  static unsigned char cc[NS*NT];memset(cc,0,sizeof cc);
  for(int r=0;r<n;r++)for(int j=0;j<NS;j++)cc[j*NT+cidx(a[r],j)]=1;
  int m=0;for(int i=0;i<NS*NT;i++)if(!cc[i])m++;return m;
}
int main(){
  signal(SIGTERM,onsig);signal(SIGINT,onsig);
  t0g=now();rs^=(uint64_t)time(NULL)*2654435761ULL;for(int i=0;i<10;i++)rnd();
  int ns=0;for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){cs[ns][0]=a;cs[ns][1]=b;cs[ns][2]=c;cs[ns][3]=d;ns++;}
  int nn[K]={0};
  for(int j=0;j<NS;j++)for(int i=0;i<4;i++){int c=cs[j][i];cj[c][nn[c]]=j;cp[c][nn[c]]=i;nn[c]++;}
  for(int u=0;u<160;u++)phi[u]=10.0*sqrt((double)u);
  const char*f=getenv("BEST_CANDIDATE");FILE*fp=f?fopen(f,"r"):NULL;
  if(!fp){fprintf(stderr,"no candidate\n");return 0;}
  int x,k=0;while(fscanf(fp,"%d",&x)==1&&inN<MAXN){INP[inN][k++]=(unsigned char)x;if(k==K){k=0;inN++;}}fclose(fp);
  N=inN;memcpy(G,INP,(size_t)K*N);build();
  fprintf(stderr,"input N=%d unc=%d\n",N,unc);
  if(unc!=0){for(int r=0;r<inN;r++)for(int c=0;c<K;c++)printf("%d%c",INP[r][c],c==K-1?'\n':' ');return 0;}
  bestN=N;memcpy(BEST,G,(size_t)K*N);memcpy(SNAP,G,(size_t)K*N);snapN=N;
  runA(150,3.0,0.3);
  while(!stopf&&now()-t0g<TL){
    memcpy(SNAP,G,(size_t)K*N);snapN=N;
    runB(240);
    if(stopf||now()-t0g>TL)break;
    N=snapN;memcpy(G,SNAP,(size_t)K*N);build();
    if(unc!=0){N=bestN;memcpy(G,BEST,(size_t)K*N);build();}
    runA(90,2.0,0.3);
  }
  fprintf(stderr,"end bestN=%d\n",bestN);
  if(verify(BEST,bestN)==0&&bestN<=inN){for(int r=0;r<bestN;r++)for(int c=0;c<K;c++)printf("%d%c",BEST[r][c],c==K-1?'\n':' ');}
  else{for(int r=0;r<inN;r++)for(int c=0;c<K;c++)printf("%d%c",INP[r][c],c==K-1?'\n':' ');}
  return 0;
}
