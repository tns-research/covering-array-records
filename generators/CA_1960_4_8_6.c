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
static void evalmove(int r,int c,int b,int*du,double*dp){
  int a=G[r][c];nd=0;int d=0;
  for(int t=0;t<35;t++){int j=cj[c][t],p=cp[c][t];int o=idx[r][j],n=o+(b-a)*pw[p];
    int ido=j*NT+o,idn=j*NT+n;
    if(cnt[ido]==1){d++;addd(r,-1);}else if(cnt[ido]==2)addd(xr[ido]^r,1);
    if(cnt[idn]==0){d--;addd(r,1);}else if(cnt[idn]==1)addd(xr[idn],-1);}
  double s=0;for(int i=0;i<nd;i++){int u=uq[ds[i]];s+=phi[u+dv[i]]-phi[u];}
  *du=d;*dp=s;
}
static void runC(double dur,double T0,double T1,double L,double per){
  double ts=now(),last=ts,lastlog=ts;long it=0;double T=T0;
  freeDeletes();if(unc==0)pickdelete();
  int bu=unc;
  while(!stopf){
    if((++it&4095)==0){double tn=now();if(tn-t0g>TL||tn-last>dur)break;
      double fr=fmod(tn-ts,per)/per;T=T0*pow(T1/T0,fr);
      if(tn-lastlog>20){lastlog=tn;fprintf(stderr,"C t=%.0f N=%d T=%.2f L=%.1f unc=%d bu=%d\n",tn-t0g,N,T,L,unc,bu);}}
    if(unc==0){saveBest();freeDeletes();pickdelete();bu=unc;last=now();continue;}
    int r,c=-1,b=0;
    if((rnd()&15)==0){
      int id=ulist[rnd()%unc];int j=id/NT,t=id%NT;int v[4];for(int i=0;i<4;i++)v[i]=(t/pw[i])%V;
      int n1=0,n2=0;
      for(int q=0;q<N;q++){int d=0;for(int i=0;i<4;i++)d+=(G[q][cs[j][i]]!=v[i]);
        if(d==1)cr[n1++]=q;else if(d==2)cbest[n2++]=q;}
      if(n1>0)r=cr[rnd()%n1];else if(n2>0)r=cbest[rnd()%n2];else continue;
      int m=0;for(int i=0;i<4;i++)if(G[r][cs[j][i]]!=v[i]){m++;if(rnd()%m==0){c=cs[j][i];b=v[i];}}
      if(c<0)continue;
    } else {r=rnd()%N;c=rnd()%K;b=rnd()%(V-1);if(b>=G[r][c])b++;}
    int du;double dp;evalmove(r,c,b,&du,&dp);
    double dc=L*du+dp;
    if(dc>0&&urand()>=exp(-dc/T))continue;
    setcell(r,c,b);
    if(unc<bu)bu=unc;
  }
}
static void runB(double dur,int del){
  double last=now();
  memset(tabu,0,sizeof tabu);
  if(del||unc==0)pickdelete();
  int bu=unc;double T=0.6;
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
#define MAXH 1500
#define SIDE 60
static unsigned char colmask[NS],um[MAXN][70];static int umn[MAXN];
static int hN,hid[MAXH],hord[MAXH];
static unsigned char hc[MAXH][4],hv[MAXH][4];
static uint64_t cm[K][V],emptyq;
static signed char prow[64][K];
static long nodes,nodelim=3000,statH=0,statAbort=0;
static int Rl[64],inR[MAXN],rstamp=0,F[5000],rcs[NS*NT];static uint16_t rc[NS*NT];
static void buildUM(){
  for(int r=0;r<N;r++){int n=0;for(int j=0;j<NS;j++)if(cnt[j*NT+idx[r][j]]==1)um[r][n++]=colmask[j];umn[r]=n;}
}
static int dfs(int depth){
  if(depth==hN)return 1;
  if(++nodes>nodelim)return -1;
  int bp=-1,bn=1000;uint64_t bd=0;
  uint64_t lowe=emptyq&(~emptyq+1),keep=~emptyq|lowe;
  for(int p=depth;p<hN;p++){int i=hord[p];const unsigned char*c=hc[i],*v=hv[i];
    uint64_t d=cm[c[0]][v[0]]&cm[c[1]][v[1]]&cm[c[2]][v[2]]&cm[c[3]][v[3]]&keep;
    if(!d)return 0;
    int n=__builtin_popcountll(d);
    if(n<bn){bn=n;bd=d;bp=p;if(n==1)break;}}
  int i=hord[bp];hord[bp]=hord[depth];hord[depth]=i;
  int cl[64],sc[64],nc=0;
  while(bd){int q=__builtin_ctzll(bd);bd&=bd-1;int s=0;
    for(int t=0;t<4;t++)if(prow[q][hc[i][t]]<0)s++;
    s=s*16+(int)(rnd()&15);
    int p=nc++;while(p>0&&sc[p-1]>s){sc[p]=sc[p-1];cl[p]=cl[p-1];p--;}
    sc[p]=s;cl[p]=q;}
  for(int a=0;a<nc;a++){int q=cl[a];uint64_t bit=1ULL<<q;int b=0;
    for(int t=0;t<4;t++){int c=hc[i][t],v=hv[i][t];
      if(prow[q][c]<0){prow[q][c]=(signed char)v;b|=1<<t;for(int w=0;w<V;w++)if(w!=v)cm[c][w]&=~bit;}}
    uint64_t oe=emptyq;emptyq&=~bit;
    int res=dfs(depth+1);
    if(res==1)return 1;
    for(int t=0;t<4;t++)if((b>>t)&1){int c=hc[i][t],v=hv[i][t];
      prow[q][c]=-1;for(int w=0;w<V;w++)if(w!=v)cm[c][w]|=bit;}
    emptyq=oe;
    if(res<0)return -1;
  }
  return 0;
}
static void tryPlace(int id){
  int j=id/NT,t=id%NT;int c[4],v[4];uint64_t d=~0ULL;
  for(int i=0;i<4;i++){c[i]=cs[j][i];v[i]=(t/pw[i])%V;d&=cm[c[i]][v[i]];}
  if(!d)return;
  int q=__builtin_ctzll(d);uint64_t bit=1ULL<<q;
  for(int i=0;i<4;i++)if(prow[q][c[i]]<0){prow[q][c[i]]=(signed char)v[i];
    for(int w=0;w<V;w++)if(w!=v[i])cm[c[i]][w]&=~bit;}
}
static int trial(int seed,int m,int drop){
  rstamp++;int nR=0,nF=0,x=-1,nx=0;
  F[nF++]=seed;
  while(nR<m){
    int e=-1;
    if(nF>0){int p=rnd()%nF;e=F[p];F[p]=F[--nF];}
    else{for(int a=0;a<8;a++){int r0=Rl[rnd()%nR];int jj=rnd()%NS;e=jj*NT+idx[r0][jj];if(cnt[e]==1)break;}}
    int j=e/NT,t=e%NT;int c[4],v[4];
    for(int i=0;i<4;i++){c[i]=cs[j][i];v[i]=(t/pw[i])%V;}
    int br=-1,bs=1<<30;
    for(int r=0;r<N;r++){
      if(inR[r]==rstamp)continue;
      const unsigned char*g=G[r];int M=0,d=0;
      for(int i=0;i<4;i++)if(g[c[i]]!=v[i]){M|=1<<c[i];d++;}
      if(d>2||d==0)continue;
      int cost=0;for(int k=0;k<umn[r];k++)cost+=(um[r][k]&M)!=0;
      int s=cost*4+(int)(rnd()%8);
      if(s<bs){bs=s;br=r;}
    }
    if(br<0){do{br=rnd()%N;}while(inR[br]==rstamp);}
    inR[br]=rstamp;Rl[nR]=br;
    int M=0;for(int i=0;i<4;i++)if(G[br][c[i]]!=v[i])M|=1<<c[i];
    for(int jj=0;jj<NS;jj++){int id=jj*NT+idx[br][jj];
      if(cnt[id]==1&&(colmask[jj]&M)){F[nF++]=id;
        if(drop&&nR==0){nx++;if(rnd()%nx==0)x=id;}}}
    nR++;
  }
  hN=0;
  for(int a=0;a<nR;a++){int r=Rl[a];
    for(int j=0;j<NS;j++){int id=j*NT+idx[r][j];
      if(rcs[id]!=rstamp){rcs[id]=rstamp;rc[id]=0;}
      rc[id]++;
      if(rc[id]==cnt[id]){if(hN>=MAXH-1)return 0;hid[hN++]=id;}}}
  hid[hN++]=seed;
  if(x>=0){int f=-1;for(int i=0;i<hN;i++)if(hid[i]==x){f=i;break;}
    if(f>=0)hid[f]=hid[--hN];else x=-1;}
  statH+=hN;
  int perj[NS];memset(perj,0,sizeof perj);
  for(int i=0;i<hN;i++){if(++perj[hid[i]/NT]>nR)return 0;}
  for(int i=0;i<hN;i++){int id=hid[i],j=id/NT,t=id%NT;
    for(int q=0;q<4;q++){hc[i][q]=(unsigned char)cs[j][q];hv[i][q]=(unsigned char)((t/pw[q])%V);}
    hord[i]=i;}
  for(int i=hN-1;i>0;i--){int o=rnd()%(i+1);int tmp=hord[i];hord[i]=hord[o];hord[o]=tmp;}
  uint64_t all=(1ULL<<nR)-1;
  for(int c=0;c<K;c++)for(int v=0;v<V;v++)cm[c][v]=all;
  memset(prow,-1,sizeof prow);
  emptyq=all;nodes=0;
  int res=dfs(0);
  if(res<0)statAbort++;
  if(res!=1)return 0;
  for(int u=0;u<unc;u++)if(ulist[u]!=seed)tryPlace(ulist[u]);
  if(x>=0)tryPlace(x);
  int used[64];memset(used,0,sizeof used);
  for(int q=0;q<nR;q++){
    int nb=0;for(int c=0;c<K;c++)if(prow[q][c]>=0)nb++;
    if(!nb)continue;
    int bi=-1,ba=-1;
    for(int a=0;a<nR;a++){if(used[a])continue;int ag=0;
      for(int c=0;c<K;c++)if(prow[q][c]>=0&&G[Rl[a]][c]==prow[q][c])ag++;
      if(ag>ba){ba=ag;bi=a;}}
    used[bi]=1;
    for(int c=0;c<K;c++)if(prow[q][c]>=0)setcell(Rl[bi],c,prow[q][c]);
  }
  return 1;
}
static void runD(double dur){
  for(int j=0;j<NS;j++){colmask[j]=0;for(int i=0;i<4;i++)colmask[j]|=(unsigned char)(1<<cs[j][i]);}
  if(unc==0){saveBest();freeDeletes();pickdelete();}
  buildUM();
  int bu=unc;double last=now(),lastlog=last;long tr=0,succ=0,side=0,fails=0;statH=0;statAbort=0;
  while(!stopf){
    if((tr&15)==0){double tn=now();if(tn-t0g>TL||tn-last>dur)break;
      if(tn-lastlog>20){lastlog=tn;fprintf(stderr,"D t=%.0f N=%d unc=%d bu=%d tr=%ld succ=%ld side=%ld avgH=%.0f abort=%ld\n",
        tn-t0g,N,unc,bu,tr,succ,side,tr?(double)statH/tr:0.0,statAbort);}}
    tr++;
    if(unc==0){saveBest();freeDeletes();pickdelete();buildUM();bu=unc;last=now();continue;}
    int seed=ulist[rnd()%unc];
    int mx=(rnd()&1)?16:48;int m=3+(int)(rnd()%(mx-2));
    int drop=fails>=SIDE;
    int ok=trial(seed,m,drop);
    if(ok){buildUM();fails=0;if(drop)side++;else succ++;
      if(unc<bu){bu=unc;last=now();fprintf(stderr,"D improved unc=%d N=%d m=%d t=%.0f\n",unc,N,m,last-t0g);}}
    else fails++;
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
  int round=0;
  while(!stopf&&now()-t0g<TL){
    memcpy(SNAP,G,(size_t)K*N);snapN=N;
    if(round%2==1)runB(60,1);
    runD(240);
    round++;
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
