#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#define K 8
#define V 6
#define NT 1296
#define MAXN 2100
static int N;
static unsigned char G[MAXN][K], BEST[MAXN][K];
static int bestN;
static int cs[70][4], ncs=0;
static int cj[K][35], cp[K][35];
static int pw[4]={1,6,36,216};
static uint16_t cnt[70][NT];
static int unc;
static int ulist[70*NT], upos[70*NT];
static int tabu[MAXN][K];
static long step;
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static int listmode=1;
static inline void addU(int id){if(!listmode)return;upos[id]=unc;ulist[unc]=id;}
static inline void remU(int id){if(!listmode)return;int p=upos[id];int last=ulist[unc];ulist[p]=last;upos[last]=p;}
/* note: caller updates unc after add/rem calls ordering below */
static inline int tidx(int r,int j){int s=0;for(int i=0;i<4;i++)s+=G[r][cs[j][i]]*pw[i];return s;}
static void setcell(int r,int c,int b){
  int a=G[r][c]; if(a==b)return;
  for(int t=0;t<35;t++){int j=cj[c][t],p=cp[c][t];
    int o=tidx(r,j); int nw=o+(b-a)*pw[p];
    int id=j*NT+o;
    if(--cnt[j][o]==0){ if(listmode){upos[id]=unc;ulist[unc]=id;} unc++; }
    id=j*NT+nw;
    if(cnt[j][nw]++==0){ unc--; if(listmode){int pp=upos[id];int last=ulist[unc];ulist[pp]=last;upos[last]=pp;} }
  }
  G[r][c]=b;
}
static void build(){
  memset(cnt,0,sizeof cnt);
  for(int r=0;r<N;r++)for(int j=0;j<70;j++)cnt[j][tidx(r,j)]++;
  unc=0;for(int j=0;j<70;j++)for(int t=0;t<NT;t++)if(!cnt[j][t]){int id=j*NT+t;upos[id]=unc;ulist[unc++]=id;}
}
static int uniqcount(int r){int u=0;for(int j=0;j<70;j++)if(cnt[j][tidx(r,j)]==1)u++;return u;}
static void delrow(int r){
  for(int j=0;j<70;j++){int o=tidx(r,j);if(--cnt[j][o]==0){int id=j*NT+o;upos[id]=unc;ulist[unc++]=id;}}
  memcpy(G[r],G[N-1],K);memcpy(tabu[r],tabu[N-1],sizeof(int)*K);N--;
}
static void pickdelete(){
  int u[MAXN],mn=1<<30;
  for(int r=0;r<N;r++){u[r]=uniqcount(r);if(u[r]<mn)mn=u[r];}
  int lim=mn+(rnd()%3==0?2:0);
  int c[MAXN],nc=0;for(int r=0;r<N;r++)if(u[r]<=lim)c[nc++]=r;
  delrow(c[rnd()%nc]);
}
int main(){
  double t0=now(),TL=3480;
  rs^=(uint64_t)time(NULL)*2654435761ULL;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){cs[ncs][0]=a;cs[ncs][1]=b;cs[ncs][2]=c;cs[ncs][3]=d;ncs++;}
  int nn[K]={0};
  for(int j=0;j<70;j++)for(int i=0;i<4;i++){int c=cs[j][i];cj[c][nn[c]]=j;cp[c][nn[c]]=i;nn[c]++;}
  const char*f=getenv("BEST_CANDIDATE");FILE*fp=f?fopen(f,"r"):NULL;
  if(!fp){fprintf(stderr,"no candidate\n");return 0;}
  N=0;int x,k=0;while(fscanf(fp,"%d",&x)==1&&N<MAXN){G[N][k++]=x;if(k==K){k=0;N++;}}fclose(fp);
  build();
  if(unc==0){bestN=N;memcpy(BEST,G,sizeof(unsigned char)*K*N);} else {bestN=0;fprintf(stderr,"input invalid unc=%d\n",unc);}
  if(bestN==0){for(int r=0;r<N;r++){for(int c=0;c<K;c++)printf("%d%c",G[r][c],c==K-1?'\n':' ');}return 0;}
  memset(tabu,0,sizeof tabu);
  pickdelete();
  int bu=unc; long lastimp=0; double T=0.6; long stall=3000000;
  fprintf(stderr,"start N=%d unc=%d\n",N,unc);
  int cr[MAXN];int cbest[MAXN];
  while(1){
    step++;
    if((step&1023)==0){if(now()-t0>TL)break;}
    if(unc==0){
      bestN=N;memcpy(BEST,G,(size_t)K*N);
      fprintf(stderr,"FOUND N=%d at %.1fs\n",N,now()-t0);
      pickdelete();bu=unc;lastimp=step;continue;
    }
    if(step-lastimp>stall){
      N=bestN;memcpy(G,BEST,(size_t)K*N);memset(tabu,0,sizeof tabu);build();pickdelete();bu=unc;lastimp=step;
      fprintf(stderr,"restart unc=%d t=%.0f\n",unc,now()-t0);continue;
    }
    int id=ulist[rnd()%unc];int j=id/NT,t=id%NT;
    int v[4];for(int i=0;i<4;i++)v[i]=(t/pw[i])%V;
    int nc=0;
    for(int r=0;r<N;r++){int d=0;for(int i=0;i<4;i++)d+=(G[r][cs[j][i]]!=v[i]);if(d<=2)cr[nc++]=r;}
    int bd=1<<30,nb=0;
    listmode=0;
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
    listmode=1;
    int r;
    if(nb==0||urand()<0.02){ if(nc==0)continue; r=cr[rnd()%nc]; }
    else { if(bd>0&&urand()>exp(-bd/T))continue; r=cbest[rnd()%nb]; }
    for(int i=0;i<4;i++){int c=cs[j][i];if(G[r][c]!=v[i]){setcell(r,c,v[i]);tabu[r][c]=step+8+rnd()%8;}}
    if(unc<bu){bu=unc;lastimp=step;}
  }
  fprintf(stderr,"end bestN=%d steps=%ld\n",bestN,step);
  for(int r=0;r<bestN;r++)for(int c=0;c<K;c++)printf("%d%c",BEST[r][c],c==K-1?'\n':' ');
  return 0;
}
