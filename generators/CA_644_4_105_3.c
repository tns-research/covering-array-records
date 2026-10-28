#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#define K 105
#define MAXR 720
static unsigned char A[MAXR][K], SNAP[MAXR][K];
static int snapN, active[MAXR], act[MAXR], nact;
static uint16_t *cnt; static uint32_t *sm;
static int u[MAXR][K]; static long long tabu[MAXR][K];
static int B1[K+1],B2[K+1],B3[K+1],B4[K+1];
typedef struct{int rank,code;unsigned char cols[4],vals[4];}Hole;
static Hole *holes; static int nh,caph;
static long long iter; static uint64_t rs=88172645463325252ULL;
static inline uint64_t rng(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double T0;
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9-T0;}
static inline void hole_add(int rank,int code,const int*cols){
  if(nh==caph){caph=caph?caph*2:1024;holes=realloc(holes,sizeof(Hole)*caph);}
  Hole*h=&holes[nh++];h->rank=rank;h->code=code;
  for(int t=0;t<4;t++)h->cols[t]=cols[t];
  h->vals[0]=code/27;h->vals[1]=code/9%3;h->vals[2]=code/3%3;h->vals[3]=code%3;}
static inline void hole_del(int rank,int code){
  for(int i=0;i<nh;i++)if(holes[i].rank==rank&&holes[i].code==code){holes[i]=holes[--nh];return;}}
static unsigned char Ac[K][MAXR];
static void recompute(){
  nact=0;for(int r=0;r<snapN;r++)if(active[r])act[nact++]=r;
  for(int c=0;c<K;c++)for(int i=0;i<nact;i++)Ac[c][i]=A[act[i]][c];
  memset(u,0,sizeof(u));nh=0;
  static int pre[MAXR];int lc[81];uint32_t ls[81];
  for(int c3=3;c3<K;c3++)for(int c2=2;c2<c3;c2++)for(int c1=1;c1<c2;c1++){
    for(int i=0;i<nact;i++)pre[i]=9*Ac[c1][i]+3*Ac[c2][i]+Ac[c3][i];
    for(int c0=0;c0<c1;c0++){
      int rank=B1[c0]+B2[c1]+B3[c2]+B4[c3];
      memset(lc,0,sizeof lc);memset(ls,0,sizeof ls);
      unsigned char*a0=Ac[c0];
      for(int i=0;i<nact;i++){int code=pre[i]+27*a0[i];lc[code]++;ls[code]+=act[i];}
      long base=(long)rank*81;int cols[4]={c0,c1,c2,c3};
      for(int t=0;t<81;t++){cnt[base+t]=lc[t];sm[base+t]=ls[t];
        if(lc[t]==1){int r=ls[t];u[r][c0]++;u[r][c1]++;u[r][c2]++;u[r][c3]++;}
        else if(lc[t]==0)hole_add(rank,t,cols);}
    }}
}
static void change(int r,int c,int nv){
  int ov=A[r][c];if(ov==nv)return;
  int o[K],n=0;for(int x=0;x<K;x++)if(x!=c)o[n++]=x;
  unsigned char*row=A[r];static const int W[4]={27,9,3,1};
  for(int i=0;i<n;i++){int a=o[i];for(int j=i+1;j<n;j++){int b=o[j];for(int k=j+1;k<n;k++){int d=o[k];
    int cols[4],p;
    if(c<a){cols[0]=c;cols[1]=a;cols[2]=b;cols[3]=d;p=0;}
    else if(c<b){cols[0]=a;cols[1]=c;cols[2]=b;cols[3]=d;p=1;}
    else if(c<d){cols[0]=a;cols[1]=b;cols[2]=c;cols[3]=d;p=2;}
    else{cols[0]=a;cols[1]=b;cols[2]=d;cols[3]=c;p=3;}
    int rank=B1[cols[0]]+B2[cols[1]]+B3[cols[2]]+B4[cols[3]];
    int base=0;for(int t=0;t<4;t++)if(t!=p)base+=row[cols[t]]*W[t];
    int co=base+ov*W[p],cn=base+nv*W[p];
    long io=(long)rank*81+co,in=(long)rank*81+cn;
    cnt[io]--;sm[io]-=r;
    if(cnt[io]==0){for(int t=0;t<4;t++)u[r][cols[t]]--;hole_add(rank,co,cols);}
    else if(cnt[io]==1){int r2=sm[io];for(int t=0;t<4;t++)u[r2][cols[t]]++;}
    int pc=cnt[in];
    if(pc==0){hole_del(rank,cn);for(int t=0;t<4;t++)u[r][cols[t]]++;}
    else if(pc==1){int r2=sm[in];for(int t=0;t<4;t++)u[r2][cols[t]]--;}
    cnt[in]++;sm[in]+=r;
  }}}
  row[c]=nv;
}
static int gain(int r,int c,int v){
  int g=0;for(int i=0;i<nh;i++){Hole*h=&holes[i];int t=-1;
    for(int s=0;s<4;s++)if(h->cols[s]==c){t=s;break;}
    if(t<0||h->vals[t]!=v)continue;int ok=1;
    for(int s=0;s<4;s++)if(s!=t&&A[r][h->cols[s]]!=h->vals[s]){ok=0;break;}
    g+=ok;}
  return g;}
static void apply(int r,int c,int v){change(r,c,v);tabu[r][c]=iter+8+rng()%12;}
static void printsnap(){
  for(int r=0;r<snapN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",SNAP[r][c]);printf("\n");}
  fflush(stdout);}
int main(){
  T0=0;T0=now();
  for(int x=0;x<=K;x++){B1[x]=x;B2[x]=x*(x-1)/2;B3[x]=x*(x-1)*(x-2)/6;B4[x]=(int)((long)x*(x-1)*(x-2)*(x-3)/24);}
  const char*fn=getenv("BEST_CANDIDATE");FILE*f=fn?fopen(fn,"r"):NULL;
  if(!f){fprintf(stderr,"no candidate\n");return 0;}
  int v,nv=0;static int buf[MAXR*K];
  while(nv<MAXR*K&&fscanf(f,"%d",&v)==1)buf[nv++]=v;fclose(f);
  snapN=nv/K;for(int r=0;r<snapN;r++)for(int c=0;c<K;c++)SNAP[r][c]=(unsigned char)(buf[r*K+c]%3);
  long NS=B4[K];cnt=malloc(sizeof(uint16_t)*NS*81);sm=malloc(sizeof(uint32_t)*NS*81);
  const double LIMIT=3420.0;
  int order[MAXR],epoch=0,needOrder=1;
  while(now()<LIMIT-60){
    memcpy(A,SNAP,sizeof(A));
    if(needOrder){
      for(int r=0;r<snapN;r++)active[r]=1;recompute();
      fprintf(stderr,"[%.0f] N=%d holes=%d\n",now(),snapN,nh);
      if(nh>0){fprintf(stderr,"input not valid\n");break;}
      static long uq[MAXR];
      for(int r=0;r<snapN;r++){long s=0;for(int c=0;c<K;c++)s+=u[r][c];uq[r]=s/4;order[r]=r;}
      for(int i=1;i<snapN;i++){int x=order[i],j=i;while(j>0&&uq[order[j-1]]>uq[x]){order[j]=order[j-1];j--;}order[j]=x;}
      for(int i=0;i<8&&i<snapN;i++)fprintf(stderr," row %d uniq %ld\n",order[i],uq[order[i]]);
      needOrder=0;epoch=0;
    }
    int del=order[epoch%snapN];epoch++;
    for(int r=0;r<snapN;r++)active[r]=(r!=del);
    memcpy(A,SNAP,sizeof(A));recompute();
    memset(tabu,0,sizeof(tabu));iter=0;
    int besth=nh;double tstart=now(),budget=(epoch==1?1100:600);
    double T=1.0;
    fprintf(stderr,"[%.0f] epoch %d delete row %d holes %d\n",now(),epoch,del,nh);
    while(nh>0&&now()<LIMIT&&now()-tstart<budget){
      iter++;
      Hole h=holes[rng()%nh];
      int bestd=1<<30,br=-1,bc=0,bv=0,ties=0;
      for(int i=0;i<nact;i++){int r=act[i];int mm=0,t=-1;
        for(int s=0;s<4;s++)if(A[r][h.cols[s]]!=h.vals[s]){mm++;t=s;if(mm>1)break;}
        if(mm!=1)continue;
        int c=h.cols[t],vv=h.vals[t];int d=u[r][c]-gain(r,c,vv);
        if(tabu[r][c]>iter&&!(nh+d<besth))continue;
        if(d<bestd){bestd=d;br=r;bc=c;bv=vv;ties=1;}
        else if(d==bestd&&rng()%(++ties)==0){br=r;bc=c;bv=vv;}}
      if(br>=0&&(bestd<=0||(double)(rng()%1000000)/1e6<exp(-bestd/T))){apply(br,bc,bv);}
      else{
        int cr=-1,cc=0,cv=0,nt=0;
        if(rng()%10){
          for(int i=0;i<nact;i++){int r=act[i];int mm=0,ts[4];
            for(int s=0;s<4;s++)if(A[r][h.cols[s]]!=h.vals[s])ts[mm++]=s;
            if(mm!=2)continue;
            for(int q=0;q<2;q++){int c=h.cols[ts[q]];
              if(u[r][c]==0&&tabu[r][c]<=iter&&rng()%(++nt)==0){cr=r;cc=c;cv=h.vals[ts[q]];}}}
        }
        if(cr<0){for(int tries=0;tries<500;tries++){int r=act[rng()%nact],c=rng()%K;
            if(u[r][c]==0&&tabu[r][c]<=iter){cr=r;cc=c;cv=(A[r][c]+1+rng()%2)%3;break;}}}
        if(cr>=0)apply(cr,cc,cv);
      }
      if(nh<besth){besth=nh;fprintf(stderr,"[%.0f] it %lld holes %d\n",now(),iter,nh);}
      if(iter%2000==0)fprintf(stderr,"[%.0f] it %lld holes %d best %d\n",now(),iter,nh,besth);
    }
    if(nh==0){
      int n=0;for(int r=0;r<snapN;r++)if(active[r])memcpy(SNAP[n++],A[r],K);
      snapN=n;needOrder=1;
      fprintf(stderr,"[%.0f] SUCCESS N=%d\n",now(),snapN);
    }
  }
  printsnap();
  return 0;
}
