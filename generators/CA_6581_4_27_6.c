#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define K 27
#define V 6
#define NS 17550
#define NT 1296
#define NG (NS*NT)
#define CAP 7000
#define PC 2600
static int N;
static unsigned char A[CAP][K], B[CAP][K], I0[CAP][K];
static unsigned short cnt[NG];
static int w[NG];
static int mis[NG],pos[NG],M=0;
static int hv[NG],hpos[NG],H=0;
static long long tb[CAP][K];
static int scol[NS][4];
static int csid[K][PC],cot[K][PC][3],cw[K][PC],cwo[K][PC][3];
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(int r,int s){unsigned char*x=A[r];return x[scol[s][0]]*216+x[scol[s][1]]*36+x[scol[s][2]]*6+x[scol[s][3]];}
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  unsigned char*x=A[r];
  for(int j=0;j<PC;j++){
    int base=x[cot[c][j][0]]*cwo[c][j][0]+x[cot[c][j][1]]*cwo[c][j][1]+x[cot[c][j][2]]*cwo[c][j][2];
    int g=csid[c][j]*NT;
    int o=g+base+a*cw[c][j],n=g+base+b*cw[c][j];
    if(--cnt[o]==0)addm(o);
    if(cnt[n]++==0)delm(n);
  }
  A[r][c]=b;
}
static inline int delta(int r,int c,int b,int*du){
  int a=A[r][c],d=0,u=0;
  unsigned char*x=A[r];
  for(int j=0;j<PC;j++){
    int base=x[cot[c][j][0]]*cwo[c][j][0]+x[cot[c][j][1]]*cwo[c][j][1]+x[cot[c][j][2]]*cwo[c][j][2];
    int g=csid[c][j]*NT;
    int o=g+base+a*cw[c][j],n=g+base+b*cw[c][j];
    if(cnt[o]==1){d+=w[o];u++;}
    if(cnt[n]==0){d-=w[n];u--;}
  }
  *du=u;return d;
}
static void rebuild(){
  memset(cnt,0,sizeof cnt);
  for(int i=0;i<NG;i++)pos[i]=-1;
  M=0;
  for(int r=0;r<N;r++)for(int s=0;s<NS;s++)cnt[s*NT+tidx(r,s)]++;
  for(int i=0;i<NG;i++)if(cnt[i]==0)addm(i);
}
static void resetw(){
  for(int i=0;i<H;i++){w[hv[i]]=1;hpos[hv[i]]=-1;}
  H=0;
}
static void decayw(int num,int den){
  for(int i=0;i<H;){
    int g=hv[i];
    w[g]=1+(int)((long long)(w[g]-1)*num/den);
    if(w[g]<=1){w[g]=1;hpos[g]=-1;int l=hv[--H];hv[i]=l;hpos[l]=i;}
    else i++;
  }
}
static int uniq(int r){int u=0;for(int s=0;s<NS;s++)if(cnt[s*NT+tidx(r,s)]==1)u++;return u;}
static void delrow(int br){
  for(int s=0;s<NS;s++){int g=s*NT+tidx(br,s);if(--cnt[g]==0)addm(g);}
  N--; if(br!=N)memcpy(A[br],A[N],K);
}
static int pickrow(int rank){
  static int us[CAP],ord[CAP];
  for(int r=0;r<N;r++){us[r]=uniq(r);ord[r]=r;}
  for(int i=N-1;i>0;i--){int j=rnd()%(i+1);int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  for(int i=0;i<=rank&&i<N;i++){
    int bi=i;for(int j=i+1;j<N;j++)if(us[ord[j]]<us[ord[bi]])bi=j;
    int t=ord[i];ord[i]=ord[bi];ord[bi]=t;
  }
  int r=ord[rank<N?rank:0];
  fprintf(stderr,"delete row %d unique=%d\n",r,us[r]);
  return r;
}
static void output(int n){
  for(int r=0;r<n;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",B[r][c]);printf("\n");}
  fflush(stdout);
}
int main(){
  double t0=now(),TL=3440,STALL=500;
  const char*fn=getenv("BEST_CANDIDATE");
  N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;if(N>=CAP)break;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 1;}
  int N0=N; memcpy(I0,A,sizeof(A[0])*N);
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={216,36,6,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){int c=scol[s][i];int j=cc[c]++;csid[c][j]=s;cw[c][j]=wts[i];int q=0;for(int k=0;k<4;k++)if(k!=i){cot[c][j][q]=scol[s][k];cwo[c][j][q]=wts[k];q++;}}
  for(int i=0;i<NG;i++){w[i]=1;hpos[i]=-1;}
  rebuild();
  fprintf(stderr,"input N=%d missing=%d\n",N,M);
  int bestN=-1;
  if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
  static int cr[CAP],ccol[CAP],cb[CAP];
  long long it=0,bumps=0; int bestMs=M; double tlast=now(); int restarts=0;
  memset(tb,0,sizeof tb);
  while(1){
    if((it&63)==0){
      double tn=now();
      if(tn-t0>TL)break;
      if(bestN>0 && tn-tlast>STALL && M>0){
        restarts++;
        N=bestN; memcpy(A,B,sizeof(A[0])*N); rebuild(); resetw();
        memset(tb,0,sizeof tb);
        int br=pickrow((int)(rnd()%(3+restarts)));
        delrow(br); bestMs=M; tlast=now();
        fprintf(stderr,"restart %d N=%d missing=%d t=%.0f\n",restarts,N,M,tn-t0);
        continue;
      }
      if((it&((1<<16)-1))==0)fprintf(stderr,"it=%lld N=%d M=%d bestMs=%d H=%d t=%.0f\n",it,N,M,bestMs,H,tn-t0);
    }
    it++;
    if(M==0){
      if(bestN<0||N<bestN){bestN=N;memcpy(B,A,sizeof(A[0])*N);fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);}
      int br=pickrow(0);
      delrow(br);
      memset(tb,0,sizeof tb);
      decayw(1,2);
      bestMs=M; tlast=now();
      fprintf(stderr,"deleted row, N=%d missing=%d\n",N,M);
      continue;
    }
    int g=mis[rnd()%M]; s=g/NT; int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int *sc=scol[s]; int nc=0;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];
      int mm=0,mc=0;
      if(x[sc[0]]!=tv[0]){mm++;mc=0;}
      if(x[sc[1]]!=tv[1]){if(++mm>1)continue;mc=1;}
      if(x[sc[2]]!=tv[2]){if(++mm>1)continue;mc=2;}
      if(x[sc[3]]!=tv[3]){if(++mm>1)continue;mc=3;}
      if(mm==1){cr[nc]=r;ccol[nc]=sc[mc];cb[nc]=tv[mc];nc++;}
    }
    if(nc==0){
      int r=rnd()%N;
      for(int i=0;i<4;i++){change(r,sc[i],tv[i]);tb[r][sc[i]]=it+10;}
      continue;
    }
    int pick=-1,bd=1<<30,ties=0;
    if(rnd()%1000<3){
      pick=rnd()%nc; int du; bd=delta(cr[pick],ccol[pick],cb[pick],&du);
    } else {
      for(int i=0;i<nc;i++){
        int du;int d=delta(cr[i],ccol[i],cb[i],&du);
        int tabu=tb[cr[i]][ccol[i]]>it;
        if(tabu && !(M+du<bestMs))continue;
        if(d<bd){bd=d;pick=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)pick=i;}
      }
      if(pick<0){pick=rnd()%nc;int du;bd=delta(cr[pick],ccol[pick],cb[pick],&du);}
    }
    int apply=0;
    if(bd<=0)apply=1;
    else{
      for(int i=0;i<M;i++){int q=mis[i];if(w[q]==1){hpos[q]=H;hv[H++]=q;}w[q]++;}
      bumps++;
      if(bumps%2000==0)decayw(4,5);
      if(rnd()%100<10)apply=1;
    }
    if(apply){
      int r=cr[pick],c=ccol[pick];
      change(r,c,cb[pick]);
      tb[r][c]=it+4+(long long)(rnd()%10);
      if(M<bestMs){bestMs=M;tlast=now();}
    }
  }
  if(bestN<0){
    if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
    else {bestN=N0;memcpy(B,I0,sizeof(A[0])*N0);}
  }
  fprintf(stderr,"final N=%d iters=%lld\n",bestN,it);
  output(bestN);
  return 0;
}
