#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define K 12
#define V 6
#define NS 495
#define NT 1296
static int N;
static unsigned char (*A)[K];
static unsigned short *cnt;
static int *mis,*pos,M=0;
static int scol[NS][4];
static int csid[K][165],cot[K][165][3],cw[K][165],cwo[K][165][3];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(int r,int s){unsigned char*x=A[r];return x[scol[s][0]]*216+x[scol[s][1]]*36+x[scol[s][2]]*6+x[scol[s][3]];}
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  for(int j=0;j<165;j++){
    int base=A[r][cot[c][j][0]]*cwo[c][j][0]+A[r][cot[c][j][1]]*cwo[c][j][1]+A[r][cot[c][j][2]]*cwo[c][j][2];
    int g=csid[c][j]*NT;
    int o=g+base+a*cw[c][j],n=g+base+b*cw[c][j];
    if(--cnt[o]==0)addm(o);
    if(cnt[n]++==0)delm(n);
  }
  A[r][c]=b;
}
static inline int delta(int r,int c,int b){
  int a=A[r][c],d=0;
  for(int j=0;j<165;j++){
    int base=A[r][cot[c][j][0]]*cwo[c][j][0]+A[r][cot[c][j][1]]*cwo[c][j][1]+A[r][cot[c][j][2]]*cwo[c][j][2];
    int g=csid[c][j]*NT;
    if(cnt[g+base+a*cw[c][j]]==1)d++;
    if(cnt[g+base+b*cw[c][j]]==0)d--;
  }
  return d;
}
int main(){
  double t0=now(),TL=3480;
  const char*fn=getenv("BEST_CANDIDATE");
  int cap=5000; A=malloc(sizeof(*A)*cap); N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;if(N>=cap)break;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 1;}
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={216,36,6,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){int c=scol[s][i];int j=cc[c]++;csid[c][j]=s;cw[c][j]=wts[i];int q=0;for(int k=0;k<4;k++)if(k!=i){cot[c][j][q]=scol[s][k];cwo[c][j][q]=wts[k];q++;}}
  cnt=calloc(NS*NT,2);mis=malloc(sizeof(int)*NS*NT);pos=malloc(sizeof(int)*NS*NT);
  for(int i=0;i<NS*NT;i++)pos[i]=-1;
  for(int r=0;r<N;r++)for(s=0;s<NS;s++)cnt[s*NT+tidx(r,s)]++;
  for(int i=0;i<NS*NT;i++)if(cnt[i]==0)addm(i);
  int bestN=-1;unsigned char (*B)[K]=malloc(sizeof(*B)*cap);
  if(M==0){bestN=N;memcpy(B,A,sizeof(*A)*N);}
  int *cr=malloc(sizeof(int)*cap),*ccol=malloc(sizeof(int)*cap),*cb=malloc(sizeof(int)*cap),*cd=malloc(sizeof(int)*cap);
  double T=0.35; long long it=0;
  while(1){
    if((it&1023)==0&&now()-t0>TL)break;
    it++;
    if(M==0){
      if(bestN<0||N<bestN){bestN=N;memcpy(B,A,sizeof(*A)*N);fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);}
      int br=-1,bu=1<<30;
      for(int r=0;r<N;r++){int u=0;for(s=0;s<NS;s++)if(cnt[s*NT+tidx(r,s)]==1)u++; if(u<bu){bu=u;br=r;}}
      for(s=0;s<NS;s++){int g=s*NT+tidx(br,s);if(--cnt[g]==0)addm(g);}
      N--; if(br!=N)memcpy(A[br],A[N],K);
      fprintf(stderr,"deleted row, N=%d missing=%d\n",N,M);
      continue;
    }
    int g=mis[rnd()%M]; s=g/NT; int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int *sc=scol[s]; int nc=0;
    for(int r=0;r<N;r++){
      int mm=0,mc=0;
      for(int i=0;i<4;i++)if(A[r][sc[i]]!=tv[i]){mm++;mc=i;if(mm>1)break;}
      if(mm==1){cr[nc]=r;ccol[nc]=sc[mc];cb[nc]=tv[mc];nc++;}
    }
    if(nc==0){
      int r=rnd()%N;
      for(int i=0;i<4;i++)change(r,sc[i],tv[i]);
      continue;
    }
    int bd=1<<30,bi=-1,ties=0;
    for(int i=0;i<nc;i++){int d=delta(cr[i],ccol[i],cb[i]);cd[i]=d;
      if(d<bd){bd=d;bi=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)bi=i;}}
    int pick=bi;
    if(rnd()%100==0)pick=rnd()%nc;
    int d=cd[pick];
    if(d<=0||(double)(rnd()%1000000)/1e6<exp(-d/T))change(cr[pick],ccol[pick],cb[pick]);
  }
  if(bestN<0){bestN=N;memcpy(B,A,sizeof(*A)*N);}
  fprintf(stderr,"final N=%d iters=%lld\n",bestN,it);
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",B[r][c]);printf("\n");}
  return 0;
}
