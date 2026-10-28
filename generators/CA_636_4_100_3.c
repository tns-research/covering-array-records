#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define K 100
#define MAXR 700
static int N;
static unsigned char A[MAXR][K];
static unsigned char BEST[MAXR][K]; static int NB;
static uint32_t Bn[101][5];
static uint8_t *cnt;
static uint32_t NT;
static uint32_t *hl; static size_t hlen=0, hcap=0;
static long nunc=0;
static double t0;
static double now(){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+ts.tv_nsec*1e-9;}
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static inline void push(uint32_t x){ if(hlen==hcap){hcap=hcap?hcap*2:1024;hl=realloc(hl,hcap*4);} hl[hlen++]=x; }

/* enumerate 4-sets containing column c; macro body gets rank R, tuple tu (row r), weight w */
static long evalc(int r,int c,int nv,long cutoff){
  unsigned char *row=A[r]; int o=row[c]; int dv=nv-o;
  long loss=0,gain=0;
  for(int x=0;x<K;x++){ if(x==c)continue;
    for(int y=x+1;y<K;y++){ if(y==c)continue;
      for(int z=y+1;z<K;z++){ if(z==c)continue;
        uint32_t R; int tu,w;
        if(c>z){R=Bn[x][1]+Bn[y][2]+Bn[z][3]+Bn[c][4];tu=27*row[x]+9*row[y]+3*row[z]+o;w=1;}
        else if(c>y){R=Bn[x][1]+Bn[y][2]+Bn[c][3]+Bn[z][4];tu=27*row[x]+9*row[y]+3*o+row[z];w=3;}
        else if(c>x){R=Bn[x][1]+Bn[c][2]+Bn[y][3]+Bn[z][4];tu=27*row[x]+9*o+3*row[y]+row[z];w=9;}
        else {R=Bn[c][1]+Bn[x][2]+Bn[y][3]+Bn[z][4];tu=27*o+9*row[x]+3*row[y]+row[z];w=27;}
        uint8_t *p=cnt+(size_t)R*81;
        if(p[tu]==1)loss++;
        if(p[tu+dv*w]==0)gain++;
      }
    }
    if(loss-gain-nunc>cutoff) return 1L<<40;
  }
  return loss-gain;
}
static void applyc(int r,int c,int nv){
  unsigned char *row=A[r]; int o=row[c]; int dv=nv-o; if(!dv)return;
  for(int x=0;x<K;x++){ if(x==c)continue;
    for(int y=x+1;y<K;y++){ if(y==c)continue;
      for(int z=y+1;z<K;z++){ if(z==c)continue;
        uint32_t R; int tu,w;
        if(c>z){R=Bn[x][1]+Bn[y][2]+Bn[z][3]+Bn[c][4];tu=27*row[x]+9*row[y]+3*row[z]+o;w=1;}
        else if(c>y){R=Bn[x][1]+Bn[y][2]+Bn[c][3]+Bn[z][4];tu=27*row[x]+9*row[y]+3*o+row[z];w=3;}
        else if(c>x){R=Bn[x][1]+Bn[c][2]+Bn[y][3]+Bn[z][4];tu=27*row[x]+9*o+3*row[y]+row[z];w=9;}
        else {R=Bn[c][1]+Bn[x][2]+Bn[y][3]+Bn[z][4];tu=27*o+9*row[x]+3*row[y]+row[z];w=27;}
        uint8_t *p=cnt+(size_t)R*81;
        if(--p[tu]==0){nunc++;push(R*81+tu);}
        if(p[tu+dv*w]++==0)nunc--;
      }
    }
  }
  row[c]=nv;
}
static unsigned char colv[K][MAXR];
static void buildcol(){for(int r=0;r<N;r++)for(int c=0;c<K;c++)colv[c][r]=A[r][c];}
static void buildcounts(){
  memset(cnt,0,(size_t)NT*81); buildcol();
  static int q3[MAXR],q2[MAXR];
  uint32_t T=0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++){
    for(int r=0;r<N;r++)q3[r]=3*colv[c][r]+colv[d][r];
    for(int b=1;b<c;b++){
      for(int r=0;r<N;r++)q2[r]=9*colv[b][r]+q3[r];
      for(int a=0;a<b;a++,T++){
        uint8_t *p=cnt+(size_t)T*81; unsigned char *ca=colv[a];
        for(int r=0;r<N;r++)p[27*ca[r]+q2[r]]++;
      }
    }
  }
  nunc=0; hlen=0;
  for(size_t i=0;i<(size_t)NT*81;i++) if(!cnt[i]){nunc++;push((uint32_t)i);}
}
static long uniq[MAXR];
static void computeuniq(){
  buildcol(); memset(uniq,0,sizeof(uniq));
  static int q3[MAXR],q2[MAXR]; uint32_t T=0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++){
    for(int r=0;r<N;r++)q3[r]=3*colv[c][r]+colv[d][r];
    for(int b=1;b<c;b++){
      for(int r=0;r<N;r++)q2[r]=9*colv[b][r]+q3[r];
      for(int a=0;a<b;a++,T++){
        uint8_t *p=cnt+(size_t)T*81; unsigned char *ca=colv[a];
        for(int r=0;r<N;r++)if(p[27*ca[r]+q2[r]]==1)uniq[r]++;
      }
    }
  }
}
static void deleterow(int r0){
  unsigned char *row=A[r0]; uint32_t T=0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++)for(int b=1;b<c;b++)for(int a=0;a<b;a++,T++){
    int tu=27*row[a]+9*row[b]+3*row[c]+row[d];
    uint8_t *p=cnt+(size_t)T*81;
    if(--p[tu]==0){nunc++;push(T*81+tu);}
  }
  if(r0!=N-1)memcpy(A[r0],A[N-1],K);
  N--;
}
static void unrank(uint32_t T,int s[4]){
  for(int k=4;k>=1;k--){int n=k-1; while(n+1<K && Bn[n+1][k]<=T)n++; s[k-1]=n; T-=Bn[n][k];}
}
static void savebest(){NB=N;memcpy(BEST,A,sizeof(A[0])*N);}
static void output(){
  for(int r=0;r<NB;r++){for(int c=0;c<K;c++){if(c)putchar(' ');putchar('0'+BEST[r][c]);}putchar('\n');}
  fflush(stdout);
}
static int tabu[MAXR][K];
int main(){
  t0=now(); double TL=3420.0;
  for(int n=0;n<=100;n++){Bn[n][0]=1;for(int k=1;k<5;k++)Bn[n][k]=(n==0)?0:Bn[n-1][k-1]+Bn[n-1][k];}
  NT=Bn[100][4];
  const char *fn=getenv("BEST_CANDIDATE"); FILE *f=fn?fopen(fn,"r"):NULL;
  if(!f){fprintf(stderr,"no candidate\n");return 0;}
  int v,i=0; while(fscanf(f,"%d",&v)==1){ if(i/K>=MAXR)break; A[i/K][i%K]=(unsigned char)(v%3); i++; }
  fclose(f); N=i/K;
  savebest();
  cnt=malloc((size_t)NT*81);
  buildcounts();
  fprintf(stderr,"N=%d uncovered=%ld t=%.1f\n",N,nunc,now()-t0);
  if(nunc!=0){ NB=0; } /* input invalid: don't trust as best */
  int validInput = (nunc==0);
  if(!validInput){ savebest(); } /* still output it at worst */
  long iter=0; long bestn=1L<<60; double temp=1.5; long lastimp=0;
  while(now()-t0<TL){
    if(nunc==0){
      savebest(); fprintf(stderr,"VALID N=%d t=%.1f\n",N,now()-t0);
      computeuniq(); int br=0; for(int r=1;r<N;r++)if(uniq[r]<uniq[br])br=r;
      fprintf(stderr,"delete row %d uniq=%ld\n",br,uniq[br]);
      deleterow(br); bestn=nunc; lastimp=iter;
      memset(tabu,0,sizeof(tabu));
      continue;
    }
    iter++;
    if(hlen>4*(size_t)nunc+4096){size_t j=0;for(size_t k=0;k<hlen;k++)if(cnt[hl[k]]==0)hl[j++]=hl[k];hlen=j;}
    uint32_t h;
    for(;;){size_t k=rnd()%hlen; h=hl[k]; if(cnt[h]==0)break; hl[k]=hl[--hlen];}
    uint32_t T=h/81; int t=h%81; int s[4],vv[4];
    unrank(T,s); vv[0]=t/27;vv[1]=(t/9)%3;vv[2]=(t/3)%3;vv[3]=t%3;
    static int cr[MAXR],cc[MAXR],cv[MAXR]; int nc=0;
    static int d2[MAXR]; int n2=0;
    for(int r=0;r<N;r++){int dif=0,pc=-1;
      for(int k=0;k<4;k++)if(A[r][s[k]]!=vv[k]){dif++;pc=k;}
      if(dif==1){ if(tabu[r][s[pc]]<=iter){cr[nc]=r;cc[nc]=s[pc];cv[nc]=vv[pc];nc++;} }
      else if(dif==2)d2[n2++]=r;
    }
    if(nc==0){
      if(n2==0)continue;
      int r=d2[rnd()%n2];
      for(int k=0;k<4;k++)if(A[r][s[k]]!=vv[k]){applyc(r,s[k],vv[k]);tabu[r][s[k]]=iter+10;}
    } else {
      /* shuffle */
      for(int k=nc-1;k>0;k--){int j=rnd()%(k+1);int a;a=cr[k];cr[k]=cr[j];cr[j]=a;a=cc[k];cc[k]=cc[j];cc[j]=a;a=cv[k];cv[k]=cv[j];cv[j]=a;}
      long best=1L<<40; int bi=-1;
      for(int k=0;k<nc;k++){long dlt=evalc(cr[k],cc[k],cv[k],best);
        if(dlt<best){best=dlt;bi=k;}}
      if(bi>=0 && (best<=0 || urand()<exp(-(double)best/temp))){
        applyc(cr[bi],cc[bi],cv[bi]); tabu[cr[bi]][cc[bi]]=iter+10+(int)(rnd()%10);
      }
    }
    if(nunc<bestn){bestn=nunc;lastimp=iter;}
    if(iter-lastimp>3000){temp=(temp<4.0)?temp*1.2:1.0; lastimp=iter;}
    if(iter%500==0)fprintf(stderr,"it=%ld N=%d unc=%ld best=%ld temp=%.2f t=%.0f\n",iter,N,nunc,bestn,temp,now()-t0);
  }
  output();
  return 0;
}
