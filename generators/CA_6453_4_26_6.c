#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define K 26
#define V 6
#define VT 1296
#define MAXR 6600
static int N,NB,NO;
static unsigned char A[MAXR][K],BEST[MAXR][K],ORIG[MAXR][K],LBA[MAXR][K];
static uint32_t Bn[K+1][5]; static uint32_t NT;
static uint8_t *cnt; static uint16_t *X; static uint32_t *hpos;
static int32_t U[MAXR][K]; static int32_t U2[MAXR][K][K];
typedef struct{uint32_t id; uint8_t c[4],v[4];}Hole;
static Hole *H; static size_t nh=0,hcap=0;
static uint32_t *CL; static size_t clcap=0; static int cst[K+1];
static double t0;
static double now(){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+ts.tv_nsec*1e-9;}
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static inline void addhole(uint32_t id,int s0,int s1,int s2,int s3,int tu){
  if(nh==hcap){hcap=hcap?2*hcap:4096;H=realloc(H,hcap*sizeof(Hole));}
  Hole*h=&H[nh];h->id=id;h->c[0]=s0;h->c[1]=s1;h->c[2]=s2;h->c[3]=s3;
  h->v[0]=tu/216;h->v[1]=(tu/36)%6;h->v[2]=(tu/6)%6;h->v[3]=tu%6;
  hpos[id]=(uint32_t)(nh+1); nh++;
}
static inline void remhole(uint32_t id){
  size_t i=hpos[id]-1; hpos[id]=0; nh--; if(i!=nh){H[i]=H[nh];hpos[H[i].id]=(uint32_t)(i+1);}
}
static inline void uadd(int rr,int a,int b,int c,int d,int du){
  int32_t*u=U[rr]; u[a]+=du;u[b]+=du;u[c]+=du;u[d]+=du;
  U2[rr][a][b]+=du;U2[rr][a][c]+=du;U2[rr][a][d]+=du;U2[rr][b][c]+=du;U2[rr][b][d]+=du;U2[rr][c][d]+=du;
}
static void applyc(int r,int c,int nv){
  unsigned char*row=A[r]; int o=row[c]; int dv=nv-o; if(!dv)return;
  for(int z=1;z<K;z++){ if(z==c)continue;
    for(int y=0;y<z;y++){ if(y==c)continue;
      for(int x=0;x<y;x++){ if(x==c)continue;
        uint32_t R;int tu,w,s0,s1,s2,s3;
        if(c>z){R=x+Bn[y][2]+Bn[z][3]+Bn[c][4];tu=216*row[x]+36*row[y]+6*row[z]+o;w=1;s0=x;s1=y;s2=z;s3=c;}
        else if(c>y){R=x+Bn[y][2]+Bn[c][3]+Bn[z][4];tu=216*row[x]+36*row[y]+6*o+row[z];w=6;s0=x;s1=y;s2=c;s3=z;}
        else if(c>x){R=x+Bn[c][2]+Bn[y][3]+Bn[z][4];tu=216*row[x]+36*o+6*row[y]+row[z];w=36;s0=x;s1=c;s2=y;s3=z;}
        else{R=c+Bn[x][2]+Bn[y][3]+Bn[z][4];tu=216*o+36*row[x]+6*row[y]+row[z];w=216;s0=c;s1=x;s2=y;s3=z;}
        long i0=(long)R*VT+tu, i1=i0+(long)dv*w;
        uint8_t co=cnt[i0]--; X[i0]^=(uint16_t)r;
        if(co==1){uadd(r,s0,s1,s2,s3,-1);addhole((uint32_t)i0,s0,s1,s2,s3,tu);}
        else if(co==2){uadd(X[i0],s0,s1,s2,s3,1);}
        uint8_t cn=cnt[i1]++;
        if(cn==0){remhole((uint32_t)i1);uadd(r,s0,s1,s2,s3,1);}
        else if(cn==1){uadd(X[i1],s0,s1,s2,s3,-1);}
        X[i1]^=(uint16_t)r;
      }
    }
  }
  row[c]=nv;
}
static unsigned char colv[K][MAXR];
static void buildcol(unsigned char (*M)[K],int n){for(int r=0;r<n;r++)for(int c=0;c<K;c++)colv[c][r]=M[r][c];}
static void buildall(){
  memset(cnt,0,(size_t)NT*VT); memset(X,0,(size_t)NT*VT*2);
  for(size_t j=0;j<nh;j++)hpos[H[j].id]=0; nh=0;
  memset(U,0,sizeof U); memset(U2,0,sizeof U2);
  buildcol(A,N);
  static int q3[MAXR],q2[MAXR]; uint32_t T=0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++){
    for(int r=0;r<N;r++)q3[r]=6*colv[c][r]+colv[d][r];
    for(int b=1;b<c;b++){
      for(int r=0;r<N;r++)q2[r]=36*colv[b][r]+q3[r];
      for(int a=0;a<b;a++,T++){
        size_t base=(size_t)T*VT; uint8_t*p=cnt+base; uint16_t*px=X+base; unsigned char*ca=colv[a];
        for(int r=0;r<N;r++){int id=216*ca[r]+q2[r];p[id]++;px[id]^=(uint16_t)r;}
        for(int t=0;t<VT;t++){
          if(p[t]==1)uadd(px[t],a,b,c,d,1);
          else if(p[t]==0)addhole((uint32_t)(base+t),a,b,c,d,t);
        }
      }
    }
  }
}
static int verify(unsigned char (*M)[K],int n){
  memset(cnt,0,(size_t)NT*VT); buildcol(M,n);
  static int q3[MAXR],q2[MAXR]; uint32_t T=0; long miss=0;
  for(int r=0;r<n;r++)for(int c=0;c<K;c++)if(M[r][c]>=V)return 0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++){
    for(int r=0;r<n;r++)q3[r]=6*colv[c][r]+colv[d][r];
    for(int b=1;b<c;b++){
      for(int r=0;r<n;r++)q2[r]=36*colv[b][r]+q3[r];
      for(int a=0;a<b;a++,T++){
        uint8_t*p=cnt+(size_t)T*VT; unsigned char*ca=colv[a];
        for(int r=0;r<n;r++)p[216*ca[r]+q2[r]]=1;
        for(int t=0;t<VT;t++)if(!p[t])miss++;
      }
    }
  }
  return miss==0;
}
static void deleterow(int r0){
  int last=N-1; unsigned char*row=A[r0],*rl=A[last]; uint16_t rel=(uint16_t)(last^r0); uint32_t T=0;
  for(int d=3;d<K;d++)for(int c=2;c<d;c++)for(int b=1;b<c;b++)for(int a=0;a<b;a++,T++){
    size_t base=(size_t)T*VT; int tu=216*row[a]+36*row[b]+6*row[c]+row[d]; size_t i=base+tu;
    uint8_t co=cnt[i]--; X[i]^=(uint16_t)r0;
    if(co==1)addhole((uint32_t)i,a,b,c,d,tu);
    else if(co==2)uadd(X[i],a,b,c,d,1);
    if(r0!=last){int tl=216*rl[a]+36*rl[b]+6*rl[c]+rl[d];X[base+tl]^=rel;}
  }
  if(r0!=last){memcpy(A[r0],A[last],K);memcpy(U[r0],U[last],sizeof U[0]);memcpy(U2[r0],U2[last],sizeof U2[0]);}
  N--;
}
static long uniqv[MAXR]; static int ordr[MAXR];
static int cmpu(const void*a,const void*b){long x=uniqv[*(int*)a],y=uniqv[*(int*)b];return x<y?-1:x>y;}
static int pickrow(int rank){
  for(int r=0;r<N;r++){long s=0;for(int c=0;c<K;c++)s+=U[r][c];uniqv[r]=s;ordr[r]=r;}
  for(int r=N-1;r>0;r--){int j=rnd()%(r+1);int t=ordr[r];ordr[r]=ordr[j];ordr[j]=t;}
  qsort(ordr,N,sizeof(int),cmpu);
  return ordr[rank%N];
}
static void buildcl(){
  int cc[K+1]; memset(cc,0,sizeof cc);
  for(size_t j=0;j<nh;j++)for(int k=0;k<4;k++)cc[H[j].c[k]]++;
  cst[0]=0; for(int c=0;c<K;c++)cst[c+1]=cst[c]+cc[c];
  if(clcap<4*nh+4){clcap=8*nh+64;CL=realloc(CL,clcap*4);}
  int pos[K]; for(int c=0;c<K;c++)pos[c]=cst[c];
  for(size_t j=0;j<nh;j++)for(int k=0;k<4;k++)CL[pos[H[j].c[k]]++]=(uint32_t)((j<<2)|k);
}
static inline int evals(const unsigned char*row,int col,int val){
  int gg=0;
  for(int i=cst[col];i<cst[col+1];i++){uint32_t e=CL[i];const Hole*q=&H[e>>2];int p=e&3;
    if(q->v[p]!=val)continue; int ok=1;
    for(int k=0;k<4;k++)if(k!=p&&row[q->c[k]]!=q->v[k]){ok=0;break;}
    gg+=ok;
  }
  return gg;
}
static inline int evalp(const unsigned char*row,int c1,int v1,int c2,int v2){
  int gg=0;
  for(int i=cst[c1];i<cst[c1+1];i++){uint32_t e=CL[i];const Hole*q=&H[e>>2];int p=e&3;
    if(q->v[p]!=v1)continue; int ok=1;
    for(int k=0;k<4;k++){if(k==p)continue;int cc=q->c[k];int rv=(cc==c2)?v2:row[cc];if(rv!=q->v[k]){ok=0;break;}}
    gg+=ok;
  }
  for(int i=cst[c2];i<cst[c2+1];i++){uint32_t e=CL[i];const Hole*q=&H[e>>2];int p=e&3;
    if(q->v[p]!=v2)continue;
    if(q->c[0]==c1||q->c[1]==c1||q->c[2]==c1||q->c[3]==c1)continue;
    int ok=1;
    for(int k=0;k<4;k++)if(k!=p&&row[q->c[k]]!=q->v[k]){ok=0;break;}
    gg+=ok;
  }
  return gg;
}
static void out(unsigned char (*M)[K],int n){
  for(int r=0;r<n;r++){for(int c=0;c<K;c++){if(c)putchar(' ');putchar('0'+M[r][c]);}putchar('\n');}
  fflush(stdout);
}
static double ex[256];
static void settemp(double T){for(int k=0;k<256;k++)ex[k]=exp(-k/T);}
static void revert(){
  long diff=0;
  for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(A[r][c]!=LBA[r][c])diff++;
  if(diff==0)return;
  if(diff<=4000){for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(A[r][c]!=LBA[r][c])applyc(r,c,LBA[r][c]);}
  else{memcpy(A,LBA,(size_t)N*K);buildall();}
}
static int cr_[MAXR],cc1[MAXR],cv1[MAXR],cc2[MAXR],cv2[MAXR],cd[MAXR],ctp[MAXR]; static double cw[MAXR];
int main(){
  t0=now(); double TL=3350.0;
  rs^=(uint64_t)(now()*1e9); for(int i=0;i<10;i++)rnd();
  for(int n=0;n<=K;n++){Bn[n][0]=1;for(int k=1;k<5;k++)Bn[n][k]=(n==0)?0:Bn[n-1][k-1]+Bn[n-1][k];}
  NT=Bn[K][4];
  const char*fn=getenv("BEST_CANDIDATE"); FILE*f=fn?fopen(fn,"r"):NULL;
  if(!f){fprintf(stderr,"no candidate\n");return 0;}
  int v,i=0; while(fscanf(f,"%d",&v)==1){if(i/K>=MAXR)break;A[i/K][i%K]=(unsigned char)(((v%V)+V)%V);i++;}
  fclose(f); N=i/K; NO=N; memcpy(ORIG,A,(size_t)N*K); NB=N; memcpy(BEST,A,(size_t)N*K);
  cnt=malloc((size_t)NT*VT); X=malloc((size_t)NT*VT*2); hpos=calloc((size_t)NT*VT,4);
  if(!cnt||!X||!hpos){out(ORIG,NO);return 0;}
  buildall();
  fprintf(stderr,"N=%d holes=%zu t=%.1f\n",N,nh,now()-t0);
  int validseen=(nh==0);
  const double Thi=0.7,Tlo=0.12,CYC=60.0,STAG=500.0;
  double T=Thi; settemp(T);
  long lbn=(long)nh; memcpy(LBA,A,(size_t)N*K);
  double cyc0=now(),lastimp=now(),lastlog=now(); int tries=0,dirty=1; long iters=0,moves=0;
  while(1){
    if(nh==0){
      NB=N; memcpy(BEST,A,(size_t)N*K); validseen=1;
      fprintf(stderr,"VALID N=%d t=%.1f moves=%ld\n",N,now()-t0,moves);
      tries=0; int r0=pickrow(0); fprintf(stderr,"delete row %d uniq*4=%ld\n",r0,uniqv[r0]);
      deleterow(r0); dirty=1; lbn=(long)nh; memcpy(LBA,A,(size_t)N*K);
      cyc0=lastimp=now(); T=Thi; settemp(T); continue;
    }
    iters++;
    if((iters&63)==0){
      double tn=now(); if(tn-t0>TL)break;
      double fr=(tn-cyc0)/CYC;
      if(fr>=1.0){
        if((long)nh>lbn+2){revert();dirty=1;}
        cyc0=tn; fr=0;
      }
      double nT=Thi*pow(Tlo/Thi,fr); if(fabs(nT-T)>1e-3){T=nT;settemp(T);}
      if(validseen && tn-lastimp>STAG && TL-(tn-t0)>300){
        tries++; fprintf(stderr,"restart level N=%d try=%d lbn=%ld\n",NB-1,tries,lbn);
        memcpy(A,BEST,(size_t)NB*K); N=NB; buildall();
        int r0=pickrow(tries); deleterow(r0); dirty=1;
        lbn=(long)nh; memcpy(LBA,A,(size_t)N*K); cyc0=lastimp=now(); T=Thi; settemp(T); continue;
      }
      if(tn-lastlog>20){lastlog=tn;fprintf(stderr,"t=%.0f N=%d holes=%zu lbn=%ld T=%.3f moves=%ld iters=%ld\n",tn-t0,N,nh,lbn,T,moves,iters);}
    }
    if(dirty){buildcl();dirty=0;}
    Hole hh=H[rnd()%nh];
    int nc=0,dmin=0;
    int h0=hh.c[0],h1=hh.c[1],h2=hh.c[2],h3=hh.c[3];
    for(int r=0;r<N;r++){
      const unsigned char*row=A[r];
      int dif=(row[h0]!=hh.v[0])+(row[h1]!=hh.v[1])+(row[h2]!=hh.v[2])+(row[h3]!=hh.v[3]);
      if(dif>2||dif==0)continue;
      int p1=-1,p2=-1;
      for(int k=0;k<4;k++)if(row[hh.c[k]]!=hh.v[k]){if(p1<0)p1=k;else p2=k;}
      if(dif==1){
        int col=hh.c[p1],val=hh.v[p1];
        int d=U[r][col]-evals(row,col,val);
        cr_[nc]=r;ctp[nc]=1;cc1[nc]=col;cv1[nc]=val;cd[nc]=d;nc++; if(d<dmin)dmin=d;
      } else {
        int c1=hh.c[p1],v1=hh.v[p1],c2=hh.c[p2],v2=hh.v[p2];
        int d=U[r][c1]+U[r][c2]-U2[r][c1][c2]-evalp(row,c1,v1,c2,v2);
        cr_[nc]=r;ctp[nc]=2;cc1[nc]=c1;cv1[nc]=v1;cc2[nc]=c2;cv2[nc]=v2;cd[nc]=d;nc++; if(d<dmin)dmin=d;
      }
    }
    if(nc==0)continue;
    int kn=-dmin; if(kn>255)kn=255;
    double Z=ex[kn];
    for(int q=0;q<nc;q++){int k=cd[q]-dmin; if(k>255)k=255; cw[q]=ex[k]; Z+=cw[q];}
    double u=urand()*Z;
    if(u<ex[kn])continue;
    u-=ex[kn]; int sel=nc-1;
    for(int q=0;q<nc;q++){if(u<cw[q]){sel=q;break;}u-=cw[q];}
    applyc(cr_[sel],cc1[sel],cv1[sel]);
    if(ctp[sel]==2)applyc(cr_[sel],cc2[sel],cv2[sel]);
    moves++; dirty=1;
    if((long)nh<lbn){lbn=(long)nh;memcpy(LBA,A,(size_t)N*K);lastimp=now();}
  }
  fprintf(stderr,"final best N=%d, verifying t=%.1f\n",NB,now()-t0);
  if(NB<NO && verify(BEST,NB)){out(BEST,NB);}
  else out(ORIG,NO);
  return 0;
}
