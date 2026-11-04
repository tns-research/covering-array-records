#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <sys/mman.h>
#define TEMP 0.35
#define DRIFT 30
#define T 4
#define K 34
#define V 6
_Static_assert(T==4,"strength 4 code");
#define W2 (V*V)
#define W3 (V*V*V)
#define NT (V*V*V*V)
#define NS (K*(K-1)*(K-2)*(K-3)/24)
#define PC ((K-1)*(K-2)*(K-3)/6)
#define NG (NS*NT)
#define MAXN 7000
#define HMAX (1<<21)
#define TL 3450.0
#define NC 48
typedef struct{int gb;unsigned char o0,o1,o2,pad;unsigned short w0,w1,w2,wc;}E;
static volatile sig_atomic_t stopf=0;
static void onterm(int s){(void)s;stopf=1;}
static unsigned short *cnt;static int *pos;
static int hl[HMAX],holes;
static int scol[NS][4];static E ce[K][PC];
static unsigned char A[MAXN][K],Best[MAXN][K],Inp[MAXN][K];static int N,bestN,inN;
static long long tabu[MAXN][K];
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int tidx(const unsigned char*x,int s){int*c=scol[s];return s*NT+x[c[0]]*W3+x[c[1]]*W2+x[c[2]]*V+x[c[3]];}
static inline void addh(int g){if(holes<HMAX){pos[g]=holes;hl[holes++]=g;}}
static inline void delh(int g){int p=pos[g];if(p<0)return;int l=hl[--holes];hl[p]=l;pos[l]=p;pos[g]=-1;}
static int verify(unsigned char(*X)[K],int n){
  static unsigned char seen[NT];
  for(int s=0;s<NS;s++){memset(seen,0,NT);int d=0;
    for(int r=0;r<n;r++){int id=tidx(X[r],s)-s*NT;if(!seen[id]){seen[id]=1;d++;}}
    if(d!=NT)return 0;}
  return 1;
}
static void rebuild(void){
  memset(cnt,0,(size_t)NG*sizeof(*cnt));
  for(int r=0;r<N;r++)for(int s=0;s<NS;s++)cnt[tidx(A[r],s)]++;
  holes=0;for(int g=0;g<NG;g++){pos[g]=-1;if(cnt[g]==0)addh(g);}
}
static inline int udelta(const unsigned char*x,int c,int b){
  int a=x[c],d=0;E*ec=ce[c];
  for(int j=0;j<PC;j++){E*e=&ec[j];
    int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    d+=(cnt[base+a*e->wc]==1)-(cnt[base+b*e->wc]==0);}
  return d;
}
static inline int udelta2(const unsigned char*x,int c,int b,int cut){
  if(holes>4096)return udelta(x,c,b);
  int gain=0;
  for(int i=0;i<holes;i++){int g=hl[i],s=g/NT,tv=g%NT;int*cl=scol[s];
    int vv[4]={tv/W3,(tv/W2)%V,(tv/V)%V,tv%V};int hit=0,ok=1;
    for(int k=0;k<4;k++){int q=cl[k];int xv;if(q==c){hit=1;xv=b;}else xv=x[q];if(xv!=vv[k]){ok=0;break;}}
    if(ok&&hit)gain++;}
  int a=x[c],loss=0,lim=cut+gain;E*ec=ce[c];
  for(int j=0;j<PC;j++){E*e=&ec[j];
    int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    if(cnt[base+a*e->wc]==1){if(++loss>lim)return loss-gain;}}
  return loss-gain;
}
static unsigned char Ck[MAXN][K];
static inline void apply(unsigned char*x,int c,int b){
  int a=x[c];E*ec=ce[c];
  for(int j=0;j<PC;j++){E*e=&ec[j];
    int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    int ga=base+a*e->wc,gb=base+b*e->wc;
    if(--cnt[ga]==0)addh(ga);
    if(cnt[gb]++==0)delh(gb);}
  x[c]=(unsigned char)b;
}
static void restore_ck(void){
  for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(A[r][c]!=Ck[r][c])apply(A[r],c,Ck[r][c]);
}
static void delete_row(int r){
  for(int s=0;s<NS;s++){int g=tidx(A[r],s);if(--cnt[g]==0)addh(g);}
  N--;if(r!=N){memcpy(A[r],A[N],K);memcpy(tabu[r],tabu[N],sizeof tabu[0]);}
}
static int pick_weak(int rank){
  static int u[MAXN];int mn=1<<30;
  for(int r=0;r<N;r++){int c=0;
    for(int s=0;s<NS;s++)if(cnt[tidx(A[r],s)]==1){c++;if(rank==0&&c>mn)break;}
    if(c<mn)mn=c;
    u[r]=c*1024+(int)(rnd()&1023);}
  static char used[MAXN];memset(used,0,N);int pick=0;
  for(int k=0;k<=rank;k++){int b=-1;for(int r=0;r<N;r++)if(!used[r]&&(b<0||u[r]<u[b]))b=r;if(b<0)break;used[b]=1;pick=b;}
  fprintf(stderr,"delete row %d unique=%d (rank %d) N->%d\n",pick,u[pick]/1024,rank,N-1);
  return pick;
}
static void output(void){
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",Best[r][c]);printf("\n");}
  fflush(stdout);
}
int main(void){
  double t0=now(),hard=t0+TL;
  signal(SIGTERM,onterm);signal(SIGINT,onterm);
  {size_t al=(size_t)1<<21;
   size_t sc=((size_t)NG*sizeof(unsigned short)+al-1)/al*al,sp=((size_t)NG*sizeof(int)+al-1)/al*al;
   void*p1=0,*p2=0;
   if(posix_memalign(&p1,al,sc)||posix_memalign(&p2,al,sp)){fprintf(stderr,"oom\n");return 1;}
#ifdef MADV_HUGEPAGE
   madvise(p1,sc,MADV_HUGEPAGE);madvise(p2,sp,MADV_HUGEPAGE);
#endif
   cnt=(unsigned short*)p1;pos=(int*)p2;}
  {struct timespec ts;clock_gettime(CLOCK_REALTIME,&ts);rs^=(unsigned long long)ts.tv_nsec*0x2545F4914F6CDD1DULL;if(!rs)rs=1;for(int i=0;i<16;i++)rnd();}
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  int cc[K];memset(cc,0,sizeof cc);int wts[4]={W3,W2,V,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){
    int c=scol[s][i],j=cc[c]++;E*e=&ce[c][j];
    e->gb=s*NT;e->wc=(unsigned short)wts[i];
    int oc[3],ow[3],q=0;for(int k=0;k<4;k++)if(k!=i){oc[q]=scol[s][k];ow[q]=wts[k];q++;}
    e->o0=oc[0];e->o1=oc[1];e->o2=oc[2];e->w0=ow[0];e->w1=ow[1];e->w2=ow[2];}
  const char*fn=getenv("BEST_CANDIDATE");N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(N<MAXN&&fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;}}fclose(f);}}
  if(N==0||!verify(A,N)){fprintf(stderr,"no valid input (%d rows)\n",N);bestN=N;memcpy(Best,A,(size_t)N*K);output();return 1;}
  bestN=N;memcpy(Best,A,(size_t)N*K);fprintf(stderr,"input valid N=%d\n",N);
  inN=N;memcpy(Inp,A,(size_t)N*K);
  fprintf(stderr,"PROGRESS %d\n",-N);
  rebuild();
  int rank=0;double LB=300,limp=now();long long step=0;int lbest=1<<30;
  delete_row(pick_weak(0));
  static int cr[NC],ccol[NC],cv[NC],cc2[NC],cv2[NC];
  while(!stopf){
    if(holes==0){
      bestN=N;memcpy(Best,A,(size_t)N*K);fprintf(stderr,"NEW BEST N=%d t=%.0f\nPROGRESS %d\n",N,now()-t0,-N);
      rank=0;LB=300;limp=now();lbest=1<<30;
      if(now()>hard)break;
      delete_row(pick_weak(0));continue;
    }
    step++;
    if((step&255)==0){
      double tn=now();if(tn>hard)break;
      if((step&65535)==0)fprintf(stderr,"N=%d holes=%d lbest=%d step=%lld t=%.0f\n",N,holes,lbest,step,tn-t0);
      if(tn-limp>LB){
        N=bestN;memcpy(A,Best,(size_t)N*K);memset(tabu,0,sizeof tabu);rebuild();
        rank=(rank+1)%8;if(rank==0)LB*=1.5;
        limp=now();lbest=1<<30;delete_row(pick_weak(rank));continue;
      }
    }
    if(holes<lbest){lbest=holes;limp=now();memcpy(Ck,A,(size_t)N*K);}
    else if(holes>lbest+DRIFT){restore_ck();memset(tabu,0,sizeof tabu);continue;}
    int h=hl[rnd()%holes];int ss=h/NT,tv=h%NT;
    int val[4]={tv/W3,(tv/W2)%V,(tv/V)%V,tv%V};int*col=scol[ss];
    int nc=0,seen=0;
    int use2=(holes<=16)&&(rnd()&1);int lim=use2?2:1;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int mm=0,m1=0,m2=-1;
      for(int i=0;i<4;i++)if(x[col[i]]!=val[i]){if(mm==0)m1=i;else m2=i;mm++;if(mm>lim)break;}
      if(mm<1||mm>lim)continue;
      seen++;
      int slot=nc<NC?nc++:(int)(rnd()%seen);
      if(slot<NC){cr[slot]=r;ccol[slot]=col[m1];cv[slot]=val[m1];
        if(mm==2){cc2[slot]=col[m2];cv2[slot]=val[m2];}else cc2[slot]=-1;}
    }
    if(nc==0){
      int r=(int)(rnd()%N),i=(int)(rnd()%4);
      if(A[r][col[i]]!=val[i])apply(A[r],col[i],val[i]);
      continue;
    }
    int bi=-1,bd=1<<30,ties=0;
    for(int k=0;k<nc;k++){
      if(tabu[cr[k]][ccol[k]]>step)continue;
      int d;
      if(cc2[k]<0)d=udelta2(A[cr[k]],ccol[k],cv[k],bd);
      else{
        if(tabu[cr[k]][cc2[k]]>step)continue;
        int r=cr[k],a=A[r][ccol[k]],h0=holes;
        apply(A[r],ccol[k],cv[k]);
        int d1=holes-h0;
        d=d1+udelta2(A[r],cc2[k],cv2[k],bd>=(1<<29)?bd:bd-d1);
        apply(A[r],ccol[k],a);
      }
      if(d<bd){bd=d;bi=k;ties=1;}else if(d==bd&&(rnd()%(++ties))==0)bi=k;
    }
    if(bi<0)continue;
    int acc=(bd<=0)||((double)(rnd()%1000000)<exp(-bd/TEMP)*1e6);
    if(acc){
      apply(A[cr[bi]],ccol[bi],cv[bi]);
      tabu[cr[bi]][ccol[bi]]=step+10+(long long)(rnd()%20);
      if(cc2[bi]>=0){
        apply(A[cr[bi]],cc2[bi],cv2[bi]);
        tabu[cr[bi]][cc2[bi]]=step+10+(long long)(rnd()%20);
      }
    }else{
      for(int z=0;z<48;z++){
        int r=(int)(rnd()%N),c=(int)(rnd()%K),b=(int)(rnd()%V);
        if(b==A[r][c])continue;
        if(udelta2(A[r],c,b,0)<=0)apply(A[r],c,b);
      }
    }
  }
  if(bestN!=inN&&!verify(Best,bestN)){fprintf(stderr,"final verify failed, fallback\n");bestN=inN;memcpy(Best,Inp,(size_t)inN*K);}
  fprintf(stderr,"final N=%d t=%.0f\n",bestN,now()-t0);
  output();
  return 0;
}
