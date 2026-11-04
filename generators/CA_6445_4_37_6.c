#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <sys/mman.h>
#define TEMP 0.35
#define DRIFT 30
#define T 4
#define K 37
#define V 6
_Static_assert(T==4,"strength 4 code");
#define W2 (V*V)
#define W3 (V*V*V)
#define NT (V*V*V*V)
#define NS (K*(K-1)*(K-2)*(K-3)/24)
#define PC ((K-1)*(K-2)*(K-3)/6)
#define NG (NS*NT)
#define MAXN 7300
_Static_assert(MAXN<65536,"row index in 16 bits");
#define HMAX (1<<21)
#define TL 3430.0
#define NC 1024
typedef struct{int gb;unsigned char o0,o1,o2,pad;unsigned short w0,w1,w2,wc;}E;
static volatile sig_atomic_t stopf=0;
static void onterm(int s){(void)s;stopf=1;}
static uint32_t *cw;static int *pos;static unsigned char *tw;
static int hl[HMAX],holes;static unsigned short hw[HMAX];
_Static_assert(K<=64,"column mask in 64 bits");
static unsigned char hc[HMAX][4],hv[HMAX][4];static uint64_t hm[HMAX];
#define WCAP 48
static int scol[NS][4];static E ce[K][PC];
static unsigned char A[MAXN][K],Best[MAXN][K],Inp[MAXN][K];static int N,bestN,inN;
static long long tabu[MAXN][K];
static int U[MAXN][K];static unsigned short P[MAXN][K][K];
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int tidx(const unsigned char*x,int s){int*c=scol[s];return s*NT+x[c[0]]*W3+x[c[1]]*W2+x[c[2]]*V+x[c[3]];}
static inline void addh(int g){if(holes<HMAX){int i=holes++;pos[g]=i;hw[i]=tw[g];hl[i]=g;
  int s=g/NT,tv=g%NT;int*cl=scol[s];
  hc[i][0]=(unsigned char)cl[0];hc[i][1]=(unsigned char)cl[1];hc[i][2]=(unsigned char)cl[2];hc[i][3]=(unsigned char)cl[3];
  hv[i][0]=(unsigned char)(tv/W3);hv[i][1]=(unsigned char)((tv/W2)%V);hv[i][2]=(unsigned char)((tv/V)%V);hv[i][3]=(unsigned char)(tv%V);
  hm[i]=(1ULL<<cl[0])|(1ULL<<cl[1])|(1ULL<<cl[2])|(1ULL<<cl[3]);}}
static inline void delh(int g){int p=pos[g];if(p<0)return;int l=hl[--holes];hl[p]=l;hw[p]=hw[holes];
  memcpy(hc[p],hc[holes],4);memcpy(hv[p],hv[holes],4);hm[p]=hm[holes];pos[l]=p;pos[g]=-1;}
static inline void upd(int row,int a,int b,int c,int d,int dl){
  int*u=U[row];u[a]+=dl;u[b]+=dl;u[c]+=dl;u[d]+=dl;
  unsigned short(*p)[K]=P[row];
  p[a][b]+=dl;p[b][a]+=dl;p[a][c]+=dl;p[c][a]+=dl;p[a][d]+=dl;p[d][a]+=dl;
  p[b][c]+=dl;p[c][b]+=dl;p[b][d]+=dl;p[d][b]+=dl;p[c][d]+=dl;p[d][c]+=dl;
}
static int verify(unsigned char(*X)[K],int n){
  static unsigned char seen[NT];
  for(int s=0;s<NS;s++){memset(seen,0,NT);int d=0;
    for(int r=0;r<n;r++){int id=tidx(X[r],s)-s*NT;if(!seen[id]){seen[id]=1;d++;}}
    if(d!=NT)return 0;}
  return 1;
}
static void rebuild(void){
  memset(cw,0,(size_t)NG*sizeof(*cw));
  for(int r=0;r<N;r++){uint32_t rr=(uint32_t)r<<16;for(int s=0;s<NS;s++){int g=tidx(A[r],s);cw[g]=(cw[g]+1)^rr;}}
  memset(U,0,sizeof U);memset(P,0,sizeof P);
  for(size_t g=0;g<(size_t)NG;g++)tw[g]=(unsigned char)(1+(tw[g]-1)/2);
  holes=0;
  for(int s=0;s<NS;s++){int*cl=scol[s];
    for(int t=0;t<NT;t++){int g=s*NT+t;pos[g]=-1;uint32_t v=cw[g],n=v&0xFFFFu;
      if(n==0)addh(g);else if(n==1)upd((int)(v>>16),cl[0],cl[1],cl[2],cl[3],tw[g]);}}
}
#define PF 16
#define PFB(f,x) ((f)->gb+(x)[(f)->o0]*(f)->w0+(x)[(f)->o1]*(f)->w1+(x)[(f)->o2]*(f)->w2)
static inline int gainh(const unsigned char*x,int c1,int b1,int c2,int b2){
  int gain=0;uint64_t m=1ULL<<c1;if(c2>=0)m|=1ULL<<c2;
  for(int i=0;i<holes;i++){if(!(hm[i]&m))continue;
    const unsigned char*cl=hc[i],*vv=hv[i];int ok=1;
    for(int k=0;k<4;k++){int q=cl[k];int xv;
      if(q==c1)xv=b1;else if(q==c2)xv=b2;else xv=x[q];
      if(xv!=vv[k]){ok=0;break;}}
    if(ok)gain+=hw[i];}
  return gain;
}
static inline int delta1(int r,int c,int b){return U[r][c]-gainh(A[r],c,b,-1,0);}
static inline int delta2(int r,int c1,int b1,int c2,int b2){
  return U[r][c1]+U[r][c2]-(int)P[r][c1][c2]-gainh(A[r],c1,b1,c2,b2);
}
static unsigned char Ck[MAXN][K];
static inline void apply(int r,int c,int b){
  unsigned char*x=A[r];int a=x[c];E*ec=ce[c];uint32_t rr=(uint32_t)r<<16;
  for(int j=0;j<PF&&j<PC;j++){E*f=&ec[j];int bb=PFB(f,x);__builtin_prefetch(&cw[bb+a*f->wc],1);__builtin_prefetch(&cw[bb+b*f->wc],1);}
  for(int j=0;j<PC;j++){E*e=&ec[j];
    if(j+PF<PC){E*f=&ec[j+PF];int bb=PFB(f,x);__builtin_prefetch(&cw[bb+a*f->wc],1);__builtin_prefetch(&cw[bb+b*f->wc],1);}
    int base=e->gb+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;
    int ga=base+a*e->wc,gb=base+b*e->wc;
    uint32_t v=cw[ga],n=v&0xFFFFu;v=(v-1)^rr;cw[ga]=v;
    if(n==1){addh(ga);upd(r,c,e->o0,e->o1,e->o2,-(int)tw[ga]);}
    else if(n==2)upd((int)(v>>16),c,e->o0,e->o1,e->o2,tw[ga]);
    v=cw[gb];n=v&0xFFFFu;
    if(n==0){delh(gb);upd(r,c,e->o0,e->o1,e->o2,tw[gb]);}
    else if(n==1)upd((int)(v>>16),c,e->o0,e->o1,e->o2,-(int)tw[gb]);
    cw[gb]=(v+1)^rr;}
  x[c]=(unsigned char)b;
}
static void restore_ck(void){
  for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(A[r][c]!=Ck[r][c])apply(r,c,Ck[r][c]);
}
static void delete_row(int r){
  uint32_t rr=(uint32_t)r<<16;
  for(int s=0;s<NS;s++){int g=tidx(A[r],s);int*cl=scol[s];
    uint32_t v=cw[g],n=v&0xFFFFu;v=(v-1)^rr;cw[g]=v;
    if(n==1){addh(g);upd(r,cl[0],cl[1],cl[2],cl[3],-(int)tw[g]);}
    else if(n==2)upd((int)(v>>16),cl[0],cl[1],cl[2],cl[3],tw[g]);}
  N--;
  if(r!=N){uint32_t m=(uint32_t)(r^N)<<16;
    for(int s=0;s<NS;s++)cw[tidx(A[N],s)]^=m;
    memcpy(A[r],A[N],K);memcpy(tabu[r],tabu[N],sizeof tabu[0]);
    memcpy(U[r],U[N],sizeof U[0]);memcpy(P[r],P[N],sizeof P[0]);}
}
#define NCAND 64
#define NOFIX 8
static int excl=-1;
static int pick_weak(int rank){
  static int u[MAXN],rk[MAXN];static char used[MAXN];
  for(int r=0;r<N;r++){int c=0;for(int q=0;q<K;q++)c+=U[r][q];u[r]=c/4;rk[r]=(int)(rnd()&1023);}
  int M=NCAND<N-1?NCAND:N-1;static int cand[NCAND];static long scr[NCAND];
  memset(used,0,N);if(excl>=0&&excl<N)used[excl]=1;
  for(int k=0;k<M;k++){int b=-1;
    for(int r=0;r<N;r++)if(!used[r]&&(b<0||u[r]*1024+rk[r]<u[b]*1024+rk[b]))b=r;
    used[b]=1;cand[k]=b;}
  for(int k=0;k<M;k++){int w=cand[k];const unsigned char*xw=A[w];long tot=0;
    for(int s=0;s<NS;s++){int g=tidx(xw,s);if((cw[g]&0xFFFFu)!=1)continue;
      int*cl=scol[s];int best=NOFIX;
      for(int r=0;r<N&&best>0;r++){if(r==w)continue;const unsigned char*x=A[r];int mm=0,mc=0;
        for(int i=0;i<4;i++)if(x[cl[i]]!=xw[cl[i]]){mc=cl[i];if(++mm>1)break;}
        if(mm==1&&U[r][mc]<best)best=U[r][mc];}
      tot+=best;}
    scr[k]=tot*4096+(long)u[w]*4+(long)(rnd()&3);}
  for(int i=1;i<M;i++){long s0=scr[i];int c0=cand[i];int j=i-1;
    while(j>=0&&scr[j]>s0){scr[j+1]=scr[j];cand[j+1]=cand[j];j--;}scr[j+1]=s0;cand[j+1]=c0;}
  int pick=cand[rank%M];
  fprintf(stderr,"delete row %d unique=%d fixcost=%ld (rank %d) N->%d\n",pick,u[pick],scr[rank%M]/4096,rank,N-1);
  return pick;
}
/* drain: cover victim's unique tuples by free (U==0) single-entry changes in other rows */
static int drain(int w){
  int moved=0;
  for(int pass=0;pass<4;pass++){int mp=0;
    for(int s=0;s<NS;s++){int g=tidx(A[w],s);if((cw[g]&0xFFFFu)!=1)continue;
      int*cl=scol[s];const unsigned char*xw=A[w];int off=(int)(rnd()%N);
      for(int q=0;q<N;q++){int r=q+off;if(r>=N)r-=N;if(r==w)continue;
        const unsigned char*x=A[r];int mm=0,mc=-1;
        for(int i=0;i<4;i++)if(x[cl[i]]!=xw[cl[i]]){mc=cl[i];if(++mm>1)break;}
        if(mm==1&&U[r][mc]==0){apply(r,mc,xw[mc]);mp++;break;}}}
    moved+=mp;if(mp==0)break;}
  return moved;
}
static void dodelete(int w){
  int u0=0;for(int q=0;q<K;q++)u0+=U[w][q];
  int m=drain(w);int u1=0;for(int q=0;q<K;q++)u1+=U[w][q];
  fprintf(stderr,"drain: %d moves, unique %d -> %d\n",m,u0/4,u1/4);
  delete_row(w);
}
static int kick(long long step){
  int req[K];for(int c=0;c<K;c++)req[c]=-1;
  int off=holes?(int)(rnd()%holes):0,used=0;
  for(int q=0;q<holes;q++){int i=(q+off)%holes;int ok=1;
    for(int k=0;k<4;k++){int c=hc[i][k];if(req[c]>=0&&req[c]!=hv[i][k]){ok=0;break;}}
    if(!ok)continue;for(int k=0;k<4;k++)req[hc[i][k]]=hv[i][k];used++;}
  int br=-1,ties=0;long bc=1L<<60;
  for(int r=0;r<N;r++){long c=0;for(int q=0;q<K;q++)if(req[q]>=0&&A[r][q]!=req[q])c+=U[r][q]+1;
    if(c<bc){bc=c;br=r;ties=1;}else if(c==bc&&(rnd()%(++ties))==0)br=r;}
  if(br<0)return 0;
  for(int q=0;q<K;q++)if(req[q]>=0&&A[br][q]!=req[q]){apply(br,q,req[q]);tabu[br][q]=step+40+(long long)(rnd()%40);}
  fprintf(stderr,"kick: %d holes merged, row %d cost %ld -> holes=%d\n",used,br,bc,holes);
  return used;
}
static void add_row(const unsigned char*x){
  int r=N;memcpy(A[r],x,K);memset(U[r],0,sizeof U[0]);memset(P[r],0,sizeof P[0]);memset(tabu[r],0,sizeof tabu[0]);
  uint32_t rr=(uint32_t)r<<16;
  for(int s=0;s<NS;s++){int g=tidx(x,s);int*cl=scol[s];uint32_t v=cw[g],n=v&0xFFFFu;
    if(n==0){delh(g);upd(r,cl[0],cl[1],cl[2],cl[3],tw[g]);}
    else if(n==1)upd((int)(v>>16),cl[0],cl[1],cl[2],cl[3],-(int)tw[g]);
    cw[g]=(v+1)^rr;}
  N++;
}
static int exchange(void){
  int req[K];for(int c=0;c<K;c++)req[c]=-1;
  int off=holes?(int)(rnd()%holes):0,used=0;
  for(int q=0;q<holes;q++){int i=(q+off)%holes;int ok=1;
    for(int k=0;k<4;k++){int c=hc[i][k];if(req[c]>=0&&req[c]!=hv[i][k]){ok=0;break;}}
    if(!ok)continue;for(int k=0;k<4;k++)req[hc[i][k]]=hv[i][k];used++;}
  unsigned char x[K];for(int c=0;c<K;c++)x[c]=(unsigned char)(req[c]>=0?req[c]:(int)(rnd()%V));
  int h0=holes;add_row(x);int h1=holes;
  excl=N-1;int w=pick_weak(0);excl=-1;
  dodelete(w);
  fprintf(stderr,"exchange: %d holes merged into new row, holes %d -> %d -> %d\n",used,h0,h1,holes);
  return used;
}
#define KMAX 4
#define KLB 150.0
static void output(void){
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",Best[r][c]);printf("\n");}
  fflush(stdout);
}
int main(void){
  double t0=now(),hard=t0+TL;
  signal(SIGTERM,onterm);signal(SIGINT,onterm);
  {size_t al=(size_t)1<<21;
   size_t sc=((size_t)NG*sizeof(uint32_t)+al-1)/al*al,sp=((size_t)NG*sizeof(int)+al-1)/al*al;
   void*p1=0,*p2=0;
   if(posix_memalign(&p1,al,sc)||posix_memalign(&p2,al,sp)){fprintf(stderr,"oom\n");
     const char*fn=getenv("BEST_CANDIDATE");if(fn){FILE*f=fopen(fn,"r");if(f){int ch;while((ch=fgetc(f))!=EOF)putchar(ch);fclose(f);}}
     return 1;}
#ifdef MADV_HUGEPAGE
   madvise(p1,sc,MADV_HUGEPAGE);madvise(p2,sp,MADV_HUGEPAGE);
#endif
   cw=(uint32_t*)p1;pos=(int*)p2;
   tw=(unsigned char*)malloc((size_t)NG);
   if(!tw){fprintf(stderr,"oom tw\n");
     const char*fn2=getenv("BEST_CANDIDATE");if(fn2){FILE*f=fopen(fn2,"r");if(f){int ch;while((ch=fgetc(f))!=EOF)putchar(ch);fclose(f);}}
     return 1;}
   memset(tw,1,(size_t)NG);}
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
  int rank=0;double LB=300,limp=now();long long step=0;int lbest=1<<30;int kicks=0;
  dodelete(pick_weak(0));
  static int cr[NC],ccol[NC],cv[NC],cc2[NC],cv2[NC];
#define NNR 64
  static int nr[NNR],n1[NNR],nv1[NNR],n2[NNR],nv2[NNR];
  while(!stopf){
    if(holes==0){
      bestN=N;memcpy(Best,A,(size_t)N*K);fprintf(stderr,"NEW BEST N=%d t=%.0f step=%lld\nPROGRESS %d\n",N,now()-t0,step,-N);
      rank=0;LB=300;limp=now();lbest=1<<30;kicks=0;
      if(now()>hard)break;
      rebuild();
      dodelete(pick_weak(0));continue;
    }
    step++;
    if((step&255)==0){
      double tn=now();if(tn>hard)break;
      if((step&((1<<20)-1))==0)fprintf(stderr,"N=%d holes=%d lbest=%d step=%lld t=%.0f\n",N,holes,lbest,step,tn-t0);
      if(tn-limp>KLB&&kicks<KMAX&&holes<=12){
        kicks++;if(kicks&1)kick(step);else exchange();
        limp=now();lbest=holes;memcpy(Ck,A,(size_t)N*K);continue;
      }
      if(tn-limp>LB){
        kicks=0;
        N=bestN;memcpy(A,Best,(size_t)N*K);memset(tabu,0,sizeof tabu);rebuild();
        rank=(rank+1)%8;if(rank==0)LB*=1.5;
        limp=now();lbest=1<<30;dodelete(pick_weak(rank));continue;
      }
    }
    if(holes<lbest){lbest=holes;limp=now();memcpy(Ck,A,(size_t)N*K);}
    else if(holes>lbest+DRIFT){restore_ck();memset(tabu,0,sizeof tabu);continue;}
    int hi=(int)(rnd()%holes);int h=hl[hi];int ss=h/NT,tv=h%NT;
    int val[4]={tv/W3,(tv/W2)%V,(tv/V)%V,tv%V};int*col=scol[ss];
    int nc=0,seen=0,nn=0,nseen=0;
    int use2=(holes<=4)||((holes<=16)&&(rnd()&1));int lim=use2?2:1;int brk=(lim==1)?2:lim;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int mm=0,m1=0,m2=-1;
      for(int i=0;i<4;i++)if(x[col[i]]!=val[i]){if(mm==0)m1=i;else if(mm==1)m2=i;mm++;if(mm>brk)break;}
      if(lim==1&&mm==2){nseen++;int sl=nn<NNR?nn++:(int)(rnd()%nseen);
        if(sl<NNR){nr[sl]=r;n1[sl]=col[m1];nv1[sl]=val[m1];n2[sl]=col[m2];nv2[sl]=val[m2];}
        continue;}
      if(mm<1||mm>lim)continue;
      seen++;
      int slot=nc<NC?nc++:(int)(rnd()%seen);
      if(slot<NC){cr[slot]=r;ccol[slot]=col[m1];cv[slot]=val[m1];
        if(mm==2){cc2[slot]=col[m2];cv2[slot]=val[m2];}else cc2[slot]=-1;}
    }
    if(nc==0){
      int r=(int)(rnd()%N),i=(int)(rnd()%4);
      if(A[r][col[i]]!=val[i])apply(r,col[i],val[i]);
      continue;
    }
    int bi=-1,bd=1<<30,ties=0;
    for(int k=0;k<nc;k++){
      if(tabu[cr[k]][ccol[k]]>step)continue;
      int d;
      if(cc2[k]<0)d=delta1(cr[k],ccol[k],cv[k]);
      else{
        if(tabu[cr[k]][cc2[k]]>step)continue;
        d=delta2(cr[k],ccol[k],cv[k],cc2[k],cv2[k]);
      }
      if(d<bd){bd=d;bi=k;ties=1;}else if(d==bd&&(rnd()%(++ties))==0)bi=k;
    }
    if(bd>0&&hw[hi]<WCAP){hw[hi]++;tw[h]=(unsigned char)hw[hi];}
    if(bi<0)continue;
    int acc=(bd<=0)||((double)(rnd()%1000000)<exp(-bd/TEMP)*1e6);
    if(acc){
      apply(cr[bi],ccol[bi],cv[bi]);
      tabu[cr[bi]][ccol[bi]]=step+10+(long long)(rnd()%20);
      if(cc2[bi]>=0){
        apply(cr[bi],cc2[bi],cv2[bi]);
        tabu[cr[bi]][cc2[bi]]=step+10+(long long)(rnd()%20);
      }
    }else{
      int done=0;
      if(nn>0){int off=(int)(rnd()%nn);
        for(int q=0;q<nn&&done<2;q++){int k=(q+off)%nn;int r=nr[k];
          int c=n1[k],b=nv1[k];if(rnd()&1){c=n2[k];b=nv2[k];}
          if(A[r][c]==b||tabu[r][c]>step)continue;
          if(delta1(r,c,b)<=0){apply(r,c,b);tabu[r][c]=step+10+(long long)(rnd()%20);done++;}
        }
      }
      for(int z=0;z<16;z++){
        int r=(int)(rnd()%N),c=(int)(rnd()%K),b=(int)(rnd()%V);
        if(b==A[r][c])continue;
        if(delta1(r,c,b)<=0)apply(r,c,b);
      }
    }
  }
  if(bestN!=inN&&!verify(Best,bestN)){fprintf(stderr,"final verify failed, fallback\n");bestN=inN;memcpy(Best,Inp,(size_t)inN*K);}
  fprintf(stderr,"final N=%d t=%.0f steps=%lld\n",bestN,now()-t0,step);
  output();
  return 0;
}
