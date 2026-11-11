#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <sys/mman.h>
#include <dirent.h>
#include <sys/stat.h>
// DOMAIN: k=30..50
#define TEMP 0.35
#define DRIFT 30
#define T 4
#define K 38
#define V 6
_Static_assert(T==4,"strength 4 code");
#define W2 (V*V)
#define W3 (V*V*V)
#define NT (V*V*V*V)
#define NS (K*(K-1)*(K-2)*(K-3)/24)
#define PC ((K-1)*(K-2)*(K-3)/6)
#define NG (NS*NT)
#define MAXN 9000
_Static_assert(MAXN<65536,"row index in 16 bits");
#define HMAX (1<<21)
#define TL 3420.0
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
static unsigned char Tmp[MAXN][K];static int nfiles=0;static double scan_end;
static int try_file(const char*p){
  struct stat st;if(stat(p,&st)||!S_ISREG(st.st_mode))return 0;
  if(st.st_size<(off_t)NT*K*2||st.st_size>(off_t)16000000)return 0;
  if(++nfiles>6000)return 0;
  FILE*f=fopen(p,"r");if(!f)return 0;
  int c=0,cols=-1,n=0,v=-1,ch,bad=0;
  while((ch=fgetc(f))!=EOF){
    if(ch>='0'&&ch<='9'){v=(v<0?0:v*10)+(ch-'0');if(v>=V){bad=1;break;}}
    else if(ch==' '||ch=='\t'||ch=='\n'||ch=='\r'||ch==','){
      if(v>=0){if(c<K)Tmp[n][c]=(unsigned char)v;c++;v=-1;}
      if(ch=='\n'&&c>0){
        if(cols<0){cols=c;if(cols<K){bad=1;break;}}else if(c!=cols){bad=1;break;}
        n++;c=0;if(n>=MAXN){bad=1;break;}}
    }else{bad=1;break;}
  }
  fclose(f);
  if(bad)return 0;
  if(v>=0){if(c<K)Tmp[n][c]=(unsigned char)v;c++;}
  if(c>0){if(cols<0)cols=c;if(c!=cols)return 0;n++;}
  if(cols<K||n<NT||n>=N)return 0;
  if(!verify(Tmp,n))return 0;
  return n;
}
static int dirok(const char*nm,int depth){
  const char*q=strstr(nm,"ca-");
  if(q){int t,k,v;if(sscanf(q,"ca-%d-%d-%d",&t,&k,&v)==3)return t>=T&&k>=K&&v==V;}
  return depth>=2;
}
static void scan(const char*d,int depth){
  if(now()>scan_end)return;
  DIR*D=opendir(d);if(!D)return;struct dirent*e;
  while((e=readdir(D))&&now()<=scan_end){
    if(e->d_name[0]=='.')continue;
    char p[4096];snprintf(p,sizeof p,"%s/%s",d,e->d_name);
    struct stat st;if(stat(p,&st))continue;
    if(S_ISDIR(st.st_mode)){if(depth>0&&dirok(e->d_name,depth))scan(p,depth-1);}
    else if(S_ISREG(st.st_mode)){int n=try_file(p);
      if(n>0&&n<N){memcpy(A,Tmp,(size_t)n*K);N=n;fprintf(stderr,"projection seed %s -> N=%d\n",p,n);}}
  }
  closedir(D);
}
static void wide_seed(void){
  scan_end=now()+40;
  static char r1[4096],r2[4200],r3[4300],r4[4400];
  const char*bc=getenv("BEST_CANDIDATE");
  if(bc){snprintf(r1,sizeof r1,"%s",bc);char*sl=strrchr(r1,'/');if(sl){if(sl==r1)sl[1]=0;else *sl=0;}else strcpy(r1,".");}
  else strcpy(r1,".");
  snprintf(r2,sizeof r2,"%s/..",r1);snprintf(r3,sizeof r3,"%s/../..",r1);snprintf(r4,sizeof r4,"%s/../../..",r1);
  scan(r1,1);scan(r2,2);scan(r3,2);scan(".",1);scan("..",2);scan("../..",2);scan(r4,3);
  fprintf(stderr,"wide scan: %d files tried, N=%d\n",nfiles,N);
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
static void reweight(void){
  memset(U,0,sizeof U);memset(P,0,sizeof P);
  for(size_t g=0;g<(size_t)NG;g++)tw[g]=(unsigned char)(1+(tw[g]-1)/2);
  for(int s=0;s<NS;s++){int*cl=scol[s];
    for(int t=0;t<NT;t++){int g=s*NT+t;uint32_t v=cw[g],n=v&0xFFFFu;
      if(n==1)upd((int)(v>>16),cl[0],cl[1],cl[2],cl[3],tw[g]);}}
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
#define SMAX 8192
static unsigned char Sc[SMAX][4],Sv[SMAX][4],Sm[SMAX];static int Sw[SMAX];static uint64_t Smk[SMAX];
static long long rsv_try=0,rsv_acc=0;
static int resolve_row(int r,int hi,long long step){
  static int gl[NS];const unsigned char*xr=A[r];int n=0;
  if(holes>2000)return 0;
  rsv_try++;
  for(int s=0;s<NS;s++)gl[s]=tidx(xr,s);
  for(int s=0;s<NS;s++){if(s+16<NS)__builtin_prefetch(&cw[gl[s+16]]);
    int g=gl[s];if((cw[g]&0xFFFFu)!=1)continue;if(n>=SMAX-holes)break;
    int*cl=scol[s];uint64_t m=0;for(int k=0;k<4;k++){Sc[n][k]=(unsigned char)cl[k];Sv[n][k]=xr[cl[k]];m|=1ULL<<cl[k];}
    Smk[n]=m;Sw[n]=tw[g];n++;}
  int base=0;for(int i=0;i<n;i++)base+=Sw[i];
  for(int i=0;i<holes&&n<SMAX;i++){memcpy(Sc[n],hc[i],4);memcpy(Sv[n],hv[i],4);Smk[n]=hm[i];Sw[n]=hw[i];n++;}
  unsigned char x[K];memcpy(x,xr,K);
  for(int k=0;k<4;k++)x[hc[hi][k]]=hv[hi][k];
  int score=0;
  for(int i=0;i<n;i++){int m=0;for(int k=0;k<4;k++)m+=(x[Sc[i][k]]!=Sv[i][k]);Sm[i]=(unsigned char)m;if(!m)score+=Sw[i];}
  for(int sw=0;sw<10;sw++){int imp=0;int off=(int)(rnd()%K);
    for(int cc=0;cc<K;cc++){int c=(cc+off)%K;uint64_t bit=1ULL<<c;int gain[V];for(int b=0;b<V;b++)gain[b]=0;
      for(int i=0;i<n;i++){if(!(Smk[i]&bit))continue;int p=0;while(Sc[i][p]!=c)p++;
        int m=Sm[i]-(x[c]!=Sv[i][p]);if(m)continue;gain[Sv[i][p]]+=Sw[i];}
      int ob=x[c],bb=ob;for(int b=0;b<V;b++)if(gain[b]>gain[bb])bb=b;
      if(bb!=ob){
        for(int i=0;i<n;i++){if(!(Smk[i]&bit))continue;int p=0;while(Sc[i][p]!=c)p++;
          Sm[i]=(unsigned char)(Sm[i]+(bb!=Sv[i][p])-(ob!=Sv[i][p]));}
        score+=gain[bb]-gain[ob];x[c]=(unsigned char)bb;imp=1;}
    }
    if(!imp)break;}
  if(!(score>base||(score==base&&memcmp(x,xr,K)&&(rnd()&3)==0)))return 0;
  rsv_acc++;
  for(int c=0;c<K;c++)if(x[c]!=A[r][c]){apply(r,c,x[c]);tabu[r][c]=step+10+(long long)(rnd()%20);}
  return 1;
}
static unsigned char Src[MAXN][K],Out[MAXN][K];static int srcN=0,srcCols=0,outN=0;static long srckey=1L<<40;
#define CPEN 500
#define EXTLS 520.0
#define WFR 0.55
#define WCAPX 250
#define NEXT 2
#define MARGIN 600
#define SAT0 1.2
#define SAT1 0.08
static int verify_sub(unsigned char(*X)[K],int n,int kk){
  static unsigned char seen[NT];
  for(int s=0;s<NS;s++){if(scol[s][3]>=kk)continue;memset(seen,0,NT);int d=0;
    for(int r=0;r<n;r++){int id=tidx(X[r],s)-s*NT;if(!seen[id]){seen[id]=1;d++;}}
    if(d!=NT)return 0;}
  return 1;
}
static int try_sub(const char*p){
  struct stat st;if(stat(p,&st)||!S_ISREG(st.st_mode))return 0;
  if(st.st_size<(off_t)NT*(K-3)*2||st.st_size>(off_t)16000000)return 0;
  if(++nfiles>6000)return 0;
  FILE*f=fopen(p,"r");if(!f)return 0;
  int c=0,cols=-1,n=0,v=-1,ch,bad=0;
  while((ch=fgetc(f))!=EOF){
    if(ch>='0'&&ch<='9'){v=(v<0?0:v*10)+(ch-'0');if(v>=V){bad=1;break;}}
    else if(ch==' '||ch=='\t'||ch=='\n'||ch=='\r'||ch==','){
      if(v>=0){if(c<K)Tmp[n][c]=(unsigned char)v;c++;v=-1;}
      if(ch=='\n'&&c>0){
        if(cols<0){cols=c;if(cols<K-3||cols>=K){bad=1;break;}}else if(c!=cols){bad=1;break;}
        n++;c=0;if(n>=MAXN){bad=1;break;}}
    }else{bad=1;break;}
  }
  fclose(f);
  if(bad)return 0;
  if(v>=0){if(c<K)Tmp[n][c]=(unsigned char)v;c++;}
  if(c>0){if(cols<0)cols=c;if(c!=cols)return 0;n++;}
  if(cols<K-3||cols>=K||n<NT)return 0;
  long key=n+(long)CPEN*(K-cols);if(key>=srckey)return 0;
  if(!verify_sub(Tmp,n,cols))return 0;
  srckey=key;srcN=n;srcCols=cols;memcpy(Src,Tmp,(size_t)n*K);
  fprintf(stderr,"extension source %s: %d rows x %d cols\n",p,n,cols);
  return n;
}
static int dirok2(const char*nm,int depth){
  const char*q=strstr(nm,"ca-");
  if(q){int t,k,v;if(sscanf(q,"ca-%d-%d-%d",&t,&k,&v)==3)return t>=T&&k>=K-3&&k<K&&v==V;}
  return depth>=2;
}
static void scan2(const char*d,int depth){
  if(now()>scan_end)return;
  DIR*D=opendir(d);if(!D)return;struct dirent*e;
  while((e=readdir(D))&&now()<=scan_end){
    if(e->d_name[0]=='.')continue;
    char p[4096];snprintf(p,sizeof p,"%s/%s",d,e->d_name);
    struct stat st;if(stat(p,&st))continue;
    if(S_ISDIR(st.st_mode)){if(depth>0&&dirok2(e->d_name,depth))scan2(p,depth-1);}
    else if(S_ISREG(st.st_mode))try_sub(p);
  }
  closedir(D);
}
static void sub_seed(void){
  scan_end=now()+60;nfiles=0;
  static char r1[4096],r2[4200],r3[4300],r4[4400];
  const char*bc=getenv("BEST_CANDIDATE");
  if(bc){snprintf(r1,sizeof r1,"%s",bc);char*sl=strrchr(r1,'/');if(sl){if(sl==r1)sl[1]=0;else *sl=0;}else strcpy(r1,".");}
  else strcpy(r1,".");
  snprintf(r2,sizeof r2,"%s/..",r1);snprintf(r3,sizeof r3,"%s/../..",r1);snprintf(r4,sizeof r4,"%s/../../..",r1);
  scan2("..",2);scan2(r2,2);scan2(r1,1);scan2(r3,2);scan2("../..",2);scan2(".",1);scan2(r4,3);
  fprintf(stderr,"sub scan: %d files tried, source %d rows x %d cols\n",nfiles,srcN,srcCols);
}
static void col_ls(double tend){
  for(int sw=0;sw<400&&now()<tend&&!stopf;sw++){long imp=0;int moves=0;
    for(int r=0;r<N&&!stopf;r++){
      for(int c=srcCols;c<K;c++){unsigned char*x=A[r];int a=x[c];E*ec=ce[c];int gain[V],loss=0;
        for(int b=0;b<V;b++)gain[b]=0;
        for(int j=0;j<PF&&j<PC;j++){E*f=&ec[j];__builtin_prefetch(&cw[PFB(f,x)]);}
        for(int j=0;j<PC;j++){E*e=&ec[j];
          if(j+PF<PC){E*f=&ec[j+PF];int bb=PFB(f,x);__builtin_prefetch(&cw[bb]);if(f->wc>1)__builtin_prefetch(&cw[bb+(V-1)*f->wc]);}
          int base=PFB(e,x);
          for(int b=0;b<V;b++){uint32_t n=cw[base+b*e->wc]&0xFFFFu;if(b==a){if(n==1)loss++;}else if(n==0)gain[b]++;}}
        int bb=-1;for(int b=0;b<V;b++)if(b!=a&&(bb<0||gain[b]>gain[bb]))bb=b;
        if(gain[bb]>loss||(gain[bb]==loss&&(rnd()&7)==0)){imp+=gain[bb]-loss;apply(r,c,bb);moves++;}
      }}
    fprintf(stderr,"col_ls sweep %d: holes=%d imp=%ld moves=%d t=%.0f\n",sw,holes,imp,moves,now());
    if(imp<=holes/2000+1)break;}
}
static void col_walk(int c,double tend){
  double ts=now();if(tend<=ts+1)return;double span=tend-ts;
  static unsigned char cb[MAXN];int bh=holes,h0=holes;for(int r=0;r<N;r++)cb[r]=A[r][c];
  long sweep=0;long mv=0;
  while(now()<tend&&!stopf&&holes>0){
    double fr=(now()-ts)/span;if(fr>1)fr=1;double temp=0.5*pow(0.1,fr);
    int off=(int)(rnd()%N);
    for(int q=0;q<N&&!stopf&&holes>0;q++){int r=q+off;if(r>=N)r-=N;
      if((q&255)==0&&now()>tend)break;
      unsigned char*x=A[r];int a=x[c];E*ec=ce[c];int gain[V],loss=0;
      for(int b=0;b<V;b++)gain[b]=0;
      for(int j=0;j<PF&&j<PC;j++){E*f=&ec[j];__builtin_prefetch(&cw[PFB(f,x)]);}
      for(int j=0;j<PC;j++){E*e=&ec[j];
        if(j+PF<PC){E*f=&ec[j+PF];int bb=PFB(f,x);__builtin_prefetch(&cw[bb]);if(f->wc>1)__builtin_prefetch(&cw[bb+(V-1)*f->wc]);}
        int base=PFB(e,x);
        for(int b=0;b<V;b++){uint32_t n=cw[base+b*e->wc]&0xFFFFu;if(b==a){if(n==1)loss++;}else if(n==0)gain[b]++;}}
      int bb=-1,ties=0;for(int b=0;b<V;b++){if(b==a)continue;
        if(bb<0||gain[b]>gain[bb]){bb=b;ties=1;}else if(gain[b]==gain[bb]&&(rnd()%(++ties))==0)bb=b;}
      int d=gain[bb]-loss;int acc=d>0||(d==0&&(rnd()&1));
      if(!acc&&d<0){double u=(double)(rnd()>>11)*(1.0/9007199254740992.0);acc=u<exp(d/temp);}
      if(acc){apply(r,c,bb);mv++;
        if(holes<bh){bh=holes;for(int z=0;z<N;z++)cb[z]=A[z][c];}}
    }
    sweep++;
  }
  if(holes>bh){for(int r=0;r<N;r++)if(A[r][c]!=cb[r])apply(r,c,cb[r]);}
  fprintf(stderr,"col_walk col %d: %ld sweeps %ld moves holes %d -> %d\n",c,sweep,mv,h0,holes);
}
static unsigned char Cb[MAXN][K];
static void col_sa(double tend){
  double ts=now();if(tend<=ts+1)return;
  int bh=holes;for(int r=0;r<N;r++)memcpy(Cb[r],A[r],K);
  long sweep=0;
  while(now()<tend&&!stopf){
    double fr=(now()-ts)/(tend-ts);if(fr>1)fr=1;
    double temp=SAT0*pow(SAT1/SAT0,fr);
    int moves=0;
    for(int r=0;r<N&&!stopf;r++){
      if((r&255)==0&&now()>tend)break;
      for(int c=srcCols;c<K;c++){unsigned char*x=A[r];int a=x[c];E*ec=ce[c];int gain[V],loss=0;
        for(int b=0;b<V;b++)gain[b]=0;
        for(int j=0;j<PF&&j<PC;j++){E*f=&ec[j];__builtin_prefetch(&cw[PFB(f,x)]);}
        for(int j=0;j<PC;j++){E*e=&ec[j];
          if(j+PF<PC){E*f=&ec[j+PF];int bb=PFB(f,x);__builtin_prefetch(&cw[bb]);if(f->wc>1)__builtin_prefetch(&cw[bb+(V-1)*f->wc]);}
          int base=PFB(e,x);
          for(int b=0;b<V;b++){uint32_t n=cw[base+b*e->wc]&0xFFFFu;if(b==a){if(n==1)loss++;}else if(n==0)gain[b]++;}}
        int bb=-1,ties=0;for(int b=0;b<V;b++){if(b==a)continue;
          if(bb<0||gain[b]>gain[bb]){bb=b;ties=1;}else if(gain[b]==gain[bb]&&(rnd()%(++ties))==0)bb=b;}
        int d=gain[bb]-loss;int acc=d>0||(d==0&&(rnd()&1));
        if(!acc&&d<0){double u=(double)(rnd()>>11)*(1.0/9007199254740992.0);acc=u<exp(d/temp);}
        if(acc){apply(r,c,bb);moves++;}
      }}
    if(holes<bh){bh=holes;for(int r=0;r<N;r++)memcpy(Cb[r]+srcCols,A[r]+srcCols,K-srcCols);}
    if((sweep&7)==0)fprintf(stderr,"col_sa sweep %ld: holes=%d best=%d temp=%.3f moves=%d t=%.0f\n",sweep,holes,bh,temp,moves,now());
    sweep++;
  }
  if(holes>bh){for(int r=0;r<N;r++)for(int c=srcCols;c<K;c++)if(A[r][c]!=Cb[r][c])apply(r,c,Cb[r][c]);}
  fprintf(stderr,"col_sa done: %ld sweeps, holes=%d\n",sweep,holes);
}
static void col_wls(double tend){
  double ts=now();if(tend<=ts+1||holes==0)return;
  int bh=holes,h0=holes;for(int r=0;r<N;r++)memcpy(Cb[r]+srcCols,A[r]+srcCols,K-srcCols);
  long sweep=0;
  while(now()<tend&&!stopf&&holes>0){
    int moves=0;long imp=0;
    for(int r=0;r<N&&!stopf&&holes>0;r++){
      if((r&255)==0&&now()>tend)break;
      for(int c=srcCols;c<K;c++){unsigned char*x=A[r];int a=x[c];E*ec=ce[c];int gain[V],loss=0;
        for(int b=0;b<V;b++)gain[b]=0;
        for(int j=0;j<PF&&j<PC;j++){E*f=&ec[j];__builtin_prefetch(&cw[PFB(f,x)]);}
        for(int j=0;j<PC;j++){E*e=&ec[j];
          if(j+PF<PC){E*f=&ec[j+PF];int bb=PFB(f,x);__builtin_prefetch(&cw[bb]);if(f->wc>1)__builtin_prefetch(&cw[bb+(V-1)*f->wc]);}
          int base=PFB(e,x);
          for(int b=0;b<V;b++){int g=base+b*e->wc;uint32_t n=cw[g]&0xFFFFu;
            if(b==a){if(n==1)loss+=tw[g];}else if(n==0)gain[b]+=tw[g];}}
        int bb=-1,ties=0;for(int b=0;b<V;b++){if(b==a)continue;
          if(bb<0||gain[b]>gain[bb]){bb=b;ties=1;}else if(gain[b]==gain[bb]&&(rnd()%(++ties))==0)bb=b;}
        int d=gain[bb]-loss;
        if(d>0||(d==0&&(rnd()&3)==0)){apply(r,c,bb);moves++;imp+=d;}
      }}
    if(holes<bh){bh=holes;for(int r=0;r<N;r++)memcpy(Cb[r]+srcCols,A[r]+srcCols,K-srcCols);}
    for(int i=0;i<holes;i++){int g=hl[i];if(tw[g]<WCAPX)tw[g]++;}
    if((sweep&3)==0)fprintf(stderr,"col_wls sweep %ld: holes=%d best=%d moves=%d wimp=%ld t=%.0f\n",sweep,holes,bh,moves,imp,now());
    sweep++;
  }
  if(holes>bh){for(int r=0;r<N;r++)for(int c=srcCols;c<K;c++)if(A[r][c]!=Cb[r][c])apply(r,c,Cb[r][c]);}
  fprintf(stderr,"col_wls done: %ld sweeps, holes %d -> %d\n",sweep,h0,holes);
}
static void build_rows(void){
  static unsigned char bx[K];
  while(holes>0&&N<MAXN-2&&!stopf){
    int bestc=-1;
    for(int tr=0;tr<6;tr++){
      int req[K];for(int c=0;c<K;c++)req[c]=-1;
      int off=(int)(rnd()%holes);
      for(int q=0;q<holes;q++){int i=(q+off)%holes;int ok=1;
        for(int k=0;k<4;k++){int c=hc[i][k];if(req[c]>=0&&req[c]!=hv[i][k]){ok=0;break;}}
        if(!ok)continue;for(int k=0;k<4;k++)req[hc[i][k]]=hv[i][k];}
      unsigned char x[K];for(int c=0;c<K;c++)x[c]=(unsigned char)(req[c]>=0?req[c]:(int)(rnd()%V));
      for(int sw=0;sw<10;sw++){int ch=0;
        for(int c=0;c<K;c++){uint64_t bit=1ULL<<c;int g[V];for(int b=0;b<V;b++)g[b]=0;
          for(int i=0;i<holes;i++){if(!(hm[i]&bit))continue;int ok=1,p=0;
            for(int k=0;k<4;k++){int q=hc[i][k];if(q==c){p=k;continue;}if(x[q]!=hv[i][k]){ok=0;break;}}
            if(ok)g[hv[i][p]]++;}
          int bb=x[c];for(int b=0;b<V;b++)if(g[b]>g[bb])bb=b;
          if(bb!=x[c]){x[c]=(unsigned char)bb;ch=1;}}
        if(!ch)break;}
      int cov=0;for(int i=0;i<holes;i++){int ok=1;for(int k=0;k<4;k++)if(x[hc[i][k]]!=hv[i][k]){ok=0;break;}cov+=ok;}
      if(cov>bestc){bestc=cov;memcpy(bx,x,K);}
    }
    add_row(bx);
    if((N&15)==0)fprintf(stderr,"build_rows: N=%d holes=%d last cover %d\n",N,holes,bestc);
  }
}
static void output(void){
  for(int r=0;r<outN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",Out[r][c]);printf("\n");}
  fflush(stdout);
}
static void passthru(void){
  const char*fn=getenv("BEST_CANDIDATE");if(fn){FILE*f=fopen(fn,"r");if(f){int ch;while((ch=fgetc(f))!=EOF)putchar(ch);fclose(f);}}
}
int main(void){
  double t0=now(),hard=t0+TL;
  signal(SIGTERM,onterm);signal(SIGINT,onterm);
  {size_t al=(size_t)1<<21;
   size_t sc=((size_t)NG*sizeof(uint32_t)+al-1)/al*al,sp=((size_t)NG*sizeof(int)+al-1)/al*al;
   void*p1=0,*p2=0;
   if(posix_memalign(&p1,al,sc)||posix_memalign(&p2,al,sp)){fprintf(stderr,"oom\n");passthru();return 1;}
#ifdef MADV_HUGEPAGE
   madvise(p1,sc,MADV_HUGEPAGE);madvise(p2,sp,MADV_HUGEPAGE);
#endif
   cw=(uint32_t*)p1;pos=(int*)p2;
   tw=(unsigned char*)malloc((size_t)NG);
   if(!tw){fprintf(stderr,"oom tw\n");passthru();return 1;}
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
  if(N==0||!verify(A,N)){fprintf(stderr,"no valid input (%d rows)\n",N);outN=N;memcpy(Out,A,(size_t)N*K);output();return 1;}
  fprintf(stderr,"input valid N=%d\n",N);
  inN=N;memcpy(Inp,A,(size_t)N*K);
  wide_seed();
  outN=N;memcpy(Out,A,(size_t)N*K);
  sub_seed();
  {int useext=0;
   static unsigned char Eb[MAXN][K];int ebN=0;
   for(int ne=0;ne<NEXT&&srcN>0&&srcN<N&&!stopf;ne++){
     double te=now()+EXTLS;
     N=srcN;for(int r=0;r<N;r++){memcpy(A[r],Src[r],K);for(int c=srcCols;c<K;c++)A[r][c]=(unsigned char)(rnd()%V);}
     rebuild();fprintf(stderr,"extension %d start: %d rows, %d new cols, holes=%d\n",ne,N,K-srcCols,holes);
     if(holes<HMAX-1){col_ls(te);
       double tm=now()+(te-now())*WFR;col_sa(tm);
       fprintf(stderr,"extension %d after col_sa: holes=%d\n",ne,holes);
       col_wls(te);build_rows();
       fprintf(stderr,"extension %d done: N=%d holes=%d (input %d)\n",ne,N,holes,outN);
       if(holes==0&&!stopf&&(ebN==0||N<ebN)){ebN=N;memcpy(Eb,A,(size_t)N*K);}}
     memset(tw,1,(size_t)NG);
     N=outN;
   }
   if(ebN>0&&ebN<=outN+MARGIN&&!stopf){N=ebN;memcpy(A,Eb,(size_t)N*K);useext=1;}
   if(!useext){N=outN;memcpy(A,Out,(size_t)N*K);}
   else fprintf(stderr,"using extended array N=%d\n",N);}
  memset(tabu,0,sizeof tabu);
  bestN=N;memcpy(Best,A,(size_t)N*K);fprintf(stderr,"start N=%d\n",N);
  rebuild();
  int rank=0;double LB=300,limp=now();long long step=0;int lbest=1<<30;int kicks=0;int rper=63;
  int lvbest=1<<30;double lvimp=now();int lvN=-1;static unsigned char Lv[MAXN][K];
  double ema=75,tlastb=now(),LB0=300;
  dodelete(pick_weak(0));
  static int cr[NC],ccol[NC],cv[NC],cc2[NC],cv2[NC];
#define NNR 64
  static int nr[NNR],n1[NNR],nv1[NNR],n2[NNR],nv2[NNR];
  while(!stopf){
    if(holes==0){
      bestN=N;memcpy(Best,A,(size_t)N*K);fprintf(stderr,"NEW BEST N=%d t=%.0f step=%lld\n",N,now()-t0,step);
      if(N<outN){outN=N;memcpy(Out,A,(size_t)N*K);fprintf(stderr,"PROGRESS %d\n",-N);}
      {double tn=now(),dt=tn-tlastb;tlastb=tn;ema=0.8*ema+0.2*dt;
       LB0=4*ema;if(LB0<80)LB0=80;if(LB0>300)LB0=300;}
      rank=0;LB=LB0;limp=now();lbest=1<<30;kicks=0;lvbest=1<<30;lvimp=now();lvN=-1;
      if(now()>hard)break;
      if(holes==0)reweight();else rebuild();
      dodelete(pick_weak(0));continue;
    }
    step++;
    if((step&255)==0){
      double tn=now();if(tn>hard)break;
      rper=(tn-limp>20)?15:63;
      if((step&((1<<20)-1))==0)fprintf(stderr,"N=%d holes=%d lbest=%d step=%lld t=%.0f resolve %lld/%lld\n",N,holes,lbest,step,tn-t0,rsv_acc,rsv_try);
      if(tn-limp>KLB&&kicks<KMAX&&(holes<=12||lvN==N)){
        if(kicks>0&&lvN==N&&holes>lvbest){N=lvN;memcpy(A,Lv,(size_t)N*K);memset(tabu,0,sizeof tabu);rebuild();
          fprintf(stderr,"kick failed, restore level best %d holes=%d\n",lvbest,holes);}
        kicks++;if(kicks&1)kick(step);
        else{int hi=(int)(rnd()%holes);int c=hc[hi][rnd()&3];
          if((rnd()&1)&&(hm[hi]>>(K-1)&1ULL))c=K-1;
          double te=now()+40.0;if(te>hard)te=hard;col_walk(c,te);
          memset(tabu,0,sizeof tabu);}
        limp=now();lbest=holes;memcpy(Ck,A,(size_t)N*K);continue;
      }
      if(tn-lvimp>LB){
        kicks=0;
        N=bestN;memcpy(A,Best,(size_t)N*K);memset(tabu,0,sizeof tabu);rebuild();
        rank=(rank+1)%8;if(rank==0)LB*=1.5;
        limp=now();lbest=1<<30;lvbest=1<<30;lvimp=now();lvN=-1;dodelete(pick_weak(rank));continue;
      }
    }
    if(holes<lbest){lbest=holes;limp=now();memcpy(Ck,A,(size_t)N*K);
      if(holes<lvbest||lvN!=N){lvbest=holes;lvimp=limp;lvN=N;memcpy(Lv,A,(size_t)N*K);}}
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
    if((step&rper)==0&&(nc>0||nn>0)){
      int r;if(nn>0&&(nc==0||(rnd()&1)))r=nr[rnd()%nn];else r=cr[rnd()%nc];
      if(resolve_row(r,hi,step))continue;
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
  if(outN!=inN&&!verify(Out,outN)){fprintf(stderr,"final verify failed, fallback\n");outN=inN;memcpy(Out,Inp,(size_t)inN*K);}
  fprintf(stderr,"final N=%d t=%.0f steps=%lld\n",outN,now()-t0,step);
  output();
  return 0;
}
