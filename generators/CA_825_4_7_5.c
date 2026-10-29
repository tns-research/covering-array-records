#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#define K 7
#define NT 35
#define MAXR 1000
#define TT (NT*625)
static int cols[NT][4],smask[NT];
static unsigned short cnt[TT];static unsigned int W[TT];
static int ulist[TT],upos[TT],nu;
static unsigned char A[MAXR][K];static int nr;
static long long tabu[MAXR][K];
static unsigned char BEST[MAXR][K];static int bestN=0;
static unsigned char X[MAXR*2][K];
static double T0;
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timeval tv;gettimeofday(&tv,0);return tv.tv_sec+tv.tv_usec*1e-6;}
static inline int idx(const unsigned char*r,int t){return t*625+((r[cols[t][0]]*5+r[cols[t][1]])*5+r[cols[t][2]])*5+r[cols[t][3]];}
static inline void inc(int k){if(cnt[k]++==0){int p=upos[k];int l=ulist[--nu];ulist[p]=l;upos[l]=p;upos[k]=-1;}}
static inline void dec(int k){if(--cnt[k]==0){upos[k]=nu;ulist[nu++]=k;}}
static void build(){memset(cnt,0,sizeof cnt);nu=0;
 for(int q=0;q<nr;q++)for(int t=0;t<NT;t++)cnt[idx(A[q],t)]++;
 for(int k=0;k<TT;k++){if(!cnt[k]){upos[k]=nu;ulist[nu++]=k;}else upos[k]=-1;}}
static int remdelta(int q){int d=0;for(int t=0;t<NT;t++)if(cnt[idx(A[q],t)]==1)d++;return d;}
static void delrow(int q){for(int t=0;t<NT;t++)dec(idx(A[q],t));nr--;if(q!=nr)memcpy(A[q],A[nr],K);}
static int verifyArr(unsigned char(*Y)[K],int n){static unsigned char s[TT];memset(s,0,TT);
 for(int r=0;r<n;r++)for(int t=0;t<NT;t++)s[idx(Y[r],t)]=1;for(int i=0;i<TT;i++)if(!s[i])return 0;return 1;}
static int search(double tend){
 long long it=0;int bestnu=nu;
 for(int k=0;k<TT;k++)W[k]=1;
 memset(tabu,0,sizeof(tabu[0])*nr);
 while(nu>0){
  if((it&1023)==0&&now()>tend)break;
  it++;
  if(it%100000==0)for(int k=0;k<TT;k++)W[k]=1+(W[k]-1)*3/4;
  int k=ulist[rnd()%nu];int t=k/625,i=k%625;int v[4]={i/125,(i/25)%5,(i/5)%5,i%5};
  int maxd=(rnd()%20==0)?4:2;
  long long bw=(1LL<<62);int bq=-1,bcm=0,bc=0;
  for(int q=0;q<nr;q++){int d=0,cm=0;
   for(int j=0;j<4;j++){int c=cols[t][j];if(A[q][c]!=v[j]){d++;cm|=1<<c;}}
   if(d>maxd)continue;
   unsigned char w[K];memcpy(w,A[q],K);for(int j=0;j<4;j++)w[cols[t][j]]=v[j];
   long long dw=0;int du=0;
   for(int t2=0;t2<NT;t2++)if(smask[t2]&cm){int a=idx(A[q],t2),b=idx(w,t2);
    if(cnt[a]==1){dw+=W[a];du++;}if(cnt[b]==0){dw-=W[b];du--;}}
   int ist=0;for(int c=0;c<K;c++)if(((cm>>c)&1)&&tabu[q][c]>it){ist=1;break;}
   if(ist&&nu+du>=bestnu)continue;
   if(dw<bw){bw=dw;bq=q;bcm=cm;bc=1;}else if(dw==bw&&rnd()%(++bc)==0){bq=q;bcm=cm;}}
  if(bq<0)continue;
  if(bw>=0)for(int u=0;u<nu;u++)W[ulist[u]]++;
  unsigned char w[K];memcpy(w,A[bq],K);for(int j=0;j<4;j++)w[cols[t][j]]=v[j];
  for(int t2=0;t2<NT;t2++)if(smask[t2]&bcm){dec(idx(A[bq],t2));inc(idx(w,t2));}
  memcpy(A[bq],w,K);
  for(int c=0;c<K;c++)if((bcm>>c)&1)tabu[bq][c]=it+4+rnd()%10;
  if(nu<bestnu){bestnu=nu;}
 }
 return nu==0;}

/* ---------------- coset-union (linear code orbit) search ---------------- */
#define MAXM 64
#define MAXL 25
static int gd,ge,gLS,gm;
static int GM[3][K];
static int mulc[NT][4][2][5];
static int gcnt[NT*MAXL];static unsigned int gW[NT*MAXL];
static int gul[NT*MAXL],gup[NT*MAXL],gnu;
static unsigned char SH[MAXM][K];
static int glab[MAXM][NT];
static long long gtabu[MAXM][K];
static const int INV5[5]={0,1,3,2,4};
static int pow5(int d){int r=1;while(d--)r*=5;return r;}
static int nullsp(int P[3][4],int d,int basis[4][4]){
 int M[3][4];memcpy(M,P,sizeof(M));int r=0,pc[3];
 for(int c=0;c<4&&r<d;c++){int p=-1;for(int i=r;i<d;i++)if(M[i][c]){p=i;break;}
  if(p<0)continue;
  for(int k=0;k<4;k++){int tmp=M[p][k];M[p][k]=M[r][k];M[r][k]=tmp;}
  int iv=INV5[M[r][c]];for(int k=0;k<4;k++)M[r][k]=M[r][k]*iv%5;
  for(int i=0;i<d;i++)if(i!=r&&M[i][c]){int f=M[i][c];for(int k=0;k<4;k++)M[i][k]=((M[i][k]-f*M[r][k])%5+25)%5;}
  pc[r]=c;r++;}
 int isp[4]={0,0,0,0};for(int i=0;i<r;i++)isp[pc[i]]=1;
 int nb=0;for(int f=0;f<4;f++)if(!isp[f]){int v[4]={0,0,0,0};v[f]=1;for(int i=0;i<r;i++)v[pc[i]]=(5-M[i][f])%5;memcpy(basis[nb++],v,sizeof(v));}
 return r;}
static void setupG(int d){
 gd=d;ge=4-d;gLS=pow5(ge);
 for(;;){
  for(int i=0;i<d;i++)for(int c=0;c<K;c++)GM[i][c]=rnd()%5;
  int ok=1;
  for(int t=0;t<NT&&ok;t++){int P[3][4];memset(P,0,sizeof P);
   for(int i=0;i<d;i++)for(int j=0;j<4;j++)P[i][j]=GM[i][cols[t][j]];
   int B[4][4];memset(B,0,sizeof B);int r=nullsp(P,d,B);
   if(r<d){ok=0;break;}
   for(int j=0;j<4;j++)for(int comp=0;comp<2;comp++)for(int v=0;v<5;v++)
    mulc[t][j][comp][v]=(comp<ge)?(B[comp][j]*v)%5:0;}
  if(ok)return;}}
static inline int labOf(const unsigned char*s,int t){int l0=0,l1=0;
 for(int j=0;j<4;j++){int v=s[cols[t][j]];l0+=mulc[t][j][0][v];l1+=mulc[t][j][1][v];}
 return ge==2?(l0%5)*5+(l1%5):l0%5;}
static inline void ginc(int k){if(gcnt[k]++==0){int p=gup[k];int l=gul[--gnu];gul[p]=l;gup[l]=p;gup[k]=-1;}}
static inline void gdec(int k){if(--gcnt[k]==0){gup[k]=gnu;gul[gnu++]=k;}}
static void gbuild(){int NTL=NT*gLS;memset(gcnt,0,sizeof gcnt);gnu=0;
 for(int q=0;q<gm;q++)for(int t=0;t<NT;t++){glab[q][t]=labOf(SH[q],t);gcnt[t*gLS+glab[q][t]]++;}
 for(int k=0;k<NTL;k++){if(!gcnt[k]){gup[k]=gnu;gul[gnu++]=k;}else gup[k]=-1;}}
static int gsearch(double tend,int*bestout){
 long long it=0;int bestnu=gnu;int NTL=NT*gLS;
 for(int k=0;k<NTL;k++)gW[k]=1;
 memset(gtabu,0,sizeof(gtabu[0])*gm);
 while(gnu>0){
  if((it&255)==0&&now()>tend)break;
  it++;
  if(it%20000==0)for(int k=0;k<NTL;k++)gW[k]=1+(gW[k]-1)*3/4;
  int k=gul[rnd()%gnu];int t=k/gLS,L=k%gLS;int L0,L1;
  if(ge==2){L0=L/5;L1=L%5;}else{L0=L;L1=0;}
  long long bw=(1LL<<62);int bq=-1,bc=0,bcm=0;unsigned char brow[K];
  for(int q=0;q<gm;q++){
   int a=glab[q][t];int a0,a1;if(ge==2){a0=a/5;a1=a%5;}else{a0=a;a1=0;}
   int D0=(L0-a0+5)%5,D1=(L1-a1+5)%5;
   for(int j1=0;j1<4;j1++)for(int j2=j1;j2<4;j2++){
    int single=(j2==j1);
    for(int x1=1;x1<5;x1++)for(int x2=(single?0:1);x2<(single?1:5);x2++){
     int s0=mulc[t][j1][0][x1]+(single?0:mulc[t][j2][0][x2]);
     int s1=mulc[t][j1][1][x1]+(single?0:mulc[t][j2][1][x2]);
     if(s0%5!=D0||s1%5!=D1)continue;
     unsigned char w[K];memcpy(w,SH[q],K);
     int c1=cols[t][j1],c2=cols[t][j2];int cm=1<<c1;
     w[c1]=(w[c1]+x1)%5;if(!single){w[c2]=(w[c2]+x2)%5;cm|=1<<c2;}
     long long dw=0;int du=0;
     for(int t2=0;t2<NT;t2++)if(smask[t2]&cm){int aa=t2*gLS+glab[q][t2];int bb=t2*gLS+labOf(w,t2);
      if(aa==bb)continue;if(gcnt[aa]==1){dw+=gW[aa];du++;}if(gcnt[bb]==0){dw-=gW[bb];du--;}}
     int ist=0;for(int c=0;c<K;c++)if(((cm>>c)&1)&&gtabu[q][c]>it){ist=1;break;}
     if(ist&&gnu+du>=bestnu)continue;
     if(dw<bw){bw=dw;bq=q;bcm=cm;bc=1;memcpy(brow,w,K);}
     else if(dw==bw&&rnd()%(++bc)==0){bq=q;bcm=cm;memcpy(brow,w,K);}
    }}
  }
  if(bq<0)continue;
  if(bw>=0)for(int u=0;u<gnu;u++)gW[gul[u]]++;
  for(int t2=0;t2<NT;t2++)if(smask[t2]&bcm){int nl=labOf(brow,t2);
   if(nl!=glab[bq][t2]){gdec(t2*gLS+glab[bq][t2]);ginc(t2*gLS+nl);glab[bq][t2]=nl;}}
  memcpy(SH[bq],brow,K);
  for(int c=0;c<K;c++)if((bcm>>c)&1)gtabu[bq][c]=it+3+rnd()%6;
  if(gnu<bestnu)bestnu=gnu;
 }
 *bestout=bestnu;
 return gnu==0;}
static int expandAndStore(){
 int nu5=pow5(gd);int n=0;
 for(int q=0;q<gm;q++)for(int u=0;u<nu5;u++){int uu[3],x=u;for(int i=0;i<gd;i++){uu[i]=x%5;x/=5;}
  for(int c=0;c<K;c++){int v=SH[q][c];for(int i=0;i<gd;i++)v+=uu[i]*GM[i][c];X[n][c]=v%5;}n++;}
 if(!verifyArr(X,n)){fprintf(stderr,"expand verify FAILED\n");return 0;}
 memcpy(A,X,(size_t)n*K);nr=n;build();
 for(int q=nr-1;q>=0;q--)if(q<nr&&remdelta(q)==0)delrow(q);
 if(!verifyArr(A,nr))return 0;
 fprintf(stderr,"group solution: %d rows, pruned to %d\n",n,nr);
 if(bestN==0||nr<bestN){memcpy(BEST,A,(size_t)nr*K);bestN=nr;fprintf(stderr,"NEW BEST N=%d t=%.0f\n",nr,now()-T0);return 1;}
 return 0;}
static void tryGroup(int d,int m0,double phaseEnd,double per){
 int m=m0;int P=pow5(d);long long att=0;int overall=1<<30;
 while(now()<phaseEnd){
  if(bestN>0&&m*P>=bestN){m=(bestN-1)/P;}
  if(m<=0||m>MAXM)return;
  setupG(d);gm=m;
  for(int q=0;q<gm;q++)for(int c=0;c<K;c++)SH[q][c]=rnd()%5;
  gbuild();
  double te=now()+per;if(te>phaseEnd)te=phaseEnd;
  int bn;int ok=gsearch(te,&bn);att++;
  if(bn<overall)overall=bn;
  if(d==2||att%200==0)fprintf(stderr,"d=%d m=%d att=%lld bestnu=%d overall=%d t=%.0f\n",d,m,att,bn,overall,now()-T0);
  if(ok){expandAndStore();m--;}
 }}
int main(){T0=now();rs^=(unsigned long long)time(0)*2654435761ULL;for(int i=0;i<10;i++)rnd();
 int t=0;for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){
  cols[t][0]=a;cols[t][1]=b;cols[t][2]=c;cols[t][3]=d;smask[t]=(1<<a)|(1<<b)|(1<<c)|(1<<d);t++;}
 const char*fn=getenv("BEST_CANDIDATE");nr=0;
 if(fn){FILE*fp=fopen(fn,"r");if(fp){int x[K];while(nr<MAXR){int ok=1;for(int c=0;c<K;c++)if(fscanf(fp,"%d",&x[c])!=1){ok=0;break;}
   if(!ok)break;for(int c=0;c<K;c++)A[nr][c]=(unsigned char)(((x[c]%5)+5)%5);nr++;}fclose(fp);}}
 if(nr>0&&verifyArr(A,nr)){memcpy(BEST,A,(size_t)nr*K);bestN=nr;fprintf(stderr,"loaded N=%d\n",nr);}
 /* Phase 1: d=3, 7 cosets of a [7,3] code -> 875 rows */
 tryGroup(3,7,T0+150,1.0);
 /* Phase 2: d=2, 36 cosets of a [7,2] code -> 900 rows (then fewer) */
 tryGroup(2,36,T0+2750,180.0);
 /* Phase 3: generic remove-and-repair on best */
 double TEND=T0+3450;int attempt=0;
 while(bestN>0&&now()<TEND-30){
  nr=bestN;memcpy(A,BEST,(size_t)nr*K);build();
  int rd[MAXR],md=1<<30;for(int q=0;q<nr;q++){rd[q]=remdelta(q);if(rd[q]<md)md=rd[q];}
  int lim=md+(attempt>0&&rnd()%2?1:0);int bq=-1,bc=0;
  for(int q=0;q<nr;q++)if(rd[q]<=lim&&rnd()%(++bc)==0)bq=q;
  delrow(bq);
  double te=now()+300;if(te>TEND)te=TEND;
  attempt++;
  if(search(te)){
   for(int q=nr-1;q>=0;q--)if(q<nr&&remdelta(q)==0)delrow(q);
   if(nu==0&&verifyArr(A,nr)&&nr<bestN){memcpy(BEST,A,(size_t)nr*K);bestN=nr;
    fprintf(stderr,"new best N=%d t=%.0f\n",nr,now()-T0);attempt=0;}
  }
 }
 for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",BEST[r][c]);printf("\n");}
 fflush(stdout);return 0;}
