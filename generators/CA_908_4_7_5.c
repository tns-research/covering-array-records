#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#define K 7
#define V 5
#define NT 35
#define MAXN 1200
static int cols[NT][4], smask[NT], N;
static unsigned char A[MAXN][K], BEST[MAXN][K]; static int bestN;
static unsigned short cnt[NT][625];
static int ulist[NT*625], upos[NT*625], nu;
static int tabu[MAXN][K];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timeval t;gettimeofday(&t,0);return t.tv_sec+t.tv_usec*1e-6;}
static inline int idx(const unsigned char*r,int t){return ((r[cols[t][0]]*5+r[cols[t][1]])*5+r[cols[t][2]])*5+r[cols[t][3]];}
static inline void inc(int t,int i){if(cnt[t][i]++==0){int k=t*625+i,p=upos[k];int l=ulist[--nu];ulist[p]=l;upos[l]=p;upos[k]=-1;}}
static inline void dec(int t,int i){if(--cnt[t][i]==0){int k=t*625+i;upos[k]=nu;ulist[nu++]=k;}}
static void build(){memset(cnt,0,sizeof cnt);nu=0;
 for(int r=0;r<N;r++)for(int t=0;t<NT;t++)cnt[t][idx(A[r],t)]++;
 for(int t=0;t<NT;t++)for(int i=0;i<625;i++){int k=t*625+i;if(!cnt[t][i]){upos[k]=nu;ulist[nu++]=k;}else upos[k]=-1;}}
static void delrow(int r){for(int t=0;t<NT;t++)dec(t,idx(A[r],t));
 if(r!=N-1){memcpy(A[r],A[N-1],K);} N--;}
static int uniqcov(int r){int u=0;for(int t=0;t<NT;t++)if(cnt[t][idx(A[r],t)]==1)u++;return u;}
static void pickdel(){ // delete a low-uniqueness row (random among minima)
 int bu=1<<30,br=0,c=0;for(int r=0;r<N;r++){int u=uniqcov(r);if(u<bu){bu=u;br=r;c=1;}else if(u==bu&&rnd()%(++c)==0)br=r;}
 if(rnd()%4==0){int r=rnd()%N;if(uniqcov(r)<=bu+2)br=r;}
 delrow(br);}
static void printbest(){for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",BEST[r][c]);printf("\n");}fflush(stdout);}
static int verify(unsigned char (*X)[K],int n){static unsigned char s[625];
 for(int t=0;t<NT;t++){memset(s,0,625);for(int r=0;r<n;r++)s[idx(X[r],t)]=1;for(int i=0;i<625;i++)if(!s[i])return 0;}return 1;}
int main(){double T0=now(),TL=3450;rs^=(unsigned long long)time(0);
 int t=0;for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){cols[t][0]=a;cols[t][1]=b;cols[t][2]=c;cols[t][3]=d;smask[t]=(1<<a)|(1<<b)|(1<<c)|(1<<d);t++;}
 const char*fn=getenv("BEST_CANDIDATE");N=0;
 if(fn){FILE*f=fopen(fn,"r");if(f){int x[K];while(N<MAXN){int ok=1;for(int c=0;c<K;c++)if(fscanf(f,"%d",&x[c])!=1){ok=0;break;}if(!ok)break;for(int c=0;c<K;c++)A[N][c]=(unsigned char)(((x[c]%V)+V)%V);N++;}fclose(f);}}
 if(N==0||!verify(A,N)){ // fallback: trivial full array is too big; build OA-ish random cover greedily
  N=0;for(int i=0;i<625*5&&N<MAXN;i++){} // should not happen; build simple valid array
  // simple: all 5^4 on first 4 cols x shifts won't guarantee; use greedy random
  N=0;memset(cnt,0,sizeof cnt);build();
  while(1){build();if(nu==0)break;int k=ulist[rnd()%nu];int tt=k/625,i=k%625;unsigned char r[K];for(int c=0;c<K;c++)r[c]=rnd()%V;int v[4]={i/125,(i/25)%5,(i/5)%5,i%5};for(int j=0;j<4;j++)r[cols[tt][j]]=v[j];memcpy(A[N++],r,K);if(N>=MAXN)break;}
 }
 memcpy(BEST,A,sizeof(A[0])*N);bestN=N;build();
 static unsigned char SAV[MAXN][K];int savN;
 pickdel();memcpy(SAV,A,sizeof(A[0])*N);savN=N;int savnu=nu; // best partial at this level
 long long it=0,lastimp=0;int bestnu=nu;
 static int cand[MAXN],cd[MAXN];
 while(1){
  if((it&1023)==0&&now()-T0>TL)break;
  it++;
  if(nu==0){ if(verify(A,N)){memcpy(BEST,A,sizeof(A[0])*N);bestN=N;fprintf(stderr,"found N=%d t=%.0f it=%lld\n",N,now()-T0,it);}
   build();pickdel();memcpy(SAV,A,sizeof(A[0])*N);savN=N;bestnu=savnu=nu;lastimp=it;memset(tabu,0,sizeof tabu);continue;}
  if(it-lastimp>3000000){ // restart
   if(rnd()%3==0){memcpy(A,BEST,sizeof(A[0])*bestN);N=bestN;build();pickdel();}
   else{memcpy(A,SAV,sizeof(A[0])*savN);N=savN;build();}
   bestnu=nu;lastimp=it;memset(tabu,0,sizeof tabu);continue;}
  int k=ulist[rnd()%nu];int tt=k/625,i=k%625;int v[4]={i/125,(i/25)%5,(i/5)%5,i%5};
  int nc=0,maxd=2;
  for(int r=0;r<N;r++){int d=0;for(int j=0;j<4;j++)d+=A[r][cols[tt][j]]!=v[j];if(d<=maxd){cand[nc]=r;cd[nc++]=d;}}
  int bd=1<<30,br=-1,bc=0;
  int rmove=(rnd()%100==0);
  for(int q=0;q<nc;q++){int r=cand[q];unsigned char nr[K];memcpy(nr,A[r],K);int cm=0,istabu=0;
   for(int j=0;j<4;j++){int c=cols[tt][j];if(nr[c]!=v[j]){nr[c]=v[j];cm|=1<<c;if(tabu[r][c]>it)istabu=1;}}
   int dl=0;for(int s=0;s<NT;s++)if(smask[s]&cm){dl+=(cnt[s][idx(A[r],s)]==1)-(cnt[s][idx(nr,s)]==0);}
   if(istabu&&nu+dl>=bestnu)continue;
   if(rmove){if(rnd()%(++bc)==0){bd=dl;br=r;}continue;}
   if(dl<bd){bd=dl;br=r;bc=1;}else if(dl==bd&&rnd()%(++bc)==0)br=r;}
  if(br<0)continue;
  {int r=br;unsigned char nr[K];memcpy(nr,A[r],K);int cm=0;
   for(int j=0;j<4;j++){int c=cols[tt][j];if(nr[c]!=v[j]){nr[c]=v[j];cm|=1<<c;tabu[r][c]=it+3+rnd()%8;}}
   for(int s=0;s<NT;s++)if(smask[s]&cm){dec(s,idx(A[r],s));inc(s,idx(nr,s));}
   memcpy(A[r],nr,K);}
  if(nu<bestnu){bestnu=nu;lastimp=it;memcpy(SAV,A,sizeof(A[0])*N);savN=N;savnu=nu;}
 }
 if(!verify(BEST,bestN))fprintf(stderr,"warning: best invalid\n");
 printbest();return 0;}
