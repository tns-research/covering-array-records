#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define K 7
#define V 6
#define NS 35
#define NT (NS*1296)
#define MAXN 1900
static int sc[NS][4], csets[K][20], ncs[K];
static unsigned char A[MAXN][K], B[MAXN][K]; static int N, BN;
static unsigned short cnt[NT]; static int unc[NT], pos[NT], nunc;
static long long tabu[MAXN][K];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int idx(const unsigned char*r,int s){return s*1296+r[sc[s][0]]*216+r[sc[s][1]]*36+r[sc[s][2]]*6+r[sc[s][3]];}
static inline void inc(int t){if(cnt[t]++==0){int p=pos[t];int l=unc[--nunc];unc[p]=l;pos[l]=p;pos[t]=-1;}}
static inline void dec(int t){if(--cnt[t]==0){pos[t]=nunc;unc[nunc++]=t;}}
static void build(){memset(cnt,0,sizeof cnt);nunc=0;for(int t=0;t<NT;t++){pos[t]=nunc;unc[nunc++]=t;}
 for(int r=0;r<N;r++)for(int s=0;s<NS;s++)inc(idx(A[r],s));}
static int delta(int r,int c,int nv){unsigned char row[K];memcpy(row,A[r],K);int d=0;
 for(int j=0;j<20;j++){int s=csets[c][j];if(cnt[idx(row,s)]==1)d++;}
 row[c]=nv;for(int j=0;j<20;j++){int s=csets[c][j];if(cnt[idx(row,s)]==0)d--;}return d;}
static void setcell(int r,int c,int nv){for(int j=0;j<20;j++)dec(idx(A[r],csets[c][j]));A[r][c]=nv;
 for(int j=0;j<20;j++)inc(idx(A[r],csets[c][j]));}
static int rowloss(int r){int l=0;for(int s=0;s<NS;s++)if(cnt[idx(A[r],s)]==1)l++;return l;}
static void removerow(int r){for(int s=0;s<NS;s++)dec(idx(A[r],s));memcpy(A[r],A[N-1],K);N--;}
static int pickrow(int randomize){int best=1<<30,br=0;for(int r=0;r<N;r++){int l=rowloss(r)*16+(randomize?(int)(rnd()%(16*randomize)):(int)(rnd()%16));if(l<best){best=l;br=r;}}return br;}
int main(){
 int s=0;for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){
  sc[s][0]=a;sc[s][1]=b;sc[s][2]=c;sc[s][3]=d;csets[a][ncs[a]++]=s;csets[b][ncs[b]++]=s;csets[c][ncs[c]++]=s;csets[d][ncs[d]++]=s;s++;}
 double T0=now(),TL=3480;rs^=(unsigned long long)time(NULL);
 char*f=getenv("BEST_CANDIDATE");N=0;
 if(f){FILE*F=fopen(f,"r");if(F){int x[K];while(N<MAXN){int ok=1;for(int k=0;k<K;k++)if(fscanf(F,"%d",&x[k])!=1){ok=0;break;}if(!ok)break;for(int k=0;k<K;k++)A[N][k]=x[k]%V;N++;}fclose(F);}}
 build();
 if(nunc>0||N==0){/* fallback: shouldn't happen */ }
 memcpy(B,A,sizeof(A[0])*N);BN=N;
 removerow(pickrow(0));
 long long it=0,lastimp=0;int bestunc=nunc,restarts=0;
 while(1){
  if((it&1023)==0&&now()-T0>TL)break;
  it++;
  if(nunc==0){memcpy(B,A,sizeof(A[0])*N);BN=N;fprintf(stderr,"found N=%d at %.0fs\n",N,now()-T0);
   removerow(pickrow(0));bestunc=nunc;lastimp=it;continue;}
  if(nunc<bestunc){bestunc=nunc;lastimp=it;}
  if(it-lastimp>3000000){restarts++;N=BN;memcpy(A,B,sizeof(A[0])*N);build();removerow(pickrow(1+restarts%4));
   for(int p=0;p<20;p++){int r=rnd()%N,c=rnd()%K;setcell(r,c,rnd()%V);}
   bestunc=nunc;lastimp=it;memset(tabu,0,sizeof tabu);continue;}
  int t=unc[rnd()%nunc];int ss=t/1296,code=t%1296;int v[4]={code/216,(code/36)%6,(code/6)%6,code%6};
  int bd=1<<30,br=-1,bc=-1,nb=0;int cr[64],cc[64],nc=0;int r2=-1;
  for(int r=0;r<N;r++){int m=0,miss=-1;for(int j=0;j<4;j++){if(A[r][sc[ss][j]]==v[j])m++;else miss=j;}
   if(m==3){int c=sc[ss][miss];if(nc<64){cr[nc]=r;cc[nc]=c;nc++;}
    if(tabu[r][c]>it)continue;int d=delta(r,c,v[miss]);
    if(d<bd){bd=d;br=r;bc=c;nb=1;}else if(d==bd){nb++;if(rnd()%nb==0){br=r;bc=c;}}}
   else if(m==2&&(r2<0||rnd()%8==0))r2=r;}
  if(nc>0&&(br<0||rnd()%100<3)){int k=rnd()%nc;br=cr[k];bc=cc[k];}
  if(br>=0){int j;for(j=0;j<4;j++)if(sc[ss][j]==bc)break;setcell(br,bc,v[j]);tabu[br][bc]=it+5+rnd()%10;}
  else if(r2>=0){for(int j=0;j<4;j++){int c=sc[ss][j];if(A[r2][c]!=v[j]){setcell(r2,c,v[j]);tabu[r2][c]=it+5+rnd()%10;}}}
 }
 fprintf(stderr,"done N=%d it=%lld unc=%d\n",BN,it,nunc);
 for(int r=0;r<BN;r++){for(int k=0;k<K;k++)printf(k?" %d":"%d",B[r][k]);printf("\n");}
 return 0;}
