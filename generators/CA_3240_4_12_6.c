#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#define K 12
#define T4 1296
#define NT (495*T4)
#define MAXN 5000
#define MAXM 12
#define MAXU (MAXM*495+8)
static volatile sig_atomic_t stopf=0;
static void onterm(int s){(void)s;stopf=1;}
static int N; static unsigned char A[MAXN][K],bestA[MAXN][K],inputA[MAXN][K]; static int bestN=0,inputN=0;
static int sets[495][4],smask[495],ns=0;
static int cs[K][165],cw[K][165],ncs[K];
static uint16_t cnt[NT]; static int missList[NT],missPos[NT],nmiss=0;
static uint64_t rs=88172645463325252ULL;
static inline uint64_t rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline int tidx(int r,int s){int*c=sets[s];return s*T4+A[r][c[0]]+6*A[r][c[1]]+36*A[r][c[2]]+216*A[r][c[3]];}
static inline void addM(int t){missPos[t]=nmiss;missList[nmiss++]=t;}
static inline void remM(int t){int p=missPos[t];int l=missList[--nmiss];missList[p]=l;missPos[l]=p;}
static uint16_t xr[NT]; static uint16_t uq[MAXN][K]; static int uqt[MAXN];
static uint64_t colBits[K][6][MAXN/64+1];
static inline void own(int r,int t,int d){int*c=sets[t/T4];
  uq[r][c[0]]+=d;uq[r][c[1]]+=d;uq[r][c[2]]+=d;uq[r][c[3]]+=d;uqt[r]+=d;}
static inline void inc(int t,int r){int n=cnt[t]++;
  if(n==0){remM(t);xr[t]=(uint16_t)r;own(r,t,1);}
  else if(n==1){own(xr[t],t,-1);xr[t]^=(uint16_t)r;}
  else xr[t]^=(uint16_t)r;}
static inline void dec(int t,int r){int n=--cnt[t];xr[t]^=(uint16_t)r;
  if(n==0){addM(t);own(r,t,-1);}
  else if(n==1)own(xr[t],t,1);}
static inline void putA(int r,int c,int v){int a=A[r][c];if(a==v)return;
  colBits[c][a][r>>6]&=~(1ULL<<(r&63));colBits[c][v][r>>6]|=1ULL<<(r&63);A[r][c]=(unsigned char)v;}
static void rebuild(){memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);memset(uq,0,sizeof uq);
  memset(uqt,0,sizeof uqt);memset(colBits,0,sizeof colBits);nmiss=0;
  for(int r=0;r<N;r++)for(int c=0;c<K;c++)colBits[c][A[r][c]][r>>6]|=1ULL<<(r&63);
  for(int r=0;r<N;r++)for(int s=0;s<ns;s++){int t=tidx(r,s);cnt[t]++;xr[t]^=(uint16_t)r;}
  for(int t=0;t<NT;t++){if(!cnt[t])addM(t);else if(cnt[t]==1)own(xr[t],t,1);}}
static inline void decode(int t,int*cols,int*vals){int s=t/T4,v=t%T4;for(int i=0;i<4;i++){cols[i]=sets[s][i];vals[i]=v%6;v/=6;}}
static void delRow(int r){for(int s=0;s<ns;s++)dec(tidx(r,s),r);
  int l=N-1;
  for(int c=0;c<K;c++)colBits[c][A[r][c]][r>>6]&=~(1ULL<<(r&63));
  if(r!=l){
    for(int c=0;c<K;c++){colBits[c][A[l][c]][r>>6]|=1ULL<<(r&63);colBits[c][A[l][c]][l>>6]&=~(1ULL<<(l&63));}
    for(int s=0;s<ns;s++)xr[tidx(l,s)]^=(uint16_t)(l^r);
    memcpy(uq[r],uq[l],sizeof uq[0]);uqt[r]=uqt[l];
    memcpy(A[r],A[l],K);}
  memset(uq[l],0,sizeof uq[0]);uqt[l]=0;
  N--;}
static int weak[32],nweak=0;
static int weakest(){int br=0,bu=1<<30,ties=0;for(int r=0;r<N;r++){int u=uqt[r];
  if(u<bu){bu=u;br=r;ties=1;}else if(u==bu&&rnd()%(++ties)==0)br=r;}return br;}
static int wthr=0;
static void calcWeak(){int wu[32];nweak=0;for(int r=0;r<N;r++){int u=uqt[r];
  int p=nweak<32?nweak:32;if(p==32&&u>=wu[31])continue;if(p==32)p=31;
  while(p>0&&wu[p-1]>u){wu[p]=wu[p-1];weak[p]=weak[p-1];p--;}
  wu[p]=u;weak[p]=r;if(nweak<32)nweak++;}
  wthr=nweak?wu[nweak-1]:0;}
/* colouring */
static int ut[MAXU],uc[MAXU][4],uv[MAXU][4],nu,mcl;
static uint16_t mask_[MAXU]; static char done_[MAXU];
static signed char cls[MAXM][K]; static int cn[MAXM];
static int trU[MAXU*MAXM+16]; static uint16_t trM[MAXU*MAXM+16]; static int trail;
static long nodes,nodeLim; static int abortf,stuck,bestLeft;
static long statNodes=0;
static int drainF=0;
static int rec(int left){
  if(left==0)return 1;
  if(++nodes>nodeLim){abortf=1;return 0;}
  statNodes++;
  int u=-1,bp=99,off=rnd()%nu;
  for(int ii=0;ii<nu;ii++){int i=ii+off;if(i>=nu)i-=nu;if(done_[i])continue;
    int pc=__builtin_popcount(mask_[i]);if(pc<bp){bp=pc;u=i;if(pc==1)break;}}
  int opts=mask_[u],seenEmpty=0,st=rnd()%mcl;
  for(int kk=0;kk<mcl;kk++){int k;
    if(drainF&&mcl>1)k=kk<mcl-1?1+(kk+st)%(mcl-1):0;
    else{k=kk+st;if(k>=mcl)k-=mcl;}
    if(!(opts>>k&1))continue;
    if(cn[k]==0){if(seenEmpty)continue;seenEmpty=1;}
    signed char save[K];memcpy(save,cls[k],K);int mark=trail;
    for(int i=0;i<4;i++)cls[k][uc[u][i]]=uv[u][i];
    cn[k]++;done_[u]=1;int ok=1;
    for(int v=0;v<nu;v++){if(done_[v]||!(mask_[v]>>k&1))continue;
      int bad=0;for(int i=0;i<4;i++){int x=cls[k][uc[v][i]];if(x>=0&&x!=uv[v][i]){bad=1;break;}}
      if(bad){trU[trail]=v;trM[trail]=mask_[v];trail++;mask_[v]&=~(1u<<k);
        if(!mask_[v]){ok=0;if(left<bestLeft){bestLeft=left;stuck=v;}break;}}}
    if(ok&&rec(left-1))return 1;
    while(trail>mark){trail--;mask_[trU[trail]]=trM[trail];}
    memcpy(cls[k],save,K);cn[k]--;done_[u]=0;
    if(abortf)return 0;}
  return 0;}
static int solve(int hasHole){
  for(int i=0;i<nu;i++){mask_[i]=(uint16_t)((1u<<mcl)-1);done_[i]=0;}
  if(hasHole)mask_[0]=1;
  memset(cls,-1,sizeof cls);memset(cn,0,sizeof cn);trail=0;nodes=0;abortf=0;stuck=-1;bestLeft=1<<30;
  nodeLim=4*nu+600;
  return rec(nu);}
/* row selection */
static int R[MAXM],nR; static int inR[MAXN],inRstamp=0;
static int Q[MAXM*495+8],nQ;
static int stampT[NT],subT[NT],curStamp=0;
static void pushRow(int r,int dm){R[nR++]=r;inR[r]=inRstamp;
  if(dm)for(int s=0;s<ns;s++)if((smask[s]&dm)&&cnt[tidx(r,s)]==1)Q[nQ++]=tidx(r,s);}
static void addRowFor(int q){int cols[4],vals[4];decode(q,cols,vals);
  int c1[160],c1m[160],n1=0,c2[48],n2=0,nw=(N+63)/64,o=rnd()%nw;
  uint64_t*ba=colBits[cols[0]][vals[0]],*bb=colBits[cols[1]][vals[1]];
  uint64_t*bc=colBits[cols[2]][vals[2]],*bd=colBits[cols[3]][vals[3]];
  for(int i=0;i<nw;i++){int w=o+i;if(w>=nw)w-=nw;
    uint64_t x=ba[w],y=bb[w],z=bc[w],u=bd[w],ab=x&y,cd=z&u,xy=x|y,zu=z|u;
    uint64_t m3=(ab&zu)|(cd&xy),m2=(ab|cd|(xy&zu))&~m3;
    while(m3&&n1<112){int r=(w<<6)+__builtin_ctzll(m3);m3&=m3-1;if(inR[r]==inRstamp)continue;
      int dm=0;for(int j=0;j<4;j++)if(A[r][cols[j]]!=vals[j])dm|=1<<cols[j];
      c1[n1]=r;c1m[n1]=dm;n1++;}
    while(m2&&n2<48){int r=(w<<6)+__builtin_ctzll(m2);m2&=m2-1;if(inR[r]==inRstamp)continue;
      c2[n2++]=r;}}
  if(n1<8||rnd()%5==0){for(int i=0;i<n2&&n1<160;i++){int r=c2[i],dm=0;
      for(int j=0;j<4;j++)if(A[r][cols[j]]!=vals[j])dm|=1<<cols[j];
      c1[n1]=r;c1m[n1]=dm;n1++;}}
  int nc=n1;int*cr=c1,*cm=c1m;
  if(nc==0)return;
  int br=-1,bm=0,bs=1<<30;
  for(int i=0;i<nc;i++){int r=cr[i],dm=cm[i],sc=0;
    for(int e=dm;e;e&=e-1)sc+=uq[r][__builtin_ctz(e)];
    sc=sc*2+(int)(rnd()%6);
    if(sc<bs){bs=sc;br=r;bm=dm;}}
  pushRow(br,bm);}
static void buildR(int seed,int m){nR=0;nQ=0;inRstamp++;
  if(seed>=0)addRowFor(seed);
  else{int r=(nweak&&rnd()%2)?weak[rnd()%nweak]:(int)(rnd()%N);pushRow(r,0);
    for(int s=0;s<ns;s++)if(cnt[tidx(r,s)]==1&&nQ<495)Q[nQ++]=tidx(r,s);}
  while(nR<m){
    if(nweak&&rnd()%100<15){int r=weak[rnd()%nweak];if(inR[r]!=inRstamp){pushRow(r,0);continue;}}
    if(nQ==0){int r=rnd()%N;if(inR[r]!=inRstamp)pushRow(r,0);
      for(int s=0;s<ns&&nQ<495;s++)if(cnt[tidx(r,s)]==1)Q[nQ++]=tidx(r,s);
      if(nQ==0)continue;}
    int i=rnd()%nQ;int q=Q[i];Q[i]=Q[--nQ];
    int before=nR;addRowFor(q);
    if(nR==before&&nQ==0){int r=rnd()%N;if(inR[r]!=inRstamp)pushRow(r,0);}}}
static inline void addU(int t){int cols[4],vals[4];decode(t,cols,vals);ut[nu]=t;
  for(int i=0;i<4;i++){uc[nu][i]=cols[i];uv[nu][i]=vals[i];}nu++;}
static void buildU(int hole){nu=0;if(hole>=0)addU(hole);curStamp++;
  for(int k=0;k<nR;k++)for(int s=0;s<ns;s++){int t=tidx(R[k],s);
    if(stampT[t]!=curStamp){stampT[t]=curStamp;subT[t]=0;}
    if(++subT[t]==cnt[t])addU(t);}}
static void commit(){
  for(int k=0;k<nR;k++)for(int s=0;s<ns;s++)dec(tidx(R[k],s),R[k]);
  for(int k=0;k<nR;k++)for(int c=0;c<K;c++)putA(R[k],c,cls[k][c]>=0?cls[k][c]:(int)(rnd()%6));
  for(int k=0;k<nR;k++)for(int s=0;s<ns;s++)inc(tidx(R[k],s),R[k]);
  for(int k=0;k<nR;k++){int r=R[k],st=rnd()%K;
    for(int cc=0;cc<K;cc++){int c=cc+st;if(c>=K)c-=K;if(cls[k][c]>=0)continue;
      int a=A[r][c],sc[6]={0,0,0,0,0,0};
      for(int j=0;j<ncs[c];j++){int o=tidx(r,cs[c][j]);
        for(int b=0;b<6;b++){int n=o+(b-a)*cw[c][j];int val=(int)cnt[n]-(b==a);
          if(val==0)sc[b]+=1000;
          else if(val==1){int ow=b==a?(int)(xr[n]^(uint16_t)r):(int)xr[n];sc[b]+=uqt[ow]<=wthr?5:1;}}}
      if(rnd()%100<12)for(int b=0;b<6;b++)sc[b]/=1000;
      int bb=a,bv=-1,ties=0;for(int b=0;b<6;b++){if(sc[b]>bv){bv=sc[b];bb=b;ties=1;}else if(sc[b]==bv&&rnd()%(++ties)==0)bb=b;}
      if(bb!=a){for(int j=0;j<ncs[c];j++){int o=tidx(r,cs[c][j]);int n=o+(bb-a)*cw[c][j];dec(o,r);inc(n,r);}
        putA(r,c,bb);}}}}
static long hT[NT];
static unsigned char minA[MAXN][K];
static void setCell(int r,int c,int b){int a=A[r][c];if(a==b)return;
  for(int j=0;j<ncs[c];j++){int o=tidx(r,cs[c][j]);int n=o+(b-a)*cw[c][j];dec(o,r);inc(n,r);}
  putA(r,c,b);}
static void excursion(long att){int hole=missList[rnd()%nmiss],cols[4],vals[4];decode(hole,cols,vals);
  int br=-1,bs=1<<30;
  for(int r=0;r<N;r++){int dm=0;for(int j=0;j<4;j++)if(A[r][cols[j]]!=vals[j])dm|=1<<cols[j];
    int sc=(int)(rnd()%6);
    for(int s=0;s<ns&&sc<bs;s++)if(smask[s]&dm){int t=tidx(r,s);if(cnt[t]==1)sc+=hT[t]>att?12:4;}
    if(sc<bs){bs=sc;br=r;}}
  hT[hole]=att+20000+rnd()%20000;
  for(int j=0;j<4;j++)setCell(br,cols[j],vals[j]);}
static int verify(unsigned char(*X)[K],int n){static uint8_t seen[NT];memset(seen,0,sizeof seen);
  for(int r=0;r<n;r++)for(int s=0;s<ns;s++){int*c=sets[s];
    seen[s*T4+X[r][c[0]]+6*X[r][c[1]]+36*X[r][c[2]]+216*X[r][c[3]]]=1;}
  for(int t=0;t<NT;t++)if(!seen[t])return 0;return 1;}
static int readArr(const char*fn,unsigned char(*X)[K]){FILE*f=fn?fopen(fn,"r"):NULL;int n=0;if(!f)return 0;
  int x[K];while(n<MAXN){int ok=1;for(int i=0;i<K;i++)if(fscanf(f,"%d",&x[i])!=1){ok=0;break;}
    if(!ok)break;for(int i=0;i<K;i++)X[n][i]=x[i];n++;}fclose(f);return n;}
/* === orbit phase: union of row orbits of a symbol group, with subgroup descent === */
#define GM 2000
static int gG[120][6],gn,gcid[1296],gnc;
static uint16_t gcnt[NT];static int gh[NT],ghp[NT],gE;
static unsigned char gS[GM][K],gbS[GM][K];
static int gMask[GM],gbMask[GM],gCost[GM],gDom[GM][6],gDn[GM],gm,gbm,gbN;
static uint16_t gcode[GM][495];
static inline void gAddH(int u){ghp[u]=gE;gh[gE++]=u;}
static inline void gRemH(int u){int p=ghp[u];int l=gh[--gE];gh[p]=l;ghp[l]=p;}
static inline void gInc(int u){if(gcnt[u]++==0)gRemH(u);}
static inline void gDec(int u){if(--gcnt[u]==0)gAddH(u);}
static void gBuildGroup(int lv){static const int inv[5]={0,1,3,2,4};gn=0;
  for(int a=0;a<5;a++)for(int b=0;b<5;b++)for(int c=0;c<5;c++)for(int d=0;d<5;d++){
    int det=((a*d-b*c)%5+5)%5;if(!det)continue;
    if(lv==1&&det!=1&&det!=4)continue;
    if(lv>=2){if(c!=0)continue;int ra=a*inv[d]%5;
      if(lv==2&&ra!=1&&ra!=4)continue;
      if(lv==3&&ra!=1)continue;}
    int p[6];for(int x=0;x<5;x++){int den=(c*x+d)%5,num=(a*x+b)%5;p[x]=den?num*inv[den]%5:5;}
    p[5]=c?a*inv[c]%5:5;
    int dup=0;for(int g=0;g<gn&&!dup;g++)if(!memcmp(gG[g],p,sizeof p))dup=1;
    if(!dup){memcpy(gG[gn],p,sizeof p);gn++;}}
  for(int i=0;i<1296;i++)gcid[i]=-1;
  gnc=0;
  for(int code=0;code<1296;code++){if(gcid[code]>=0)continue;
    int d0=code%6,d1=code/6%6,d2=code/36%6,d3=code/216%6;
    for(int g=0;g<gn;g++){int*p=gG[g];gcid[p[d0]+6*p[d1]+36*p[d2]+216*p[d3]]=gnc;}
    gnc++;}}
static void gSetup(int i){int mk=gMask[i];gDn[i]=0;
  for(int b=0;b<6;b++)if(mk>>b&1)gDom[i][gDn[i]++]=b;
  int st=0;for(int g=0;g<gn;g++){int ok=1;for(int b=0;b<6;b++)if((mk>>b&1)&&gG[g][b]!=b)ok=0;st+=ok;}
  gCost[i]=gn/st;
  for(int s=0;s<ns;s++){int*c=sets[s];
    int code=gS[i][c[0]]+6*gS[i][c[1]]+36*gS[i][c[2]]+216*gS[i][c[3]];
    gcode[i][s]=(uint16_t)code;gInc(s*gnc+gcid[code]);}}
static void gLoad(){int U=495*gnc;memset(gcnt,0,sizeof(uint16_t)*U);gE=0;
  for(int u=0;u<U;u++)gAddH(u);
  for(int i=0;i<gm;i++)gSetup(i);}
static inline int gDelta(int i,int c,int b){int a=gS[i][c],d=0,dv=b-a;
  for(int j=0;j<ncs[c];j++){int s=cs[c][j],o=gcode[i][s],n=o+dv*cw[c][j];
    int k1=gcid[o],k2=gcid[n];if(k1==k2)continue;int base=s*gnc;
    if(gcnt[base+k1]==1)d++;
    if(gcnt[base+k2]==0)d--;}
  return d;}
static void gApply(int i,int c,int b){int a=gS[i][c],dv=b-a;
  for(int j=0;j<ncs[c];j++){int s=cs[c][j],o=gcode[i][s],n=o+dv*cw[c][j];
    int k1=gcid[o],k2=gcid[n];gcode[i][s]=(uint16_t)n;if(k1==k2)continue;int base=s*gnc;
    gInc(base+k2);gDec(base+k1);}
  gS[i][c]=(unsigned char)b;}
static void gRemove(int i){for(int s=0;s<ns;s++)gDec(s*gnc+gcid[gcode[i][s]]);
  int l=gm-1;
  if(i!=l){memcpy(gS[i],gS[l],K);gMask[i]=gMask[l];gCost[i]=gCost[l];gDn[i]=gDn[l];
    memcpy(gDom[i],gDom[l],sizeof gDom[0]);memcpy(gcode[i],gcode[l],sizeof gcode[0]);}
  gm--;}
static int gWeak(int div){int bi=-1;double bs=1e18;
  for(int i=0;i<gm;i++){if(div&&rnd()%2)continue;int u=0;
    for(int s=0;s<ns;s++)if(gcnt[s*gnc+gcid[gcode[i][s]]]==1)u++;
    double sc=(u+(double)(rnd()%100)*1e-3)/gCost[i];
    if(sc<bs){bs=sc;bi=i;}}
  if(bi<0)bi=(int)(rnd()%gm);
  return bi;}
static int gNow(){int n=0;for(int i=0;i<gm;i++)n+=gCost[i];return n;}
static int gCmp(const void*a,const void*b){long long x=*(const long long*)a,y=*(const long long*)b;return x<y?-1:x>y;}
static void gSaveBest(){int n=gNow();gbN=n;gbm=gm;
  memcpy(gbS,gS,(size_t)gm*K);memcpy(gbMask,gMask,sizeof(int)*gm);
  int lim=bestN>0?bestN:inputN;if(n>=lim)return;
  static long long codes[120*GM];static unsigned char X[MAXN][K];int nc=0;
  for(int i=0;i<gbm;i++)for(int g=0;g<gn;g++){long long v=0;
    for(int c=0;c<K;c++)v=v*6+gG[g][gbS[i][c]];
    codes[nc++]=v;}
  qsort(codes,nc,sizeof(long long),gCmp);
  int u=0;for(int i=0;i<nc;i++)if(i==0||codes[i]!=codes[i-1])codes[u++]=codes[i];
  if(u>=lim)return;
  for(int i=0;i<u;i++){long long v=codes[i];for(int c=K-1;c>=0;c--){X[i][c]=(unsigned char)(v%6);v/=6;}}
  if(verify(X,u)){memcpy(bestA,X,(size_t)u*K);bestN=u;fprintf(stderr,"ORBIT VALID N=%d (|G|=%d, starters=%d)\n",u,gn,gbm);}
  else fprintf(stderr,"orbit export FAILED verify, n=%d\n",u);}
static void gRun(double budget,double stall,double cap){double ts=now(),lastImp=ts;int tries=0;long moves=0;
  for(;;){
    if(stopf||now()-ts>budget||now()-lastImp>stall)break;
    if(gE==0){int n=gNow();
      if(gbN==0||n<gbN){gSaveBest();lastImp=now();
        fprintf(stderr,"|G|=%d valid cost=%d starters=%d t=%.0f moves=%ld\n",gn,n,gm,now()-ts,moves);}
      if(gm<=1)break;
      gRemove(gWeak(0));tries=0;continue;}
    if(gm==0)break;
    double Thi=0.16+0.03*tries;if(Thi>cap)Thi=cap;
    double Tlo=0.10;
    long L=(long)gm*K*5*200;if(L<1000000)L=1000000;
    L*=1+tries/2;
    uint32_t th[17];int quit=0;
    for(long st=0;st<L;st++){
      if((st&4095)==0){double T=Thi*pow(Tlo/Thi,(double)st/(double)L);
        for(int d=1;d<17;d++)th[d]=(uint32_t)(exp(-d/T)*4294967295.0);
        if((st&0xFFFFF)==0&&(stopf||now()-ts>budget)){quit=1;break;}}
      moves++;
      int i=-1,c=0,b=0;
      if((rnd()&7)==0){int u=gh[rnd()%gE],s=u/gnc,cl=u%gnc;
        for(int q=0;q<48&&i<0;q++){int ii=(int)(rnd()%gm);if(gDn[ii]<=1)continue;
          int j0=(int)(rnd()%4),code=gcode[ii][s];
          for(int jj=0;jj<4&&i<0;jj++){int j=(j0+jj)&3,w=j==0?1:j==1?6:j==2?36:216;
            int cc=sets[s][j],a=gS[ii][cc];
            for(int dk=0;dk<gDn[ii];dk++){int bb=gDom[ii][dk];if(bb==a)continue;
              if(gcid[code+(bb-a)*w]==cl){i=ii;c=cc;b=bb;break;}}}}
        if(i<0)continue;}
      else{i=(int)(rnd()%gm);if(gDn[i]<=1)continue;
        c=(int)(rnd()%K);b=gDom[i][rnd()%gDn[i]];if(b==gS[i][c])continue;}
      int d=gDelta(i,c,b);
      if(d<=0||(d<17&&(uint32_t)(rnd()>>32)<th[d])){gApply(i,c,b);if(gE==0)break;}}
    if(quit)break;
    if(gE==0)continue;
    tries++;
    if(tries>=10&&gbN>0){
      gm=gbm;memcpy(gS,gbS,(size_t)gm*K);memcpy(gMask,gbMask,sizeof(int)*gm);gLoad();
      if(gm>1)gRemove(gWeak(1));
      tries=0;}}}
static void gDescend(int lv){static int old[120][6];int oldn=gn;memcpy(old,gG,sizeof old);
  gBuildGroup(lv);
  static unsigned char nS[GM][K];static int nMask[GM];int nm=0;
  for(int i=0;i<gbm;i++){long long cd[120];unsigned char rows[120][K];int no=0;
    for(int g=0;g<oldn;g++){long long v=0;unsigned char t[K];
      for(int c=0;c<K;c++){t[c]=(unsigned char)old[g][gbS[i][c]];v=v*6+t[c];}
      int dup=0;for(int q=0;q<no;q++)if(cd[q]==v){dup=1;break;}
      if(!dup){cd[no]=v;memcpy(rows[no],t,K);no++;}}
    char vis[120];memset(vis,0,sizeof vis);
    for(int j=0;j<no&&nm<GM;j++){if(vis[j])continue;
      memcpy(nS[nm],rows[j],K);int mk=0;for(int c=0;c<K;c++)mk|=1<<rows[j][c];
      nMask[nm]=__builtin_popcount(mk)<=2?mk:63;nm++;
      for(int h=0;h<gn;h++){long long v=0;for(int c=0;c<K;c++)v=v*6+gG[h][rows[j][c]];
        for(int q=0;q<no;q++)if(cd[q]==v){vis[q]=1;break;}}}}
  gm=nm;memcpy(gS,nS,(size_t)nm*K);memcpy(gMask,nMask,sizeof(int)*nm);
  gbN=0;gLoad();
  fprintf(stderr,"descend: |G|=%d classes=%d starters=%d cost=%d holes=%d\n",gn,gnc,gm,gNow(),gE);}
static void groupPhase(){
  gBuildGroup(0);
  fprintf(stderr,"orbit phase: |G|=%d classes=%d\n",gn,gnc);
  gm=0;memset(gS[0],0,K);gMask[0]=1;gm=1;
  for(int i=0;i<16;i++){for(int c=0;c<K;c++)gS[gm][c]=(unsigned char)(rnd()>>40&1);gMask[gm]=3;gm++;}
  for(int i=0;i<40;i++){for(int c=0;c<K;c++)gS[gm][c]=(unsigned char)((rnd()>>20)%6);gMask[gm]=63;gm++;}
  gLoad();gbN=0;
  double bud[4]={300,300,450,450},stl[4]={150,150,200,200},cap[4]={0.45,0.40,0.30,0.28},lim[4]={1.15,1.08,1.03,9.0};
  for(int lv=0;lv<4;lv++){
    gRun(bud[lv],stl[lv],cap[lv]);
    fprintf(stderr,"level |G|=%d done: best cost=%d starters=%d global best=%d\n",gn,gbN,gbm,bestN);
    if(stopf||gbN==0)return;
    if(gbN>inputN*lim[lv])return;
    if(lv<3)gDescend(lv+1);}}
int main(){double t0=now(),TL=3450;signal(SIGTERM,onterm);
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){
    sets[ns][0]=a;sets[ns][1]=b;sets[ns][2]=c;sets[ns][3]=d;smask[ns]=1<<a|1<<b|1<<c|1<<d;
    int w=1;for(int i=0;i<4;i++){int col=sets[ns][i];cs[col][ncs[col]]=ns;cw[col][ncs[col]]=w;ncs[col]++;w*=6;}ns++;}
  rs^=(uint64_t)time(NULL)*2654435761ULL;for(int i=0;i<10;i++)rnd();
  inputN=readArr(getenv("BEST_CANDIDATE"),inputA);
  if(inputN==0)return 1;
  memcpy(A,inputA,sizeof(A[0])*inputN);N=inputN;
  if(verify(inputA,inputN)){memcpy(bestA,inputA,sizeof(A[0])*inputN);bestN=inputN;}
  groupPhase();
  fprintf(stderr,"orbit phase end t=%.0f bestN=%d inputN=%d\n",now()-t0,bestN,inputN);
  if(bestN>0&&bestN<inputN){memcpy(inputA,bestA,sizeof(A[0])*bestN);inputN=bestN;
    memcpy(A,inputA,sizeof(A[0])*inputN);N=inputN;}
  static unsigned char nearA[MAXN][K];int nn=readArr(getenv("NEAR_CANDIDATE"),nearA);
  if(nn>0&&nn<inputN){memcpy(A,nearA,sizeof(A[0])*nn);N=nn;rebuild();
    fprintf(stderr,"near candidate N=%d miss=%d\n",N,nmiss);
    if(nmiss>60){memcpy(A,inputA,sizeof(A[0])*inputN);N=inputN;rebuild();}}
  else rebuild();
  calcWeak();
  long att=0,okH=0,okMig=0,okPl=0,failH=0;int minMiss=1<<30;double lastLog=now();
  long lastImp=0,excs=0,restores=0;int stalls=0;
  while(!stopf){
    if((att&63)==0&&now()-t0>TL)break;
    if(nmiss==0){
      if(bestN==0||N<bestN){if(verify(A,N)){memcpy(bestA,A,sizeof(A[0])*N);bestN=N;
        fprintf(stderr,"VALID N=%d at %.0fs\n",N,now()-t0);}}
      int r=weakest();delRow(r);calcWeak();minMiss=nmiss;
      memcpy(minA,A,sizeof(A[0])*N);lastImp=att;stalls=0;
      fprintf(stderr,"deleted row: N=%d miss=%d\n",N,nmiss);
      continue;}
    att++;
    if(nmiss<minMiss){minMiss=nmiss;memcpy(minA,A,sizeof(A[0])*N);lastImp=att;stalls=0;
      fprintf(stderr,"N=%d miss=%d t=%.0f att=%ld\n",N,nmiss,now()-t0,att);}
    else if(nmiss==minMiss&&(att&1023)==0)memcpy(minA,A,sizeof(A[0])*N);
    if(att-lastImp>1000000){lastImp=att;stalls++;restores++;
      if(stalls>=8&&bestN==N+1){
        memcpy(A,bestA,sizeof(A[0])*bestN);N=bestN;rebuild();calcWeak();
        int r=weak[rnd()%(nweak<12?nweak:12)];delRow(r);calcWeak();minMiss=nmiss;
        memcpy(minA,A,sizeof(A[0])*N);stalls=0;
        fprintf(stderr,"switch deleted row: N=%d miss=%d t=%.0f\n",N,nmiss,now()-t0);
        continue;}
      else if(nmiss>minMiss){memcpy(A,minA,sizeof(A[0])*N);rebuild();calcWeak();}}
    if((att&0x3FFF)==0&&att-lastImp>32768&&nmiss<=minMiss+2){excursion(att);excs++;if(nmiss==0)continue;}
    if((att&0xFFF)==0)calcWeak();
    if((att&255)==0&&now()-lastLog>60){lastLog=now();
      fprintf(stderr,"t=%.0f N=%d miss=%d min=%d att=%ld okH=%ld mig=%ld pl=%ld fail=%ld nodes=%ld exc=%ld rest=%ld wthr=%d\n",
        now()-t0,N,nmiss,minMiss,att,okH,okMig,okPl,failH,statNodes,excs,restores,wthr);}
    int m;unsigned x=rnd()%100;
    if(x<10)m=1;else if(x<35)m=2+rnd()%2;else if(x<80)m=4+rnd()%4;else m=8+rnd()%(MAXM-7);
    mcl=m;
    if(rnd()%100<15){
      buildR(-1,m);mcl=nR;if(mcl==0)continue;buildU(-1);
      drainF=mcl>1&&uqt[R[0]]<=wthr&&rnd()%4!=0;
      int okS=nu>0&&solve(0);drainF=0;
      if(okS){commit();okPl++;}
      else if(nu==0){memset(cls,-1,sizeof cls);commit();okPl++;}
      continue;}
    int hole=missList[rnd()%nmiss];
    buildR(hole,m);mcl=nR;if(mcl==0)continue;
    buildU(hole);
    if(solve(1)){commit();okH++;hT[hole]=att+1000+rnd()%4000;continue;}
    failH++;
    if(stuck>=1&&rnd()%100<60){int t2=ut[stuck];
      if(hT[t2]<=att){
        nu--;ut[stuck]=ut[nu];memcpy(uc[stuck],uc[nu],sizeof uc[0]);memcpy(uv[stuck],uv[nu],sizeof uv[0]);
        if(solve(1)){commit();okMig++;hT[hole]=att+1000+rnd()%4000;}}}
  }
  fprintf(stderr,"end t=%.0f N=%d miss=%d min=%d att=%ld okH=%ld mig=%ld pl=%ld\n",now()-t0,N,nmiss,minMiss,att,okH,okMig,okPl);
  if(nmiss==0&&(bestN==0||N<bestN)&&verify(A,N)){memcpy(bestA,A,sizeof(A[0])*N);bestN=N;}
  if(bestN==0||!verify(bestA,bestN)){memcpy(bestA,inputA,sizeof(A[0])*inputN);bestN=inputN;}
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",bestA[r][c]);printf("\n");}
  return 0;}
