#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define K 27
#define NS 17550
#define NT 1296
#define NG (NS*NT)
#define CAP 7000
#define PC 2600
#define EVMAX 200
#define DUCAP 2
#define MARG 12
#define UPH 5
#define JCAP 100000
#define POLSEC 40.0
#define NCAP 40000
#define CHCAP 3000
#define CHLIM 2
#define DMAX 16
#define CHFRAC 0.5
#define CHMAXM 60
static int N;
static unsigned char A[CAP][K], B[CAP][K], I0[CAP][K];
static unsigned short cnt[NG], xr[NG];
static int w[NG];
static int mis[NG],pos[NG],M=0;
static int hv[NG],hpos[NG],H=0;
static int L[CAP][K],Lc[CAP][K];
static long long tb[CAP][K];
static int scol[NS][4];
static unsigned char mdc[NG][8];
typedef struct{int g,s;unsigned char o0,o1,o2,pad;unsigned short w0,w1,w2,wc;}E;
static E ce[K][PC];
static unsigned long long rs=0x9E3779B97F4A7C15ULL;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(int r,int s){unsigned char*x=A[r];return x[scol[s][0]]*216+x[scol[s][1]]*36+x[scol[s][2]]*6+x[scol[s][3]];}
static int U[CAP][NS],Uc[CAP],upos[NG];
static unsigned smask[NS];
static inline void addL(int r,int s,int dw,int dc,int g){int*sc=scol[s];
  if(dc>0){upos[g]=Uc[r];U[r][Uc[r]++]=g;}
  else if(dc<0){int p=upos[g];int l=U[r][--Uc[r]];U[r][p]=l;upos[l]=p;}
  L[r][sc[0]]+=dw;L[r][sc[1]]+=dw;L[r][sc[2]]+=dw;L[r][sc[3]]+=dw;
  Lc[r][sc[0]]+=dc;Lc[r][sc[1]]+=dc;Lc[r][sc[2]]+=dc;Lc[r][sc[3]]+=dc;}
static int oi[PC];
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  unsigned char*x=A[r]; int dd=b-a;
  E*ec=ce[c];
  for(int j=0;j<PC;j++){E*e=&ec[j];oi[j]=e->g+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2+a*e->wc;}
  for(int j=0;j<PC;j++){
    if(j+8<PC){int p=oi[j+8];int pn=p+dd*ec[j+8].wc;
      __builtin_prefetch(&cnt[p]);__builtin_prefetch(&xr[p]);__builtin_prefetch(&cnt[pn]);__builtin_prefetch(&xr[pn]);}
    int o=oi[j],n=o+dd*ec[j].wc,s=ec[j].s;
    unsigned short cO=--cnt[o]; xr[o]^=(unsigned short)r;
    if(cO==0){addm(o);addL(r,s,-w[o],-1,o);}
    else if(cO==1){addL(xr[o],s,w[o],1,o);}
    unsigned short cN=cnt[n]++;
    if(cN==0){delm(n);addL(r,s,w[n],1,n);}
    else if(cN==1){addL(xr[n],s,-w[n],-1,n);}
    xr[n]^=(unsigned short)r;
  }
  A[r][c]=b;
}
static void rebuild(){
  memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);
  for(int i=0;i<NG;i++)pos[i]=-1;
  M=0;
  for(int s=0;s<NS;s++){int*sc=scol[s];int g0=s*NT;
    for(int r=0;r<N;r++){unsigned char*x=A[r];int g=g0+x[sc[0]]*216+x[sc[1]]*36+x[sc[2]]*6+x[sc[3]];cnt[g]++;xr[g]^=(unsigned short)r;}}
  memset(L,0,sizeof L);memset(Lc,0,sizeof Lc);memset(Uc,0,sizeof Uc);
  for(int g=0;g<NG;g++){if(cnt[g]==0)addm(g);else if(cnt[g]==1)addL(xr[g],g/NT,w[g],1,g);}
}
static void resetw(){for(int i=0;i<H;i++){w[hv[i]]=1;hpos[hv[i]]=-1;}H=0;}
static void decayw(int num,int den){
  for(int i=0;i<H;){
    int g=hv[i];int ow=w[g];
    int nw=1+(int)((long long)(ow-1)*num/den); if(nw<1)nw=1;
    if(nw!=ow&&cnt[g]==1)addL(xr[g],g/NT,nw-ow,0,g);
    w[g]=nw;
    if(nw<=1){hpos[g]=-1;int l=hv[--H];hv[i]=l;hpos[l]=i;}
    else i++;
  }
}
static void delrow(int br){
  for(int s=0;s<NS;s++){int g=s*NT+tidx(br,s);unsigned short c=--cnt[g];xr[g]^=(unsigned short)br;
    if(c==0)addm(g);else if(c==1)addL(xr[g],s,w[g],1,g);}
  int last=N-1;
  if(br!=last){
    memcpy(A[br],A[last],K);memcpy(L[br],L[last],sizeof L[0]);memcpy(Lc[br],Lc[last],sizeof Lc[0]);
    memcpy(U[br],U[last],Uc[last]*sizeof(int));Uc[br]=Uc[last];
    unsigned short f=(unsigned short)(last^br);
    for(int s=0;s<NS;s++){int g=s*NT+tidx(br,s);xr[g]^=f;}
  }
  Uc[last]=0;
  N--;
}
static int pickrow(int rank){
  static int us[CAP],ord[CAP];
  for(int r=0;r<N;r++){int u=0;for(int c=0;c<K;c++)u+=Lc[r][c];us[r]=u/4;ord[r]=r;}
  for(int i=N-1;i>0;i--){int j=rnd()%(i+1);int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  if(rank>=N)rank=0;
  for(int i=0;i<=rank;i++){
    int bi=i;for(int j=i+1;j<N;j++)if(us[ord[j]]<us[ord[bi]])bi=j;
    int t=ord[i];ord[i]=ord[bi];ord[bi]=t;
  }
  int r=ord[rank];
  fprintf(stderr,"delete row %d unique=%d\n",r,us[r]);
  return r;
}
static inline int evald(int r,int c,int b,int*du){
  unsigned char*x=A[r];int gw=0,gc=0;
  for(int q=0;q<M;q++){unsigned char*u=mdc[q];
    if(u[0]!=c&&u[1]!=c&&u[2]!=c&&u[3]!=c)continue;
    int ok=1;
    for(int k=0;k<4;k++){int col=u[k],v=u[4+k];
      if(col==c){if(v!=b){ok=0;break;}}else if(x[col]!=v){ok=0;break;}}
    if(ok){gw+=w[mis[q]];gc++;}
  }
  *du=Lc[r][c]-gc;return L[r][c]-gw;
}
static inline int evalm(int r,int mask,int*sc,int*tv,int*du){
  unsigned char t[K];memcpy(t,A[r],K);
  unsigned cm=0;int n=0,c0=0;
  for(int k=0;k<4;k++)if(mask>>k&1){t[sc[k]]=(unsigned char)tv[k];cm|=1u<<sc[k];n++;c0=sc[k];}
  int loss=0,lc=0;
  if(n==1){loss=L[r][c0];lc=Lc[r][c0];}
  else{int*u=U[r];int m=Uc[r];for(int j=0;j<m;j++){int g=u[j];if(smask[g/NT]&cm){loss+=w[g];lc++;}}}
  int gw=0,gc=0;
  for(int q=0;q<M;q++){unsigned char*h=mdc[q];
    if(t[h[0]]==h[4]&&t[h[1]]==h[5]&&t[h[2]]==h[6]&&t[h[3]]==h[7]){gw+=w[mis[q]];gc++;}}
  *du=lc-gc;return loss-gw;
}
/* Called only with M==0: reassign free cells to drain uniques from weak rows, delete empty rows. */
static int pcell[CAP*K],pbase[PC];
static void polish(double t0,double TL,double budget){
  double ts=now();int stop=0,dels=0;long acc=0,evs=0;
  static int W[NS+2];
  W[0]=1<<20;for(int u=1;u<NS+2;u++)W[u]=(1<<20)/u;
  while(!stop){
    int F=0;
    for(int r=0;r<N;r++)for(int c=0;c<K;c++)if(Lc[r][c]==0)pcell[F++]=r*K+c;
    for(int i=F-1;i>0;i--){int j=rnd()%(i+1);int t=pcell[i];pcell[i]=pcell[j];pcell[j]=t;}
    int improved=0;
    for(int i=0;i<F;i++){
      if((i&15)==0){double tn=now();if(tn-t0>TL||tn-ts>budget){stop=1;break;}}
      int r=pcell[i]/K,c=pcell[i]%K;
      if(Lc[r][c]!=0)continue;
      unsigned char*x=A[r];int a=x[c];E*ec=ce[c];
      long long gn[6]={0,0,0,0,0,0},lossv=0;
      for(int j=0;j<PC;j++){E*e=&ec[j];pbase[j]=e->g+x[e->o0]*e->w0+x[e->o1]*e->w1+x[e->o2]*e->w2;}
      for(int j=0;j<PC;j++){
        if(j+8<PC){int p=pbase[j+8],wc8=ec[j+8].wc;for(int b=0;b<6;b++)__builtin_prefetch(&cnt[p+b*wc8]);}
        int p=pbase[j],wc=ec[j].wc;
        for(int b=0;b<6;b++){
          int idx=p+b*wc;int cv=cnt[idx];
          if(b==a){if(cv==2){int own=xr[idx]^r;lossv+=W[Uc[own]+1];}}
          else if(cv==1){int own=xr[idx];gn[b]+=W[Uc[own]];}
        }
      }
      evs++;
      int bb=-1;long long bg=0;
      for(int b=0;b<6;b++)if(b!=a&&gn[b]-lossv>bg){bg=gn[b]-lossv;bb=b;}
      if(bb>=0){change(r,c,bb);improved=1;acc++;}
    }
    for(int r=N-1;r>=0;r--)if(Uc[r]==0&&N>1){delrow(r);dels++;}
    if(!improved)break;
  }
  long long tu=0;int mn=1<<30,fr=0;
  for(int r=0;r<N;r++){tu+=Uc[r];if(Uc[r]<mn)mn=Uc[r];for(int c=0;c<K;c++)if(Lc[r][c]==0)fr++;}
  fprintf(stderr,"polish N=%d M=%d evals=%ld acc=%ld dels=%d uniq=%lld minU=%d free=%d t=%.0f\n",N,M,evs,acc,dels,tu,mn,fr,now()-t0);
}
/* Augmenting-chain BFS: node = tuple needing cover; edge = move covering it that exposes exactly one unique tuple. */
static int nT[NCAP],nP[NCAP],nR[NCAP],nM[NCAP],nD[NCAP];
static int vis[NG],stampv=0;
static int Jr[JCAP],Jc[JCAP],Jv[JCAP],J=0;
static long long chTry=0,chOk=0,chBad=0,chEx=0,chNodes=0;
static int chain_try(int g0,int cap){
  if(J+4*(DMAX+2)>=JCAP)return 0;
  int M0=M;
  chTry++;
  for(int q=0;q<M;q++){int h=mis[q];int hs=h/NT,hid=h%NT;unsigned char*u=mdc[q];int*hc=scol[hs];
    u[0]=hc[0];u[1]=hc[1];u[2]=hc[2];u[3]=hc[3];u[4]=hid/216;u[5]=(hid/36)%6;u[6]=(hid/6)%6;u[7]=hid%6;}
  stampv++;vis[g0]=stampv;
  int nn=1,head=0;nT[0]=g0;nP[0]=-1;nR[0]=-1;nM[0]=0;nD[0]=0;
  int term=-1,tr=-1,tm=0;
  int off=rnd()%N;
  while(head<nn&&head<cap&&term<0){
    int nd=head++;int g=nT[nd];int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int*sc=scol[g/NT];
    int c0=sc[0],c1=sc[1],c2=sc[2],c3=sc[3];
    chNodes++;
    for(int rr=0;rr<N;rr++){
      int r=rr+off;if(r>=N)r-=N;
      unsigned char*x=A[r];int*lc=Lc[r];
      int mask=0,sum=0;
      if(x[c0]!=tv[0]){sum+=lc[c0];if(sum>CHLIM)continue;mask|=1;}
      if(x[c1]!=tv[1]){sum+=lc[c1];if(sum>CHLIM)continue;mask|=2;}
      if(x[c2]!=tv[2]){sum+=lc[c2];if(sum>CHLIM)continue;mask|=4;}
      if(x[c3]!=tv[3]){sum+=lc[c3];if(sum>CHLIM)continue;mask|=8;}
      if(mask==0)continue;
      int bad=0;
      for(int p=nd;p>=0;p=nP[p])if(nR[p]==r){bad=1;break;}
      if(bad)continue;
      unsigned cm=0;
      unsigned char t[K];memcpy(t,x,K);
      for(int k=0;k<4;k++)if(mask>>k&1){cm|=1u<<sc[k];t[sc[k]]=(unsigned char)tv[k];}
      int loss=0,lg=-1;int*u=U[r];int m=Uc[r];
      for(int j=0;j<m;j++){int g2=u[j];if(smask[g2/NT]&cm){loss++;lg=g2;if(loss>CHLIM)break;}}
      if(loss>CHLIM)continue;
      int ex=0;
      for(int q=0;q<M;q++){unsigned char*h=mdc[q];
        if(t[h[0]]==h[4]&&t[h[1]]==h[5]&&t[h[2]]==h[6]&&t[h[3]]==h[7]&&mis[q]!=g)ex++;}
      if(loss-ex<=0){term=nd;tr=r;tm=mask;break;}
      if(loss==1&&vis[lg]!=stampv&&nn<NCAP&&nD[nd]<DMAX){
        vis[lg]=stampv;nT[nn]=lg;nP[nn]=nd;nR[nn]=r;nM[nn]=mask;nD[nn]=nD[nd]+1;nn++;
      }
    }
  }
  if(term<0){if(head>=nn)chEx++;return 0;}
  int J0=J,cur=term,r=tr,mk=tm;
  while(1){
    int g=nT[cur];int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int*sc=scol[g/NT];
    for(int k=0;k<4;k++)if(mk>>k&1){
      int c=sc[k];
      Jr[J]=r;Jc[J]=c;Jv[J]=A[r][c];J++;
      change(r,c,tv[k]);
    }
    if(nP[cur]<0)break;
    r=nR[cur];mk=nM[cur];cur=nP[cur];
  }
  if(M<M0){chOk++;return 1;}
  while(J>J0){J--;change(Jr[J],Jc[J],Jv[J]);}
  chBad++;
  return 0;
}
/* LNS on a reduced universe: only rows in S move; T = tuples with no coverage outside S, plus holes. */
#define LK 800
#define ZMH 100000
#define ZTC (LK*NS+ZMH)
#define ZJ 200000
#define LNSSEC 16.0
#define WALKSEC 8.0
#define ZEV 64
#define ZMARG 12
#define ZUPH 3
#define ZCOST 6
typedef struct{int t;unsigned char c0,c1,c2,v0,v1,v2;}ZI;
static int *ztid,*zg,*zw,*zhl,*zhp,*zup,*zO;
static unsigned short *zc;
static unsigned *zmask;
static unsigned char (*zm)[8];
static ZI *zIe;
static int (*zU)[NS];
static int zok=0;
static int zS[LK],zn,zh,zUc[LK],zL[LK][K],zLc[LK][K],zoff[K][6],zlen[K][6];
static unsigned char zX[LK][K];
static long long ztb[LK][K];
static int zJr[ZJ],zJc[ZJ],zJv[ZJ];
static double zs[CAP];
static int zord[CAP],zcm[LK],zcl[LK];
static void lns_init(void){
  ztid=malloc(sizeof(int)*(size_t)NG);
  zg=malloc(sizeof(int)*(size_t)ZTC);zw=malloc(sizeof(int)*(size_t)ZTC);
  zhl=malloc(sizeof(int)*(size_t)ZTC);zhp=malloc(sizeof(int)*(size_t)ZTC);
  zup=malloc(sizeof(int)*(size_t)ZTC);zO=malloc(sizeof(int)*(size_t)ZTC);
  zc=malloc(sizeof(unsigned short)*(size_t)ZTC);
  zmask=malloc(sizeof(unsigned)*(size_t)ZTC);
  zm=malloc(8*(size_t)ZTC);
  zIe=malloc(sizeof(ZI)*4*(size_t)ZTC);
  zU=malloc(sizeof(*zU)*(size_t)LK);
  if(!ztid||!zg||!zw||!zhl||!zhp||!zup||!zO||!zc||!zmask||!zm||!zIe||!zU){zok=0;fprintf(stderr,"lns alloc failed\n");return;}
  for(int i=0;i<NG;i++)ztid[i]=-1;
  zok=1;
}
static inline void zhadd(int t){zhp[t]=zh;zhl[zh++]=t;}
static inline void zhdel(int t){int p=zhp[t];int l=zhl[--zh];zhl[p]=l;zhp[l]=p;}
static inline void zuadd(int i,int t){
  zup[t]=zUc[i];zU[i][zUc[i]++]=t;
  unsigned char*m=zm[t];int wt=zw[t];
  for(int q=0;q<4;q++){zL[i][m[q]]+=wt;zLc[i][m[q]]++;}
}
static inline void zudel(int i,int t){
  int p=zup[t];int l=zU[i][--zUc[i]];zU[i][p]=l;zup[l]=p;
  unsigned char*m=zm[t];int wt=zw[t];
  for(int q=0;q<4;q++){zL[i][m[q]]-=wt;zLc[i][m[q]]--;}
}
static void zchange(int i,int c,int b){
  unsigned char*x=zX[i];int a=x[c];if(a==b)return;
  ZI*e=zIe+zoff[c][a];int n=zlen[c][a];
  for(int j=0;j<n;j++,e++){
    if(x[e->c0]!=e->v0||x[e->c1]!=e->v1||x[e->c2]!=e->v2)continue;
    int t=e->t;unsigned short cc=--zc[t];zO[t]^=i;
    if(cc==0){zudel(i,t);zhadd(t);}
    else if(cc==1)zuadd(zO[t],t);
  }
  x[c]=(unsigned char)b;
  e=zIe+zoff[c][b];n=zlen[c][b];
  for(int j=0;j<n;j++,e++){
    if(x[e->c0]!=e->v0||x[e->c1]!=e->v1||x[e->c2]!=e->v2)continue;
    int t=e->t;unsigned short cc=zc[t]++;
    if(cc==0){zhdel(t);zuadd(i,t);}
    else if(cc==1)zudel(zO[t],t);
    zO[t]^=i;
  }
}
static inline int zevalm(int i,int mask,unsigned char*m,int*du){
  unsigned char y[K];memcpy(y,zX[i],K);
  unsigned cm=0;int n=0,c0=0;
  for(int q=0;q<4;q++)if(mask>>q&1){y[m[q]]=m[4+q];cm|=1u<<m[q];n++;c0=m[q];}
  int loss=0,lc=0;
  if(n==1){loss=zL[i][c0];lc=zLc[i][c0];}
  else{int*u=zU[i];int nu=zUc[i];for(int j=0;j<nu;j++){int t=u[j];if(zmask[t]&cm){loss+=zw[t];lc++;}}}
  int gw=0,gc=0;
  for(int q=0;q<zh;q++){int t=zhl[q];unsigned char*h=zm[t];
    if(y[h[0]]==h[4]&&y[h[1]]==h[5]&&y[h[2]]==h[6]&&y[h[3]]==h[7]){gw+=zw[t];gc++;}}
  *du=lc-gc;return loss-gw;
}
static int zcmp(const void*a,const void*b){double x=zs[*(const int*)a],y=zs[*(const int*)b];return x<y?-1:(x>y?1:0);}
static void lns_run(double tend,int k){
  if(!zok||M>ZMH)return;
  if(k>LK)k=LK;
  if(k>N)k=N;
  int M0=M;
  /* choose S: weak rows with noise, plus rows that can cover a hole cheaply */
  for(int r=0;r<N;r++){zs[r]=(Uc[r]+1.0)*(0.6+0.8*(double)(rnd()%1000)/1000.0);zord[r]=r;}
  int hs=M<64?M:64;
  for(int q=0;q<hs;q++){
    int g=(M<=64)?mis[q]:mis[rnd()%M];int id=g%NT;int*sc=scol[g/NT];
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    for(int r=0;r<N;r++){
      unsigned char*x=A[r];int*lc=Lc[r];int cost=0;
      for(int j=0;j<4;j++)if(x[sc[j]]!=tv[j])cost+=lc[sc[j]];
      if(cost<=ZCOST)zs[r]*=0.3;
    }
  }
  qsort(zord,N,sizeof(int),zcmp);
  for(int i=0;i<k;i++){zS[i]=zord[i];memcpy(zX[i],A[zS[i]],K);}
  /* universe T */
  zn=0;
  for(int i=0;i<k;i++){int r=zS[i];for(int s=0;s<NS;s++)cnt[s*NT+tidx(r,s)]--;}
  for(int i=0;i<k;i++){int r=zS[i];for(int s=0;s<NS;s++){int g=s*NT+tidx(r,s);
    if(cnt[g]==0&&ztid[g]<0){ztid[g]=zn;zg[zn++]=g;}}}
  for(int q=0;q<M;q++){int g=mis[q];if(ztid[g]<0){ztid[g]=zn;zg[zn++]=g;}}
  for(int i=0;i<k;i++){int r=zS[i];for(int s=0;s<NS;s++)cnt[s*NT+tidx(r,s)]++;}
  for(int t=0;t<zn;t++){
    int g=zg[t],s=g/NT,id=g%NT;int*sc=scol[s];unsigned char*m=zm[t];
    m[0]=sc[0];m[1]=sc[1];m[2]=sc[2];m[3]=sc[3];
    m[4]=id/216;m[5]=(id/36)%6;m[6]=(id/6)%6;m[7]=id%6;
    zmask[t]=smask[s];zw[t]=1;zc[t]=0;zO[t]=0;
  }
  for(int i=0;i<k;i++){int r=zS[i];for(int s=0;s<NS;s++){int t=ztid[s*NT+tidx(r,s)];
    if(t>=0){zc[t]++;zO[t]^=i;}}}
  memset(zL,0,sizeof(zL[0])*k);memset(zLc,0,sizeof(zLc[0])*k);memset(zUc,0,sizeof(zUc[0])*k);
  memset(ztb,0,sizeof(ztb[0])*k);
  zh=0;
  for(int t=0;t<zn;t++){if(zc[t]==0)zhadd(t);else if(zc[t]==1)zuadd(zO[t],t);}
  memset(zlen,0,sizeof zlen);
  for(int t=0;t<zn;t++){unsigned char*m=zm[t];for(int q=0;q<4;q++)zlen[m[q]][m[4+q]]++;}
  {int off=0;for(int c=0;c<K;c++)for(int v=0;v<6;v++){zoff[c][v]=off;off+=zlen[c][v];zlen[c][v]=0;}}
  for(int t=0;t<zn;t++){unsigned char*m=zm[t];
    for(int q=0;q<4;q++){int c=m[q],v=m[4+q];ZI*e=&zIe[zoff[c][v]+zlen[c][v]++];
      unsigned char oc[3],ov[3];int p=0;
      for(int j=0;j<4;j++)if(j!=q){oc[p]=m[j];ov[p]=m[4+j];p++;}
      e->t=t;e->c0=oc[0];e->c1=oc[1];e->c2=oc[2];e->v0=ov[0];e->v1=ov[1];e->v2=ov[2];}
  }
  /* weighted tabu walk on the reduced problem */
  long long it=0,bumps=0;int best=zh,zj=0;
  while(zh>0){
    if((it&255)==0&&now()>tend)break;
    it++;
    int t=zhl[rnd()%zh];unsigned char*m=zm[t];
    int minub=1<<30;
    for(int i=0;i<k;i++){
      unsigned char*x=zX[i];int*Lr=zL[i];int mask=0,lb=0,ub=0;
      for(int q=0;q<4;q++)if(x[m[q]]!=m[4+q]){int l=Lr[m[q]];ub+=l;if(l>lb)lb=l;mask|=1<<q;}
      zcm[i]=mask;zcl[i]=lb;if(ub<minub)minub=ub;
    }
    int pick=-1,bd=1<<30,ties=0,ev=0;int off=rnd()%k;
    for(int ii=0;ii<k&&ev<ZEV;ii++){
      int i=ii+off;if(i>=k)i-=k;
      if(zcl[i]>minub)continue;
      int mask=zcm[i],tabu=0;
      for(int q=0;q<4;q++)if((mask>>q&1)&&ztb[i][m[q]]>it)tabu=1;
      ev++;
      int du;int d=zevalm(i,mask,m,&du);
      if(du>DUCAP)continue;
      if(tabu&&!(zh+du<best))continue;
      if(d<bd){bd=d;pick=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)pick=i;}
    }
    int apply=0;
    if(pick>=0&&bd<=0)apply=1;
    else{
      for(int q=0;q<zh;q++)zw[zhl[q]]++;
      bumps++;
      if(bumps%2000==0){
        for(int u=0;u<zn;u++){int ow=zw[u];if(ow<=1)continue;
          int nw=1+(ow-1)*4/5;
          if(zc[u]==1&&nw!=ow){int o=zO[u];unsigned char*h=zm[u];for(int q=0;q<4;q++)zL[o][h[q]]+=nw-ow;}
          zw[u]=nw;}
      }
      if(pick>=0&&rnd()%100<ZUPH)apply=1;
    }
    if(apply){
      int mask=zcm[pick];
      /* copy the target: zm[t] is stable, but t may stop being a hole during the move */
      for(int q=0;q<4;q++)if(mask>>q&1){
        int c=m[q];
        zJr[zj]=pick;zJc[zj]=c;zJv[zj]=zX[pick][c];zj++;
        zchange(pick,c,m[4+q]);
        ztb[pick][c]=it+4+(long long)(rnd()%10);
      }
      if(zh<=best){best=zh;zj=0;}
      else if(zh>best+ZMARG||zj+4>=ZJ){while(zj>0){zj--;zchange(zJr[zj],zJc[zj],zJv[zj]);}}
    }
  }
  while(zj>0){zj--;zchange(zJr[zj],zJc[zj],zJv[zj]);}
  /* commit to the global state */
  int nch=0;
  for(int i=0;i<k;i++){int r=zS[i];for(int c=0;c<K;c++)if(zX[i][c]!=A[r][c]){change(r,c,zX[i][c]);nch++;}}
  for(int t=0;t<zn;t++)ztid[zg[t]]=-1;
  fprintf(stderr,"lns k=%d T=%d M %d->%d its=%lld cells=%d\n",k,zn,M0,M,it,nch);
}
static void output(unsigned char (*X)[K],int n){
  for(int r=0;r<n;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",X[r][c]);printf("\n");}
  fflush(stdout);
}
int main(){
  double t0=now(),TL=3500,STALL=300;
  const char*fn=getenv("BEST_CANDIDATE");
  N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;if(N>=CAP)break;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 1;}
  int N0=N; memcpy(I0,A,sizeof(A[0])*N);
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  for(s=0;s<NS;s++)smask[s]=(1u<<scol[s][0])|(1u<<scol[s][1])|(1u<<scol[s][2])|(1u<<scol[s][3]);
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={216,36,6,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){int c=scol[s][i];int j=cc[c]++;E*e=&ce[c][j];e->g=s*NT;e->s=s;e->wc=wts[i];
    int oc[3],ow[3],q=0;for(int k=0;k<4;k++)if(k!=i){oc[q]=scol[s][k];ow[q]=wts[k];q++;}
    e->o0=oc[0];e->o1=oc[1];e->o2=oc[2];e->w0=ow[0];e->w1=ow[1];e->w2=ow[2];}
  for(int i=0;i<NG;i++){w[i]=1;hpos[i]=-1;}
  lns_init();
  rebuild();
  fprintf(stderr,"input N=%d missing=%d t=%.1f\n",N,M,now()-t0);
  int bestN=-1;
  if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
  static int cr[CAP],cmk[CAP],clb[CAP];
  double chT=0;
  long long it=0,bumps=0; int bestMs=M; double tlast=now(),tlog=now(),tlns=0; int restarts=0;
  memset(tb,0,sizeof tb);
  while(1){
    if((it&63)==0){
      double tn=now();
      if(tn-t0>TL)break;
      if(bestN>0 && tn-tlast>STALL && M>0){
        restarts++;
        N=bestN; memcpy(A,B,sizeof(A[0])*N); resetw(); rebuild();
        memset(tb,0,sizeof tb);
        polish(t0,TL,POLSEC);
        bestN=N; memcpy(B,A,sizeof(A[0])*N);
        int br=pickrow((int)(rnd()%(3+restarts)));
        delrow(br); bestMs=M; tlast=now(); J=0;
        fprintf(stderr,"restart %d N=%d missing=%d t=%.0f\n",restarts,N,M,tn-t0);
        continue;
      }
      if(M>0&&tn-tlns>WALKSEC){
        while(J>0){J--;change(Jr[J],Jc[J],Jv[J]);}
        double te=tn+LNSSEC,lim=t0+TL;
        if(te>lim)te=lim;
        static const int ks[3]={150,300,600};
        lns_run(te,ks[rnd()%3]);
        if(M<bestMs)tlast=now();
        bestMs=M;
        tlns=now();
        continue;
      }
      if(tn-tlog>30){tlog=tn;fprintf(stderr,"it=%lld N=%d M=%d bestMs=%d H=%d best=%d t=%.0f chain try=%lld ok=%lld bad=%lld ex=%lld nodes=%lld chT=%.0f\n",it,N,M,bestMs,H,bestN,tn-t0,chTry,chOk,chBad,chEx,chNodes,chT);}
    }
    it++;
    if(M==0){
      if(bestN<0||N<bestN){bestN=N;memcpy(B,A,sizeof(A[0])*N);fprintf(stderr,"valid N=%d t=%.0f it=%lld\n",N,now()-t0,it);}
      polish(t0,TL,POLSEC);
      if(N<bestN)fprintf(stderr,"valid N=%d after polish t=%.0f\n",N,now()-t0);
      bestN=N; memcpy(B,A,sizeof(A[0])*N);
      int br=pickrow(0);
      delrow(br);
      memset(tb,0,sizeof tb);
      decayw(1,2);
      bestMs=M; tlast=now(); J=0;
      continue;
    }
    if((it&127)==0&&M<=CHMAXM){
      double tc=now();
      if(chT<CHFRAC*(tc-t0)){
        for(int q=M;q>0&&M>0;q--)chain_try(mis[rnd()%M],CHCAP);
        chT+=now()-tc;
        if(M<=bestMs){if(M<bestMs)tlast=now();bestMs=M;J=0;}
        continue;
      }
    }
    int g=mis[rnd()%M]; s=g/NT; int id=g%NT;
    int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
    int *sc=scol[s]; int nc=0; int minub=1<<30;
    for(int r=0;r<N;r++){
      unsigned char*x=A[r]; int*Lr=L[r];
      int mask=0,lb=0,ub=0;
      for(int i=0;i<4;i++)if(x[sc[i]]!=tv[i]){int l=Lr[sc[i]];ub+=l;if(l>lb)lb=l;mask|=1<<i;}
      cr[nc]=r;cmk[nc]=mask;clb[nc]=lb;nc++;
      if(ub<minub)minub=ub;
    }
    for(int q=0;q<M;q++){int h=mis[q];int hs=h/NT,hid=h%NT;unsigned char*u=mdc[q];int*hc=scol[hs];
      u[0]=hc[0];u[1]=hc[1];u[2]=hc[2];u[3]=hc[3];u[4]=hid/216;u[5]=(hid/36)%6;u[6]=(hid/6)%6;u[7]=hid%6;}
    int pick=-1,bd=1<<30,ties=0,ev=0;
    int off=rnd()%nc;
    for(int ii=0;ii<nc&&ev<EVMAX;ii++){
      int i=ii+off; if(i>=nc)i-=nc;
      if(clb[i]>minub)continue;
      int r=cr[i],mask=cmk[i];
      int tabu=0;
      for(int k=0;k<4;k++)if((mask>>k&1)&&tb[r][sc[k]]>it)tabu=1;
      ev++;
      int du;int d=evalm(r,mask,sc,tv,&du);
      if(du>DUCAP)continue;
      if(tabu && !(M+du<bestMs))continue;
      if(d<bd){bd=d;pick=i;ties=1;}else if(d==bd){ties++;if(rnd()%ties==0)pick=i;}
    }
    int apply=0;
    if(pick>=0&&bd<=0)apply=1;
    else{
      for(int i=0;i<M;i++){int q=mis[i];if(w[q]==1){hpos[q]=H;hv[H++]=q;}w[q]++;}
      bumps++;
      if(bumps%2000==0)decayw(4,5);
      if(pick>=0&&rnd()%100<UPH)apply=1;
    }
    if(apply){
      int r=cr[pick],mask=cmk[pick];
      for(int k=0;k<4;k++)if(mask>>k&1){
        int c=sc[k];
        if(J<JCAP){Jr[J]=r;Jc[J]=c;Jv[J]=A[r][c];J++;}
        change(r,c,tv[k]);
        tb[r][c]=it+4+(long long)(rnd()%10);
      }
      if(M<=bestMs){if(M<bestMs)tlast=now();bestMs=M;J=0;}
      else if(M>bestMs+MARG||J+4>=JCAP){while(J>0){J--;change(Jr[J],Jc[J],Jv[J]);}}
    }
  }
  if(bestN<0){
    if(M==0){bestN=N;memcpy(B,A,sizeof(A[0])*N);}
    else {bestN=N0;memcpy(B,I0,sizeof(A[0])*N0);}
  }
  /* full verification */
  N=bestN; memcpy(A,B,sizeof(A[0])*N); resetw(); rebuild();
  fprintf(stderr,"final N=%d missing=%d iters=%lld\n",bestN,M,it);
  if(M==0) output(B,bestN);
  else { fprintf(stderr,"verification failed, printing input\n"); output(I0,N0); }
  return 0;
}
