#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#define K 11
#define V 6
#define NS 330
#define NT 1296
#define NG (NS*NT)
#define CAP 5000
#define HOLE 1000
#define TAILS 10
#define TAILT 25.0
static int N;
static unsigned char A[CAP][K],B[CAP][K],I0[CAP][K];
static unsigned short cnt[NG];
static int mis[NG],pos[NG],M=0;
static long long tb[CAP][K];
static int scol[NS][4];
static int csg[K][165],cot[K][165][3],cw[K][165],cwo[K][165][3];
static int remc[CAP+2],addc[CAP+2];
static unsigned long long rs=0xD1B54A32D192ED03ULL;
static volatile sig_atomic_t stopf=0;
static long long fit=0;
static unsigned short xr[NG];
static short brk[CAP][K];
static int csi[K][165];
#define BW (CAP/64+1)
static unsigned long long bset[K][V][BW];
static double expt[4096];
static double slk=0;
static double T=60;
static void onterm(int x){(void)x;stopf=1;}
static inline unsigned long long rnd(void){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static inline void addm(int g){pos[g]=M;mis[M++]=g;}
static inline void delm(int g){int p=pos[g];int l=mis[--M];mis[p]=l;pos[l]=p;pos[g]=-1;}
static inline int tidx(unsigned char*x,int s){return x[scol[s][0]]*216+x[scol[s][1]]*36+x[scol[s][2]]*6+x[scol[s][3]];}
static inline void addb(int r,int s,int d){
  short*b=brk[r];int*sc=scol[s];
  b[sc[0]]+=d;b[sc[1]]+=d;b[sc[2]]+=d;b[sc[3]]+=d;
}
static void setsoft(int s){
  memset(remc,0,sizeof remc);memset(addc,0,sizeof addc);
  slk=5.0*s;
  remc[1]=addc[0]=HOLE;
  remc[2]=addc[1]=s;
  remc[3]=addc[2]=s/3;
  remc[4]=addc[3]=s/9;
}
static inline int delta(int r,int c,int a,int b){
  unsigned char*x=A[r];int d=0;
  for(int j=0;j<165;j++){
    int g=csg[c][j]+x[cot[c][j][0]]*cwo[c][j][0]+x[cot[c][j][1]]*cwo[c][j][1]+x[cot[c][j][2]]*cwo[c][j][2];
    d+=remc[cnt[g+a*cw[c][j]]]-addc[cnt[g+b*cw[c][j]]];
  }
  return d;
}
static void change(int r,int c,int b){
  int a=A[r][c]; if(a==b)return;
  unsigned char*x=A[r];
  for(int j=0;j<165;j++){
    int g=csg[c][j]+x[cot[c][j][0]]*cwo[c][j][0]+x[cot[c][j][1]]*cwo[c][j][1]+x[cot[c][j][2]]*cwo[c][j][2];
    int o=g+a*cw[c][j],n=g+b*cw[c][j],s=csi[c][j];
    xr[o]^=(unsigned short)r;
    int q=--cnt[o];
    if(q==0){addm(o);addb(r,s,-1);}
    else if(q==1)addb(xr[o],s,1);
    q=cnt[n]++;
    if(q==0){delm(n);addb(r,s,1);}
    else if(q==1)addb(xr[n],s,-1);
    xr[n]^=(unsigned short)r;
  }
  bset[c][a][r>>6]&=~(1ULL<<(r&63));
  bset[c][b][r>>6]|=1ULL<<(r&63);
  A[r][c]=(unsigned char)b;
}
static void rebuild(void){
  memset(cnt,0,sizeof cnt);memset(xr,0,sizeof xr);memset(brk,0,sizeof brk);memset(bset,0,sizeof bset);M=0;
  for(int r=0;r<N;r++){
    for(int s=0;s<NS;s++){int g=s*NT+tidx(A[r],s);cnt[g]++;xr[g]^=(unsigned short)r;}
    for(int c=0;c<K;c++)bset[c][A[r][c]][r>>6]|=1ULL<<(r&63);
  }
  for(int g=0;g<NG;g++){
    pos[g]=-1;
    if(cnt[g]==0)addm(g);
    else if(cnt[g]==1)addb(xr[g],g/NT,1);
  }
}
static void delweak(void){
  int br=-1,bu=1<<30,ties=0;
  for(int r=0;r<N;r++){
    int u=0;
    for(int s=0;s<NS;s++)if(cnt[s*NT+tidx(A[r],s)]==1)u++;
    if(u<bu){bu=u;br=r;ties=1;}
    else if(u==bu){ties++;if(rnd()%ties==0)br=r;}
  }
  N--;
  if(br!=N){memcpy(A[br],A[N],K);memcpy(tb[br],tb[N],sizeof tb[0]);}
  rebuild();
  fprintf(stderr,"deleted row uq=%d, N=%d M=%d\n",bu,N,M);
}
/* all distance-1 rows of a random hole, noisy Boltzmann choice, short tabu */
static void focus(void){
  int g=mis[rnd()%M],s=g/NT,id=g%NT;
  int tv[4]={id/216,(id/36)%6,(id/6)%6,id%6};
  int*sc=scol[s];
  int br=-1,bc=-1,bv=-1;
  double bs=0;
  fit++;
  unsigned long long*p0=bset[sc[0]][tv[0]],*p1=bset[sc[1]][tv[1]],*p2=bset[sc[2]][tv[2]],*p3=bset[sc[3]][tv[3]];
  int W=(N+63)/64;
  for(int w=0;w<W;w++){
    unsigned long long a0=p0[w],a1=p1[w],a2=p2[w],a3=p3[w];
    unsigned long long m=((a0&a1)&(a2^a3))|((a2&a3)&(a0^a1));
    while(m){
      int r=w*64+__builtin_ctzll(m);
      m&=m-1;
      unsigned char*x=A[r];
      int mc=(x[sc[0]]!=tv[0])?0:(x[sc[1]]!=tv[1])?1:(x[sc[2]]!=tv[2])?2:3;
      int c=sc[mc];
      double th=T*expt[rnd()>>52];
      if((brk[r][c]-1)*HOLE-slk-th>=bs)continue;
      int d=delta(r,c,x[c],tv[mc]);
      if(tb[r][c]>fit&&d>-HOLE/2)continue;
      double q=d-th;
      if(q<bs){bs=q;br=r;bc=c;bv=tv[mc];}
    }
  }
  if(br>=0){
    change(br,bc,bv);
    tb[br][bc]=fit+10+(long long)(rnd()%20);
  }
}
int main(void){
  double t0=now(),TL=3480;
  signal(SIGTERM,onterm);
  const char*fn=getenv("BEST_CANDIDATE");
  N=0;
  if(fn){FILE*f=fopen(fn,"r");if(f){int v,c=0;while(fscanf(f,"%d",&v)==1){A[N][c++]=(unsigned char)v;if(c==K){c=0;N++;if(N>=CAP)break;}}fclose(f);}}
  if(N==0){fprintf(stderr,"no input\n");return 1;}
  int N0=N;memcpy(I0,A,(size_t)N*K);
  int s=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int c=b+1;c<K;c++)for(int d=c+1;d<K;d++){scol[s][0]=a;scol[s][1]=b;scol[s][2]=c;scol[s][3]=d;s++;}
  int cc[K];memset(cc,0,sizeof cc);
  int wts[4]={216,36,6,1};
  for(s=0;s<NS;s++)for(int i=0;i<4;i++){
    int c=scol[s][i],j=cc[c]++,q=0;
    csg[c][j]=s*NT;cw[c][j]=wts[i];csi[c][j]=s;
    for(int k=0;k<4;k++)if(k!=i){cot[c][j][q]=scol[s][k];cwo[c][j][q]=wts[k];q++;}
  }
  rebuild();
  for(int i=0;i<4096;i++)expt[i]=-log((i+0.5)/4096.0);
  fprintf(stderr,"input N=%d missing=%d\n",N,M);
  int bestN=-1;
  static const int CM[6]={100,50,200,100,50,25},CS[6]={30,30,30,15,60,30};
  int cycle=-1,phase=0,Mt=100,focE=16,bestMs=M,succ=0;
  double Th=450,tn=now(),phStart=tn,tlog=tn,tprev=tn,D=40,Tc0=450,Tc1=290,Dc=40;
  long long it=0,accU=0;
  setsoft(TAILS);T=TAILT;
  while(1){
    if((it&4095)==0){
      tn=now();
      if(tn-t0>TL||stopf)break;
      double el=tn-phStart,dt=tn-tprev;
      tprev=tn;
      if(tn-tlog>60){
        tlog=tn;
        long long h[6]={0,0,0,0,0,0};
        for(int g=0;g<NG;g++)h[cnt[g]<5?cnt[g]:5]++;
        fprintf(stderr,"t=%.0f cyc=%d ph=%d N=%d M=%d bestMs=%d T=%.1f Th=%.1f D=%.0f Mt=%d h0=%lld h1=%lld h2=%lld h3=%lld h4=%lld accU=%lld it=%lld\n",
          tn-t0,cycle,phase,N,M,bestMs,T,Th,D,Mt,h[0],h[1],h[2],h[3],h[4],accU,it);
      }
      if(phase==0){
        double dur=(M<=2)?30:8;
        if(el>dur){
          cycle++;
          Mt=CM[cycle%6];
          setsoft(CS[cycle%6]);
          T=Th;focE=128;phase=1;phStart=tn;
          fprintf(stderr,"cycle %d start t=%.0f N=%d M=%d bestMs=%d Mt=%d s=%d D=%.0f\n",cycle,tn-t0,N,M,bestMs,Mt,CS[cycle%6],D);
        }
      } else if(phase==1){
        T*=exp((M>Mt?-0.3:0.3)*dt);
        if(T<10)T=10;
        if(T>3000)T=3000;
        Th=T;
        if(el>(D/2>10?D/2:10)){
          phase=2;phStart=tn;focE=64;succ=0;
          Tc0=Th;Tc1=0.65*Th;Dc=D;
        }
      } else {
        double f=el/Dc;
        if(f>=1){
          if(!succ){D*=1.2;if(D>900)D=900;}
          phase=0;phStart=tn;focE=16;
          setsoft(TAILS);T=TAILT;
        } else T=Tc0*pow(Tc1/Tc0,f);
      }
    }
    it++;
    if(M==0){
      if(bestN<0||N<bestN){
        bestN=N;memcpy(B,A,(size_t)N*K);
        fprintf(stderr,"valid N=%d t=%.0f\n",N,now()-t0);
      }
      delweak();
      bestMs=M;
      tn=now();
      if(phase==2){
        /* chain: keep descending on the next row instead of restarting the cycle */
        succ=1;
        D*=0.9;if(D<15)D=15;
        Tc0=T*1.02;if(Tc0>Th)Tc0=Th;
        if(Tc0<Tc1*1.05)Tc0=Tc1*1.05;
        Dc=D*log(Tc0/Tc1)/log(1/0.65);if(Dc<5)Dc=5;
        phStart=tn;
      } else if(phase==0)phStart=tn;
      continue;
    }
    if(M<bestMs){
      bestMs=M;
      if(M<=3)fprintf(stderr,"N=%d bestMs=%d t=%.0f\n",N,M,now()-t0);
    }
    if((it&(focE-1))==0){focus();continue;}
    {
      unsigned long long u=rnd();
      int r=(int)(((u&0xFFFFFULL)*(unsigned long long)N)>>20);
      int c=(int)((((u>>20)&0xFFFFULL)*K)>>16);
      double th=T*expt[u>>52];
      if(brk[r][c]*HOLE>th+slk)continue;
      int a=A[r][c],b=(a+1+(int)((((u>>36)&0xFFFFULL)*5)>>16))%V;
      int d=delta(r,c,a,b);
      if((double)d<=th){change(r,c,b);accU++;}
    }
  }
  if(M==0&&(bestN<0||N<bestN)){bestN=N;memcpy(B,A,(size_t)N*K);}
  /* independent verification of the array to be printed */
  int ok=bestN>0;
  if(ok){
    static unsigned char seen[NG];
    memset(seen,0,sizeof seen);
    for(int r=0;r<bestN;r++)for(s=0;s<NS;s++)seen[s*NT+tidx(B[r],s)]=1;
    for(int i=0;i<NG;i++)if(!seen[i]){ok=0;break;}
  }
  if(!ok){bestN=N0;memcpy(B,I0,(size_t)N0*K);fprintf(stderr,"falling back to input\n");}
  fprintf(stderr,"final N=%d iters=%lld\n",bestN,it);
  for(int r=0;r<bestN;r++){for(int c=0;c<K;c++)printf(c?" %d":"%d",B[r][c]);printf("\n");}
  return 0;
}
