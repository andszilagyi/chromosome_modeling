#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "rand.h"

//#define FEATURE_NO_MUT
//#define FEATURE_NO_ASSORT

#define COU (10000000)
#define COMP_NUM (5000)
#define LR (20)
#define LM (80)
#define L (LM+LR)
#define EIGHT (8)
#define alpha (5)
#define beta (15)
#define RESTART (50)

#define MC_MAX (100)
#define SD_MAX (100)
#define D_MAX (10)

int SD;
int D;
int MC;
double mu;

double nu1 = 0.01;
double nu3 = 0.01;
double nu2 = 0.01;

//double replaTbl[D_MAX] = {1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8};
//double replaTbl[D_MAX] = {1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
double replaTbl[D_MAX] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
//double replaTbl[D_MAX] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9};

typedef struct _SEQ
{
  int EM[LM];
  int ER[LR];
  int type;
  int next;
  int head;
}SEQ;

typedef struct _COMP
{
  SEQ SEQ[SD_MAX+MC_MAX];
  int numtot;
  double fitn;
  int timer;
}COMP;

typedef struct _SORT
{
  int seq;
  int len;
}SORT;

COMP comp[COMP_NUM];

int METACT[D_MAX][LM];
int REPLACT[LR];
int num[MC_MAX];
double EA[LM];

void error(int ev, const char* errs)
{
  printf("%s\n", errs);
  exit(ev);
}

int sortfunction(const void *a, const void *b)
{
  return((((SORT*)b)->len) - (((SORT*)a)->len));
}

double act2(int *v1, int *v2)
{
  int hd,i;

  hd = 0;
  for(i = 0; i < LM; i++)
  {
    if(v1[i] != v2[i])
      hd++;
  }
  return(EA[hd]);
}

double fitness(int t)
{
  int i,j;
  double fit=0.0, Atemp[D];

  for(i=0;i<D;i++)
    Atemp[i]=0.0;

  for(j=0;j<comp[t].numtot;j++)
    Atemp[comp[t].SEQ[j].type]+=act2(METACT[comp[t].SEQ[j].type],comp[t].SEQ[j].EM);

  for(j=0;j<D;j++)
    if(Atemp[j]==0)
      return(0.0);

  for(j=0;j<D;j++)
    fit+=1.0/Atemp[j];

  fit=D*D/(SD*fit);

  return(fit);
}

void clear_vesicule(int t)
{
  int j,k;

  comp[t].fitn=0.0;
  comp[t].numtot=0;

  for(j=0;j<SD+MC;j++)
  {
    comp[t].SEQ[j].next=-1;
    comp[t].SEQ[j].head=0;
    for(k=0;k<LM;k++)
      comp[t].SEQ[j].EM[k]=0;
    for(k=0;k<LR;k++)
      comp[t].SEQ[j].ER[k]=0;
  }
}

void init(void)
{
  int i,j,m,t;

  for(t=0;t<LM;t++)
  {
    EA[t]=(1.0/exp(pow((t),2))-1.0/exp(pow(LM,2)))/(1.0-1.0/exp(pow(LM,2)));
    if(EA[t]<1E-4)
      EA[t]=0.0;
  }

  // nullazas
  for(t=0;t<COMP_NUM;t++)
    clear_vesicule(t);


  for(i=0;i<D;i++)
    for(j=0;j<LM;j++)
      METACT[i][j]=0;


  for(j=0;j<LR;j++)
    REPLACT[j]=0;

  // defining arbitrary pattern for replicase activity
  for(j=0;j<LR;j++)
  {
    if(j%2)
      REPLACT[j]=0;
    else
      REPLACT[j]=1;
  }

  if(EIGHT*D>LM)
    error(17, "#ERROR -- 17");

  for(i=0;i<D;i++)
    for(j=0;j<EIGHT;j++)
      METACT[i][EIGHT*i+j]=1;



  for(t=0;t<COMP_NUM;t++)
  {
    comp[t].timer = 0;
    int filled=randl(SD-1)+1;
    for(i=0;i<filled;i++)
    {
      comp[t].numtot++;
      comp[t].SEQ[i].type=i%D;
      comp[t].SEQ[i].head=1;
      comp[t].SEQ[i].next=-1;
      for(m=0;m<LR;m++)
        comp[t].SEQ[i].ER[m]=REPLACT[m];
      for(m=0;m<LM;m++)
        comp[t].SEQ[i].EM[m]=METACT[i%D][m];
    }
  }
  for(t=0;t<COMP_NUM;t++)
  {
    comp[t].fitn=fitness(t);
  }
}

int choose_according_fitnes(void)
{
  int ch=-99,t;
  double sumfit=0.0,ran,sum=0.0,sumold=0.0;

  for(t=0;t<COMP_NUM;t++)
    sumfit+=comp[t].fitn;

  ran=sumfit*randd();

  for(t=0;t<COMP_NUM;t++)
  {
    sum+=comp[t].fitn;
    if((sumold<=ran) && (sum>ran))
      break;
    sumold=sum;
  }

  ch=t;
  if(ch>=COMP_NUM)
    return(-1);
  return(ch);
}

int replicate_random_template(int t)
{
  int j,k,hd,c,replacnt=0,num;
  double ran,repla[SD_MAX+MC_MAX],sumrepla=0.0,sum=0.0;
  int rth, rtns, rtnt, rtn;
  int c1, c2, cl;
  int typ=-1;

  if(comp[t].timer > 0)
  {
    comp[t].timer--;
    return(0);
  }

  for(j=0;j<comp[t].numtot;j++)
  {
    if(comp[t].SEQ[j].head==0)
      continue;

    num=0;
    for(c=j;c!=-1;c=comp[t].SEQ[c].next)
      num++;

    hd=0;
    for(k=0;k<LR;k++)
      if(comp[t].SEQ[j].ER[k]!=REPLACT[k])
        hd++;

    if(replacnt>(SD+MC))
      error(4, "ERROR -- 4");

    repla[replacnt] = (1.0 - (pow((double)hd, alpha)) / (beta + pow((double)hd, alpha)));
    repla[replacnt] *= replaTbl[comp[t].SEQ[j].type];
    sumrepla+=repla[replacnt];
    replacnt++;
  }

  ran=sumrepla*randd();

  replacnt=0;
  for(j=0;j<comp[t].numtot;j++)
  {
    if(comp[t].SEQ[j].head==0)
      continue;
    sum+=repla[replacnt];
    replacnt++;
    if(sum > ran)
      break;
  }

   if(j == comp[t].numtot) // tulolvasast jelentene
    return(0);

  if(comp[t].SEQ[j].head!=1)
    error(18, "HIBA -- 18");

  rth = -1;
  rtns = -1;
  rtnt = -1;
  rtn = -1;
  if(randd()<nu3)
  {
    int ctmp[SD_MAX+MC_MAX];
    int ctmpcnt = 0;
    for(c=j,c1=0;c!=-1;c=comp[t].SEQ[c].next)
      c1++;
    rtnt = randl(c1);
    for(c=j,cl=0;(c!=-1) && (cl<=rtnt);c=comp[t].SEQ[c].next)
      if(cl == rtnt) typ = comp[t].SEQ[c].type;
    for(c=0;c<comp[t].numtot;c++)
    {
      if((comp[t].SEQ[c].head == 1) && (c!=j))
        ctmp[ctmpcnt++] = c;
    }
    if(ctmpcnt != 0)
    {
      rth = ctmp[randl(ctmpcnt)];
      for(c=rth,c2=0,ctmpcnt=0;c!=-1;c=comp[t].SEQ[c].next)
      {
        if(comp[t].SEQ[c].type == typ)
          ctmp[ctmpcnt++] = c2;
        c2++;
      }
      if(ctmpcnt > 1)
      {
        rtn=ctmp[randl(ctmpcnt)];
        for(c=rth,cl=0;(c!=-1) && (cl<=rtn);c=comp[t].SEQ[c].next)
        {
          rtns = c;
          cl++;
        }
      }else
      {
        rtnt = -1;
      }
    }else
    {
      rtnt = -1;
    }
  }
  if(rtnt != -1)
  {
    if(rtns == -1)
      error(99, "something is not good!");
    if(((rtnt + 1) + (c2 - rtn)) > MC) // tul hosszu lenne
    {
      printf("# too long...\n");
      rtnt = -1;
    }
  }

  for(c=j,cl=0;c!=-1;cl++)
  {
    comp[t].SEQ[comp[t].numtot].type=comp[t].SEQ[c].type;

    if(comp[t].SEQ[c].next==-1)
      comp[t].SEQ[comp[t].numtot].next=-1;
    else
      comp[t].SEQ[comp[t].numtot].next=comp[t].numtot+1;
    if(cl == 0)
      comp[t].SEQ[comp[t].numtot].head = 1;
    else
      comp[t].SEQ[comp[t].numtot].head = 0;
    for(k=0;k<LR;k++)
    {
#ifndef FEATURE_NO_MUT
      if(randd()<mu)
        comp[t].SEQ[comp[t].numtot].ER[k]=-1*comp[t].SEQ[c].ER[k]+1;
      else
        comp[t].SEQ[comp[t].numtot].ER[k]=comp[t].SEQ[c].ER[k];
#else
    comp[t].SEQ[comp[t].numtot].ER[k]=comp[t].SEQ[c].ER[k];
#endif
    }
    for(k=0;k<LM;k++)
    {
      if(randd()<mu)
        comp[t].SEQ[comp[t].numtot].EM[k]=-1*comp[t].SEQ[c].EM[k]+1;
      else
        comp[t].SEQ[comp[t].numtot].EM[k]=comp[t].SEQ[c].EM[k];
    }
    comp[t].numtot++;

    if(cl == rtnt) c=rtns; else c=comp[t].SEQ[c].next;
  }
  if(cl<1)
      error(1, "1");
  comp[t].timer = cl-1;
  comp[t].fitn=fitness(t);
  return(0);
}

int chromosomatization(int t)
{
  int h1,h2,t1=-1,c1,c2,i;
  int ctmp[SD_MAX+MC_MAX];
  int ctmpcnt = 0;

  for(i=0;i<comp[t].numtot;i++)
    if(comp[t].SEQ[i].head==1)
      ctmp[ctmpcnt++] = i;
  if(ctmpcnt < 2) return(-1);
  h1=ctmp[randl(ctmpcnt)];
  h2=ctmp[randl(ctmpcnt)];
  while(h1==h2) h2=ctmp[randl(ctmpcnt)];
  for(i=h1,c1=0;i!=-1;i=comp[t].SEQ[i].next)
  {
    c1++;
    if(comp[t].SEQ[i].next==-1)
      t1=i;
  }
  for(i=h2,c2=0;i!=-1;i=comp[t].SEQ[i].next)
    c2++;
  if((c1+c2)>MC)
    return(-1);
  comp[t].SEQ[t1].next=h2;
  comp[t].SEQ[h2].head=0;
  return(0);
}

int chrbreak(int t)
{
  int ctmp[SD_MAX+MC_MAX];
  int ctmpcnt = 0, b, c, cc, i, ii, j, cl;
  for(j=0;j<comp[t].numtot;j++)
  {
    if((comp[t].SEQ[j].head == 1) && (comp[t].SEQ[j].next != -1))
      ctmp[ctmpcnt++] = j;
  }
  if(ctmpcnt == 0) return(0);
  c=randl(ctmpcnt);
  cl=0;
  for(i=ctmp[c];i!=-1;i=comp[t].SEQ[i].next)
    cl++;
  if(cl < 2) return(-2);
  b=randl(cl-1);
  for(i=ctmp[c], cc=0;(i!=-1) && (cc < b); i=comp[t].SEQ[i].next, cc++){;}
  ii = comp[t].SEQ[i].next;
  if(ii == -1) return(-3);
  comp[t].SEQ[ii].head = 1;
  comp[t].SEQ[i].next = -1;
  return(0);
}

int pn = 0;

void division(int t)
{
    int c,j,k,u;
#ifdef FEATURE_NO_ASSORT
    int cc = randl(2);
    int gc[D];
    SORT chrs[SD_MAX+MC_MAX];
    int chrscnt = 0;
    for(int t = 0; t < D; t++)
        gc[t] = randl(2);
#endif
    COMP orig;
    memcpy(&orig,&(comp[t]),sizeof(COMP));

    u=randl(COMP_NUM);
    while(u==t) u=randl(COMP_NUM);

    clear_vesicule(t);
    clear_vesicule(u);
    for(j=0;j<orig.numtot;j++)
    {
      if(orig.SEQ[j].head != 1) continue;
#ifndef FEATURE_NO_ASSORT
      if(randd()<0.5)
      {
        for(c=j;c!=-1;c=orig.SEQ[c].next)
        {
          comp[t].SEQ[comp[t].numtot].type=orig.SEQ[c].type;
          if(orig.SEQ[c].next==-1)
            comp[t].SEQ[comp[t].numtot].next=-1;
          else
            comp[t].SEQ[comp[t].numtot].next=comp[t].numtot+1;
          comp[t].SEQ[comp[t].numtot].head=orig.SEQ[c].head;
          for(k=0;k<LR;k++)
            comp[t].SEQ[comp[t].numtot].ER[k]=orig.SEQ[c].ER[k];
          for(k=0;k<LM;k++)
            comp[t].SEQ[comp[t].numtot].EM[k]=orig.SEQ[c].EM[k];
          comp[t].numtot++;
        }
      }
      else
      {
        for(c=j;c!=-1;c=orig.SEQ[c].next)
        {
          comp[u].SEQ[comp[u].numtot].type=orig.SEQ[c].type;
          if(orig.SEQ[c].next==-1)
            comp[u].SEQ[comp[u].numtot].next=-1;
          else
            comp[u].SEQ[comp[u].numtot].next=comp[u].numtot+1;
          comp[u].SEQ[comp[u].numtot].head=orig.SEQ[c].head;
          for(k=0;k<LR;k++)
            comp[u].SEQ[comp[u].numtot].ER[k]=orig.SEQ[c].ER[k];
          for(k=0;k<LM;k++)
            comp[u].SEQ[comp[u].numtot].EM[k]=orig.SEQ[c].EM[k];
          comp[u].numtot++;
        }
      }
#else
      if(orig.SEQ[j].next == -1)
      {
        int cs = 0;
        cs = ((gc[comp[t].SEQ[j].type]++) & 1);
        if(cs == 1)
        {
          comp[t].SEQ[comp[t].numtot].type=orig.SEQ[j].type;
          comp[t].SEQ[comp[t].numtot].next=-1;
          comp[t].SEQ[comp[t].numtot].head=orig.SEQ[j].head;
          for(k=0;k<LR;k++)
            comp[t].SEQ[comp[t].numtot].ER[k]=orig.SEQ[j].ER[k];
          for(k=0;k<LM;k++)
            comp[t].SEQ[comp[t].numtot].EM[k]=orig.SEQ[j].EM[k];
          comp[t].numtot++;
        }
        else
        {
          comp[u].SEQ[comp[u].numtot].type=orig.SEQ[j].type;
          comp[u].SEQ[comp[u].numtot].next=-1;
          comp[u].SEQ[comp[u].numtot].head=orig.SEQ[j].head;
          for(k=0;k<LR;k++)
            comp[u].SEQ[comp[u].numtot].ER[k]=orig.SEQ[j].ER[k];
          for(k=0;k<LM;k++)
            comp[u].SEQ[comp[u].numtot].EM[k]=orig.SEQ[j].EM[k];
          comp[u].numtot++;
        }
      }
      else
      {
        int c = 0;
        for(int i = j; i != -1; i = orig.SEQ[i].next) c++;
        chrs[chrscnt].seq = j;
        chrs[chrscnt].len = c;
        chrscnt++;
      }
#endif
    }
#ifdef FEATURE_NO_ASSORT
    if(chrscnt > 0)
    {
      for(int j = 0; j < chrscnt; j++)
      {
        if((cc & 1) == 0)
        {
          for(c=chrs[j].seq; c!=-1; c=orig.SEQ[c].next)
          {
            comp[t].SEQ[comp[t].numtot].type=orig.SEQ[c].type;
            if(orig.SEQ[c].next==-1)
              comp[t].SEQ[comp[t].numtot].next=-1;
            else
              comp[t].SEQ[comp[t].numtot].next=comp[t].numtot+1;
            comp[t].SEQ[comp[t].numtot].head=orig.SEQ[c].head;
            for(k=0;k<LR;k++)
              comp[t].SEQ[comp[t].numtot].ER[k]=orig.SEQ[c].ER[k];
            for(k=0;k<LM;k++)
              comp[t].SEQ[comp[t].numtot].EM[k]=orig.SEQ[c].EM[k];
            comp[t].numtot++;
          }
        }
        else
        {
          for(c=chrs[j].seq; c!=-1; c=orig.SEQ[c].next)
          {
            comp[u].SEQ[comp[u].numtot].type=orig.SEQ[c].type;
            if(orig.SEQ[c].next==-1)
              comp[u].SEQ[comp[u].numtot].next=-1;
            else
              comp[u].SEQ[comp[u].numtot].next=comp[u].numtot+1;
            comp[u].SEQ[comp[u].numtot].head=orig.SEQ[c].head;
            for(k=0;k<LR;k++)
              comp[u].SEQ[comp[u].numtot].ER[k]=orig.SEQ[c].ER[k];
            for(k=0;k<LM;k++)
              comp[u].SEQ[comp[u].numtot].EM[k]=orig.SEQ[c].EM[k];
            comp[u].numtot++;
          }
        }
        cc++;
    }
  }
#endif
  comp[u].fitn=fitness(u);
  comp[t].fitn=fitness(t);

  if(comp[t].numtot >= SD)
  {
    clear_vesicule(t);
  }
  if(comp[u].numtot >= SD)
  {
    clear_vesicule(u);
  }
}

double aver_fitn(void)
{
  int t;
  double avf=0.0;
  int lc = 0;

  for(t=0;t<COMP_NUM;t++)
  {
    avf+=comp[t].fitn;
    if(comp[t].fitn > 0.0) lc++;
  }
  if(lc==0)
    return(0.0);
  avf=avf/(double)COMP_NUM;

  return(avf);
}

double chrom_div(void)
{
  int a[D];
  int t, i,j;
  double h=0.0,temp;

  for(j=0;j<COMP_NUM;j++)
  {
    for(t=0;t<D;t++)
      a[t]=0;
    for(i=0;i<comp[j].numtot;i++)
      a[comp[j].SEQ[i].type]++;
    temp=0.0;
    for(t=0;t<D;t++)
    {
      if(a[t]>0)
        temp+=-1.0*(double)a[t]/(double)comp[j].numtot*log((double)a[t]/(double)comp[j].numtot);
    }
    h+=temp;
  }

  return(h/log((double)D)/(double)COMP_NUM);
}

long int space_usage(void)
{
  int t;
  long int sum = 0;

  for(t=0;t<COMP_NUM;t++)
    sum += comp[t].numtot;
  return(sum);
}

void chr_num(void)
{
  int i,j,t,c;

  memset(num,0,sizeof(num));

  for(t=0;t<COMP_NUM;t++)
  {
    for(j=0;j<comp[t].numtot;j++)
    {
      if(comp[t].SEQ[j].head==1)
      {
        c=0;
        for(i=j;i!=-1;i=comp[t].SEQ[i].next)
          c++;
        if(c!=0)
           num[c-1]++;
      }
    }
  }
}

double avg_chr_num(void)
{
  double a = 0;
  double b = 0;
  int i;

  chr_num();
  for(i = 0; i < MC; i++)
  {
    a += num[i] * (i+1.0);
    b += num[i];
  }
  return(a/((double)b));
}

double non_zero_fitness(void)
{
  int t,nu=0;

  for(t=0;t<COMP_NUM;t++)
    if(fitness(t)>0.01)
      nu++;

  return((double)nu/(double)COMP_NUM);
}

void count_chr(double result[MC_MAX], int result2[MC_MAX][D_MAX])
{
  int c, g, d, t, co, gc;
  int tmp[D_MAX];
  int gcsum[MC_MAX];
  int gccnt[MC_MAX];

  for(t = 0; t < MC; t++)
  {
    result[t] = 0.0;
    gcsum[t] = 0;
    gccnt[t] = 0;
    for(d = 0; d < D; d++)
      result2[t][d] = 0;
  }

  for(c = 0; c < COMP_NUM; c++)
  {
    for(g = 0; g < comp[c].numtot; g++)
    {
      if(comp[c].SEQ[g].head == 1)
      {
        for(t = 0; t < D; t++)
          tmp[t] = 0;
        for(t = g, co = 0; t != -1; t = comp[c].SEQ[t].next, co++)
        {
          if(act2(METACT[comp[c].SEQ[t].type],comp[c].SEQ[t].EM)>0.0) tmp[comp[c].SEQ[t].type]++;
        }
        for(t = 0, gc = 0; t < D; t++)
        {
          if(tmp[t] != 0)
            gc++;
          result2[co-1][t] += tmp[t];
        }
        gcsum[co-1] += gc;
        gccnt[co-1]++;
      }
    }
  }
  for(t = 0; t < MC; t++)
    if(gccnt[t] != 0)
      result[t] = ((double)gcsum[t]) / ((double)gccnt[t]);
}

void mutatedness(double *replication, double *metabolic)
{
  int nrc = 0, nmc = 0, nr = 0, nm = 0;
  int t, i, m;

  *replication = 0.0;
  *metabolic = 0.0;
  for(t=0;t<COMP_NUM;t++)
  {
    for(i=0;i<comp[t].numtot;i++)
    {
      if(comp[t].SEQ[i].head == 1)
      {
        for(m=0;m<LR;m++)
          nr += comp[t].SEQ[i].ER[m] ^ REPLACT[m];
        nrc += LR;
      }
      for(m=0;m<LM;m++)
        nm += (comp[t].SEQ[i].EM[m] ^ METACT[comp[t].SEQ[i].type][m]);
      nmc += LM;
    }
  }
  if(nrc != 0) *replication = (double)nr/(double)nrc;
  if(nmc != 0) *metabolic = (double)nm/(double)nmc;
}

int main(int argc, char* argv[])
{

  int t,totnum,cou;
  double res[MC_MAX];
  int res2[MC_MAX][D_MAX];
  double rep, met;
  char str[255];
  int res3[D_MAX];
  int i, j, k;
  int rndseed = 42;
  double res4[D_MAX];
  int res5[D_MAX];
  FILE *outd;
  FILE *outs;
  FILE *outc;
  FILE *outh;
  FILE *outr;
  int hd;
  double repla;

  if(argc==5)
  {
    sscanf(argv[1], "%d", &D);
    sscanf(argv[2], "%d", &SD);
    sscanf(argv[3], "%lf", &mu);
    sscanf(argv[4], "%d", &rndseed);
  }else if(argc == 1)
  {
    D = 3;
    SD = 30;
    mu = 0.001;
    rndseed = 12249;
  }else
  {
    error(3, "ERROR -- 3: bad argument");
  }

  MC=SD-1;

  printf("seed=%d\n", rndseed);
  seed(rndseed);

  if(SD<=D)
    error(2, "ERROR -- 2: division densitiy too low");

  snprintf(str,sizeof(str),"s.%02d.%02d.%02d",D,SD,(int)(mu*1000));
  outs=fopen(str,"wt");
  str[0]='d';
  outd=fopen(str,"wt");
  str[0]='c';
  outc=fopen(str,"wt");
  str[0]='h';
  outh=fopen(str,"wt");
  str[0]='r';
  outr=fopen(str,"wt");
  if((outs == NULL) || (outd == NULL) || (outc == NULL) || (outh == NULL) || (outr == NULL))
    error(1, "output file open error");

  init();

  for(cou=0;cou<COU;cou++)
  {
    t=choose_according_fitnes();
    if(t==-1)
    {
      fclose(outd);
      fclose(outs);
      fclose(outc);
      str[0]='s';
      unlink(str);
      str[0]='d';
      unlink(str);
      str[0]='c';
      unlink(str);
      //error(1, "#kihalas");
      printf("died out...\n");
      return(1);
    }

    if(randd()<nu1)
      chromosomatization(t);
    if(randd()<nu2)
    {
      int r = chrbreak(t);
      if(r != 0)
        error(1, "# break error!");
    }

    while(replicate_random_template(t)!=0){;}
    if(comp[t].numtot>=SD)
        division(t);

    if((cou%1000)==0)
    {
      totnum = 0;
      for(t=0;t<COMP_NUM;t++)
        totnum += comp[t].numtot;

      for(t=0;t<MC;t++)
        num[t]=0.0;
      chr_num();
      fprintf(outd,"%d\t",cou/1000);
      for(t=0;t<MC;t++)
        fprintf(outd,"%lf\t",(double)num[t]*(t+1)/(double)totnum);
      fprintf(outd,"\n");
      fflush(outd);
      mutatedness(&rep, &met);
      fprintf(outs,"%d\t%le\t%lf\t%lf\t%lf\t%lf\t%lf\n",cou/1000,aver_fitn(),avg_chr_num(),chrom_div(),non_zero_fitness(), rep, met);
      fflush(outs);
      count_chr(res, res2);
      fprintf(outc,"%d\t",cou/1000);
      for(t = 0; t < MC; t++)
        fprintf(outc,"%lf ", res[t]);
      fprintf(outc,"\n");
      fflush(outc);

      for(t = 0; t < D; t++)
        res3[t] = 0;
      fprintf(outh,"%d\t",cou/1000);
      for(i = 0; i < COMP_NUM; i++)
        for(j = 0; j < comp[i].numtot; j++)
          if((comp[i].SEQ[j].head == 1) && (comp[i].SEQ[j].next != -1))
            res3[comp[i].SEQ[j].type]++;
      for(t = 0; t < D; t++)
        fprintf(outh, "%d\t", res3[t]);
      fprintf(outh, "\n");
      fflush(outh);

      for(t = 0; t < D; t++)
      {
        res4[t] = 0;
        res5[t] = 0;
      }
      fprintf(outr, "%d\t",cou/1000);
      for(i = 0; i < COMP_NUM; i++)
        for(j = 0; j < comp[i].numtot; j++)
          if((comp[i].SEQ[j].head == 1) && (comp[i].SEQ[j].next != -1))
          {
            hd=0;
            for(k=0;k<LR;k++)
              if(comp[i].SEQ[j].ER[k]!=REPLACT[k])
                hd++;
            repla = (1.0 - (pow((double)hd, alpha)) / (beta + pow((double)hd, alpha)));
            repla *= replaTbl[comp[i].SEQ[j].type];
            res4[comp[i].SEQ[j].type] += repla;
            res5[comp[i].SEQ[j].type]++;
          }
      for(t = 0; t < D; t++)
      {
        if(res5[t] > 0)
          fprintf(outr, "%lf\t", res4[t]/((double)res5[t]));
        else
          fprintf(outr, "0\t");
      }
      fprintf(outr, "\n");
      fflush(outr);
    }
  }

  if(cou==COU)
  {
    fprintf(outd,"#normal termination\n");
    fprintf(outs,"#normal termination\n");
    fprintf(outc,"#normal termination\n");
    fprintf(outh,"#normal termination\n");
    fprintf(outr,"#normal termination\n");
  }
  fclose(outd);
  fclose(outs);
  fclose(outc);
  fclose(outh);
  fclose(outr);

  return(0);
}
