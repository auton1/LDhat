#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#ifndef LENTYPE
#define LENTYPE (50) /* max number of different types min(Lm,(K+1)*(K+1)-1) */
#endif 


#ifndef RHO
#define RHO_Q (0) /*recombination rate for Qb matrices*/
#endif 


#ifndef EXP_TIMES
#define EXP_TIMES (25*(10)) /* terms in series approx of EXP */
#define QB_MAX (500) /*dim of Qb*/
#endif 


extern unsigned long seed[2];
extern void seed_set(FILE *a);
extern void seed_put(FILE *a);
extern double runif(void);
extern void priorpr_l_imp(double *logp, double *rem_p,int nt, int *n, int *rec,double *m,int n_pts,double *rho,double theta[2]);
extern void add_new(int *type,int *nd,int *new);
extern void remov(int *type, int *nd, int c);
extern double *dou_vec_init(const size_t num);
extern short *sh_vec_init(const size_t num);
extern void Qmat(double ****Q, double *****Qb, double theta[],int K,int L,int M, double *mu, double **P);

extern double **P;
extern long BS_LEN;



unsigned long seed[2];
void micro_satt(void);
void micro_satt_g(void);
void seed_set(FILE *a);
void seed_put(FILE *a);


double *dou_vec_init(const size_t num);
short *sh_vec_init(const size_t num);
double runif(void);
void priorpr_l_imp(double *logp, double *rem_p,int nt, int *n, int *rec,double *m,int n_pts,double *rho,double theta[2]);
void add_new(int *type,int *nd,int *new);
void remov(int *type, int *nd, int c);


/*S+ runif routine*/

double runif(void)
{
	unsigned long n,x, lambda=69069;
	*seed = *seed * lambda;
	*(seed+1) ^= *(seed+1) >> 15;
	*(seed+1) ^= *(seed+1) << 17;
	n = *seed ^ *(seed+1);
	x = ((n>>1) & 017777777777);
	return( (double)(x+0.5) / 2147483648. );
}


/***********************************
 *                                 *
 *  improved procal routine        *
 *  no logs prob= remp *log(logp)  *
 *                                 *
 ***********************************/

void priorpr_l_imp(double *logp, double *rem_p,int nt, int *n, int *rec,double *m,int n_pts,double *rho,double theta[2])
{ 
  int i,j,event;
  double norm,mut,reco,rate;
  double EXP25=72004899337.3859;
  mut=0;
  reco=0;
  if(*rec>0) reco=*rec;
  for(j=0;j<2;j++) mut += *(n+j)* *(theta+j);
    if(*(m+1)==0){
      event=0;
      rate=*m;
    }
    if(*(m+1)==1){
      event=1;
      rate=*m;
    }
    if(*(m+1)>=2){
      event=2;
      rate=*m* *(theta+(int)*(m+1)-2);
    }
   
  for(i=0;i<n_pts;i++){
    /* calculate denominator of prior prob = total rate out of state n(n-1) n.theta + rec.rho*/
    norm=nt*(nt-1)+*(rho+i)*reco+mut;
    
    /* rate =rate*v(event)/norm */
    if(event!=1)
      *(rem_p+i)*=(rate)/norm;
    else
      *(rem_p+i)*= rate* *(rho+i)/norm;
    if(*(rem_p+i)>0){
      while(*(rem_p+i)*EXP25<1){
	*(rem_p+i)*=EXP25;
	*(logp+i)-=25.0;
      }
    }
  }
}



/*********************************
 *                               *
 * Routine to read seed for file *
 *                               *
 *********************************/

void seed_set(FILE *file)
{
  long i,temp;
  unsigned long max=42949672;
  char ch;
 
  for(i=0;i<2;i++){
    /*ignore white space*/
    
    while(1){
      ch = fgetc(file);
      if(ch != '\n' && ch != '\n' && ch != '\t') break;
    }
    
    temp= ch-'0';
    
    while(1){
      ch = fgetc(file);
      if(ch == '\n' || ch == '\n' || ch == '\t' || ch ==EOF) break;
      temp= (10*temp+ch-'0') % max;
    }
    *(seed+i)=temp;
  }
}

/***********************
 *                     *
 * put seed in file    *
 * file_r and file_w   *
 * point to same file  *
 * but for read (r)    *
 * and write (w)       *
 *                     *
 ***********************/

void seed_put(FILE *file_w)
{
	(void)fprintf(file_w,"%lu\n%lu", *seed,  *(seed+1));
}   


/*add element new to sample*/
void add_new(int *type,int *nd,int *new)
{
	int i,j,c,flag;
	
	c=-1;
	for(i=0;i<*nd && c==-1;i++){
	  flag=1;
	  for(j=0;j<2 && flag==1;j++) flag= (*(type+(3)*i+j)==new[j]);
	  if(flag==1) c=i;
			}
	if(c>-1) *(type+(3)*c+2)+=1;
	else{
	  for(j=0;j<2;j++) *(type+(3)* *nd+j)=new[j];
	  *(type+(3)* *nd+2)=1;
	  *nd+=1;
	}
	return;
	
	
	
}

/*remove element *(type+3*i+(0,1,2) ) */
void remov(int *type, int *nd, int c)
{
	int i,j;
	
	if(*(type+(3)*c+2)>1){
	  *(type+(3)*c+2)-=1;
	  return;
	}
	else{
	  for(i=c;i<(*nd-1);i++){
	    for(j=0;j<3;j++) *(type+(3)*i+j)=*(type+(3)*(i+1)+j);
	  }
	  *nd-=1;
	  return;
	}
}

double *dou_vec_init(const size_t num)
{
  double *p;
  if((p = (double *) calloc(num, sizeof(double))) == NULL)
    perror("initialise_array():calloc failed!");



  return(p);
}




short *sh_vec_init(const size_t num)
{
  short *p;
  if((p = (short *) calloc(num, sizeof(short))) == NULL)
    perror("initialise_array():calloc failed!");



  return(p);
}

/**********************************************************

 * program to calcalate Q matrices for the K-allele case * 

 * L-loci, simple version                                *

 *********************************************************/

/*include NAG routines for inverting a matrix*/


#define ERROR_MIN (0.000000001)
#define TINY 1.0e-20;

void Qmat (double ****Q, double *****Qb, double theta[],int K,int L,int M, double *mu, double **P);
void matrix_inverse_new(double **a,long n);
void exp_matrix(double **out, double a,int K, double **P);


void Qmat (double ****Q, double *****Qb, double theta[],int K,int L,int M, double *mu, double **P)
{
  long i,j,l,k,ii,jj;
  double lambda,con, **Id,**temp;
  double s[]={0.3225476896192312,1.74576110115834658,4.53662029692112798,9.39507091230113313};
  
  
  Id = (double **)calloc(K,sizeof(double *));
  temp = (double **)calloc(K,sizeof(double *));
  for (i=0;i<K;i++){
    temp[i]=dou_vec_init(K);
    Id[i]= dou_vec_init(K);
  }		 

  
  /*generate Identity  matrix*/
  for(i=0;i<K;i++){
    for(j=0;j<K;j++){
      Id[i][j]= (i==j);
    }
  }
  
  /* do first entry- matrix of rows of mu */

  for(i=0;i<K;i++){
    for(j=0;j<K;j++){
      for(l=0;l<L;l++) Q[i][j][l][0]= mu[j];
    }
  }
 
  /* do other M entries */
  
  for(i=1;i<(M+1);i++){
    for(l=0;l<L;l++){
      lambda= *(theta+l)/((double)i+*(theta+l));
      
      /* calculate (I-lambdaP) */
      
      for(j=0;j<K;j++){
	for(k=0;k<K;k++){
	  temp[k][j]=Id[k][j]-lambda*P[k][j];
	}
      }
      
      matrix_inverse_new(temp,K); /* invert - to temp */    
      
       /* update Q - (1-lambda)*temp */
      for(j=0;j<K;j++){
	for(k=0;k<K;k++){
	  Q[k][j][l][i]=(1.0-lambda)*temp[k][j];
	}
      }
    }
  }
   
  /*calculate Qb matrix*/
  for(i=0;i<(QB_MAX);i++){
    for(j=0;j<4;j++){
      for(l=0;l<L;l++){
	lambda = (i+1); /*should include a rho-which one (all different
			  possibiliities?!)*/
	/* calculate exp_matrix( - theta s / lambda *(I- P)) */
	exp_matrix(temp,*(theta+l)*s[j]/lambda,K,P);
	con=exp(-1* (*(theta+l)*s[j]/lambda));
	for(ii=0;ii<K;ii++){
	  for(jj=0;jj<K;jj++){
	    Qb[ii][jj][l][j][i]= temp[ii][jj]*con;
	  }
	}
	
      }
    }
    /*(void)fprintf(stderr,"Calc matrix: %d",i+1);*/
  }
  for(i=0;i<K;i++){
    free(Id[i]);
    free(temp[i]);
  }
  free(Id);
  free(temp);

} 


		       
	
void exp_matrix(double **out, double a,int K, double **P)
{
  double **temp,**mult,**temp2, **out_temp,err,sum;
  long i,j,l,ii;

  temp =  (double **)calloc(K,sizeof(double *));
  for (i=0;i<K;i++) temp[i]= dou_vec_init(K);

  temp2 =  (double **)calloc(K,sizeof(double *));
  for (i=0;i<K;i++) temp2[i]= dou_vec_init(K);

  mult =  (double **)calloc(K,sizeof(double *));
  for (i=0;i<K;i++) mult[i]= dou_vec_init(K);

  out_temp  = (double **)calloc(K,sizeof(double *));
  for (i=0;i<K;i++) out_temp[i]= dou_vec_init(K);
 
  /* initial temp is the identity */
  for (i=0;i<K;i++){
    for (j=0;j<K;j++){
      temp[i][j]= (i==j);
      out_temp[i][j]= temp[i][j];
    }
  }

  /* mult is the matrix in the exponent */

  for (i=0;i<K;i++){
    for (j=0;j<K;j++){
      mult[i][j]=a* P[i][j];
      /*mult[i][j]*= -a;*/
    }
  }

  /*approx exponent*/

  err=1;
  for(l=1;l<EXP_TIMES && err>(ERROR_MIN);l++){
    /* calculate mult*temp/l =mult^l/l!*/
    for (i=0;i<K;i++){
      for (j=0;j<K;j++){
      temp2[i][j]=0;
      for( ii=0;ii<K;ii++){
	temp2[i][j]+=mult[i][ii]*temp[ii][j]/((double)l);
      }
      }
    }
    err*= a/((double)l);
    /* set temp=temp2 out +=temp2/(l!) */
     for (i=0;i<K;i++){
      for (j=0;j<K;j++){
	temp[i][j]=temp2[i][j];
	out_temp[i][j] +=temp2[i][j];
	}
     } 
  }

  
  for (i=0;i<K;i++){
    for (j=0;j<K;j++){
      out[i][j] =(double)out_temp[i][j];
    }
    
  }
  
  
  
  for(i=0;i<K;i++){
    free(temp[i]);
    free(temp2[i]);
    free(out_temp[i]);
    free(mult[i]);
  }
  free(temp);
  free(temp2);
  free(out_temp);
  free(mult);
}

/*inefficient Gauss-Jordan elimination*/
void matrix_inverse_new(double **a,long n)
{
  long i,j,k,imax;
  double big,temp;
  double **I;

  /*set-up I*/
  I = (double **)calloc(n,sizeof(double *));
  for(i=0;i<n;i++) I[i]=dou_vec_init(n);

  for(i=0;i<n;i++) for(j=0;j<n;j++) I[i][j]=0.0;
  for(i=0;i<n;i++) I[i][i]=1.0;

  /*loop over columns*/

  for(j=0;j<n;j++){
    /*firstly find largest a[i][j], i>=j*/
    big=0.0;
    for(i=j;i<n;i++){
      if(fabs(a[i][j])>big){
	big=fabs(a[i][j]);
	imax=i;
      }
    }
    if(big<=0.0){
      (void)fprintf(stderr,"Error in matrix_inverse routine: singular matrix? \n");
      abort();
    }

    /*swap rows*/
    for(k=0;k<n;k++){
      temp=a[j][k];
      a[j][k]=a[imax][k];
      a[imax][k]=temp;
      temp=I[j][k];
      I[j][k]=I[imax][k];
      I[imax][k]=temp;
    }

    /*divide row by a[j][j]*/
    temp=a[j][j];
     for(k=0;k<n;k++){
       a[j][k]/=temp;
       I[j][k]/=temp;
     }

     /*calculate jth column*/
     for(i=0;i<n;i++){
       if(j!=i){
	 temp=a[i][j];
	 for(k=0;k<n;k++){
	   a[i][k]-=temp*a[j][k];
	   I[i][k]-=temp*I[j][k];
	 }
       }
     }
  }

  /*replace A by I*/
  for(i=0;i<n;i++) for(j=0;j<n;j++) a[i][j]=I[i][j];

  /*free I*/
  for(i=0;i<n;i++) free(I[i]);
  free(I);
}


double ***pi; /* prob -factor into pis */
int kall_l(int *type, int *nd,double *rho,double theta[2], double *logq, double *logp,double *rem_p,int k, int n_pts,int K,int L,int M,double *mu,double **P, double ****Q, double *****Qb);
void update(int *type, double rho, double theta[2],int k, int *nd, int *nt, int *n,double *logq,double *m,int *rec,int K,int L, double **P,double ****Q, double *****Qb);
void procal_l(int *type, int nd, int k, int *n, int K, int L, double *mu,double ****Q, double *****Qb);
double pical(int *type,int k, int *new, int nt, int nd, double rho,int L,double ****Q, double *****Qb);


int kall_l(int *type, int *nd,double *rho,double theta[2], double *logq, double *logp,double *rem_p,int kk, int n_pts,int K,int L,int M,double *mu,double **P, double ****Q, double *****Qb)
{
  /* Remember P and Q matrices global. K, L and M are defined */
  
  short int flag=0; /* continue loop? */
  int min,max;
  int *n,*rec,i,j,nt=M; /*number of individuals ancestral at loci 1,2,..*/
  int *n_st,*rec_st,nt_st;
  double u,tot,m[2],fac; /* runif, sum of `probabilities', m[0] is prior rate, m[1] is type of transition */
  double pro[LENTYPE]; /*probalities -event happening to ind, and approx cond post (loci (L) changed to allele (K)) */
  unsigned long seed_store[2]; /* store initial seed- print out if error */
  double time;/*,length,time_sq,length_sq*/
  
  n=(int *)calloc(L,sizeof(int));
  rec=(int *)calloc(L-1,sizeof(int));
  n_st=(int *)calloc(L,sizeof(int));
  rec_st=(int *)calloc(L-1,sizeof(int));
  pi = (double ***) calloc(K,sizeof(double **));
  for(i=0;i<K;i++){
    pi[i]= (double **) calloc(L,sizeof(double *));
    for(j=0;j<L;j++) pi[i][j]=(double *)calloc(4,sizeof(double));
  }


  
  /* set logq and logp to 0 rem_p to 1*/
  *logq=0;
  for(i=0;i<n_pts;i++){
    *(logp+i)=0;
    *(rem_p+i)=1;
  }

  *seed_store = *seed;
  *(seed_store+1) = *(seed+1);

  /* initiate n and rec */
  for(i=0;i<L;i++) *(n+i)=M;
    
  for(i=0;i<L-1;i++) *(rec+i)=M;
    
#ifdef TEST
  fprintf(stderr,"\n seed %u \t %u",*(seed),*(seed+1));
#endif /*TEST*/
  
  while(!flag){
     
    /*store nt,n,rec for use in priorpr_l*/
    for(i=0;i<L;i++) *(n_st+i)=*(n+i);
    for(i=0;i<L-1;i++) *(rec_st+i)=*(rec+i);
    nt_st=nt;
    /*calculate probability of event happening to 1,...*ndth type*/
    tot=0;
    for(i=0;i<*nd;i++){
      *(pro+i)= *(type+(L+1)*i+L); /*number*/
      fac= (nt-1); /*coalescence*/
      min=-1;
      max=-1;
      for(j=0;j<L;j++){
	if( *(type+(L+1)*i+j)>0){
	  fac+= *(theta+j);
	}
      }
      if( *(type+(L+1)*i)>0 && *(type+(L+1)*i+1)>0)
	fac+= rho[kk];
      *(pro+i) *= fac;
      tot+= *(pro+i);
    }
   /* error */
    if(tot<=0){
      (void)fprintf(stdout,"Error log of non-positive no (tot) . \n");
      (void)fprintf(stdout,"initial seed %lu \t %lu \n\n",*(seed_store), *(seed_store+1));
      return(1);
    }
    
    /* choose individual. (type+3*i+(0,1,2,..L)) is individual chosen */ 
    u=tot*runif();
    i=0;
    while(u>0){
      u-=*(pro+i);
      i++;
    }
    i--;

    /* update logq */
    if(*(pro+i)>0){
      *logq+=log(*(pro+i)/(tot));
    }else{
      (void)fprintf(stdout,"Error log of non-positive no. (pro) \n");
      (void)fprintf(stdout,"initial seed %lu \t %lu \n\n",*(seed_store), *(seed_store+1));
      return(1);
    }
    
    /* calculate approximation of conditional posteriors-single loci */
    procal_l(type,*nd,i,n,K,L,mu,Q,Qb);
    /*choose event and update type nd m  and logq */
    update(type,rho[kk],theta,i,nd,&nt,n,logq,m,rec,K,L,P,Q,Qb);
    
    /*error*/
    if(*logq==2 || *logq==1){
      (void)fprintf(stdout,"Error in update logq = %f \n seed %lu %lu \n n %d %d \n",*logq,*seed_store,*(seed_store+1),*n,*(n+1));
      return(1);
      
    }
    
    /*update logp*/
     priorpr_l_imp(logp,rem_p,nt_st,n_st,rec_st,m,n_pts,rho,theta);
    
#ifdef TEST
    fprintf(stderr,"\n1 %d \t2 %d \tancestral %d \tdouble %d ",*n,*(n+1),nt,*rec);
      fprintf(stderr,"logq %f logp %f",*logq,*(logp)+log(*(rem_p)));
    
#endif /*TEST*/

    if(*nd>LENTYPE) return(1);
    /* update flag */
    flag=1;
    for(j=0;j<L;j++){
      flag *= (*(n+j)==1);
    }
  }
  /* update logq to take account of mutation stationary distance*/
  for(i=0;i<*nd;i++){
    for(j=0;j<L;j++){
      if(*(type+(L+1)*i+j)>0) *logq  -= log(mu[ *(type+(L+1)*i+j)-1 ]);
    }
  }
  /*TMRCA=time;*/
  free(n);
  free(rec);
  free(n_st);
  free(rec_st);
  for(i=0;i<K;i++){
    for(j=0;j<L;j++) free(pi[i][j]);
    free(pi[i]);
  }
  free(pi);
  return(0);
}

/* calculate pi[K][L]- independent loci. input types, number distinct,
   */
void procal_l(int *type, int nd, int k, int *n, int K, int L, double *mu,double ****Q, double *****Qb)
{
  int i,j,m,l, *n_copy;

  n_copy=(int *)calloc(L,sizeof(int));

  /*copy n and omit individual k from the sample */
  for(i=0;i<L;i++){
    n_copy[i]=n[i];
    if (*(type+(L+1)*k+i)>0) n_copy[i]-=1;
  }

  *(type+(L+1)*k+L)-=1;

  for(i=0;i<L;i++){
    /*init pi*/
    for(j=0;j<K;j++) for(l=0;l<4;l++) pi[j][i][l]=0;

    if(n_copy[i]>=QB_MAX) n_copy[i]=QB_MAX-1;
    
    /* providing there are some ancestral types at locus i- calc pi */
    if(n_copy[i]>0){
      for(m=0;m<nd;m++){
	if(*(type+(L+1)*m+i)>0){
	  for(j=0;j<K;j++){
	    for(l=0;l<4;l++){
	      pi[j][i][l]+= *(type+(L+1)*m+L) * Qb[*(type+(L+1)*m+i)-1][j][i][l][n_copy[i]];
	    }
	  }
	}
      }
      for(j=0;j<K;j++) for(l=0;l<4;l++) pi[j][i][l]/= (double)(n_copy[i]);
    }
    else{ /*else use stationary probs */
      for(j=0;j<K;j++) for(l=0;l<4;l++) pi[j][i][l]= mu[j];
    }
  }

  /* add individual k back to sample */
  *(type+(L+1)*k+L)+=1;

  free(n_copy);
}

void update(int *type, double rho, double theta[2],int k, int *nd, int *nt, int *n,double *logq,double *m,int *rec,int K,int L, double **P,double ****Q, double *****Qb)
{
  double *pr,u,tot,pitemp,pitemp2; /*event prob*/
  int i,j,a,e,*new,min,max,omin,omax,*rec1,*rec2,*old;
  
  rec2=(int *)calloc(L,sizeof(int));
  rec1=(int *)calloc(L,sizeof(int));
  new=(int *)calloc(L,sizeof(int));
  old=(int *)calloc(L,sizeof(int));
  

  pr=dou_vec_init(K*L+L-1+LENTYPE);
  
  for(i=0;i<L;i++) {
    new[i]=*(type+(L+1)*k+i);
    rec1[i]=0;
    rec2[i]=new[i];
  }

  /*calc prob of each event prob*pi(type) */

  /*mutations*/
  for(i=0;i<L;i++){ /*Locus*/
    if(*(type+(L+1)*k+i)>0){
      a=*(type+(L+1)*k+i)-1;
      for(j=0;j<K;j++){
	if(P[j][a]>0){
	  new[i]=j+1;
	  pitemp=pical(type,k,new,*nt,*nd,rho,L,Q,Qb);
	  new[i]=a+1;
	  pr[K*i+j]= *(theta+i) * pitemp *P[j][a];
	  }
	else
	  pr[K*i+j]=0;
	
      }
    }
    else
      for(j=0;j<K;j++) pr[K*i+j]=0;
  }

  /*recombination*/
  min=-1;
  max=-1;
  for(i=0;i<L;i++){
    pr[K*L+i]=0;
    if(*(type+(L+1)*k+i)>0){
      max=i;
      if(min<0) min=i;
    }
  }
  
    
 
    for(i=min;i<max;i++){
    rec1[i]=rec2[i];
    rec2[i]=0;
    pitemp=pical(type,k,rec1,*nt,*nd,rho,L,Q,Qb);
    pitemp2=pical(type,k,rec2,*nt,*nd,rho,L,Q,Qb);
    pr[K*L+i]=rho*pitemp*pitemp2;
    }

    /*coalescence*/
  
  for(i=0;i<*nd;i++){
    if(i!=k){
      pr[K*L+L-1+i]=0;
      for(j=0;j<L;j++){
	if(*(type+(L+1)*i+j)==*(type+(L+1)*k+j))
	  new[j]=*(type+(L+1)*i+j);
	else{
	  if(*(type+(L+1)*i+j)>0 && *(type+(L+1)*k+j)>0) j=L+1;
	  else{
	    if(*(type+(L+1)*i+j)>0) new[j]=*(type+(L+1)*i+j);
	    else 
	      new[j]=*(type+(L+1)*k+j);
	  }
	 }
	if(j<L) old[j]=*(type+(L+1)*i+j);
      }
      if(j==L){
	pitemp=pical(type,k,new,*nt,*nd,rho,L,Q,Qb);
	pitemp2=pical(type,k,old,*nt,*nd,rho,L,Q,Qb);
	pr[K*L+L-1+i]= *(type+(L+1)*i+L) * pitemp / pitemp2;
      }
    }
    else pr[K*L+L-1+i]= *(type+(L+1)*i+L)-1;
  }

  
  /*update sample and logq*/
   tot=0;
  for(e=0;e<(L*K+L+*nd-1);e++){
    tot+=*(pr+e);
  }
  
  if (tot<=0){
    *logq=1;
    return;}
  
  u=runif() * tot;
  e=0;
  while(u >= 0){
    u-= *(pr+e);
    e++;
  }
  
  if (*(pr+e-1)>0) {*logq+= log( *(pr+e-1)/(tot));
  }else{
    (void)fprintf(stdout,"\n e %d prob %f tot %f\n ",e,*(pr+e-1),tot);
    (void)fprintf(stdout,"\n types \n \n");
    for(i=0;i<*nd;i++){
      for(j=0;j<(L+1);j++) (void)fprintf(stdout,"%d ",*(type+(L+1)*i+j));
      (void)fprintf(stdout,"\n");
    }
    
    *logq= 2;
    return;
  }


  /*updates for each event, update type, nd, n, nt*/

  if(e<= K*L){ /*mutation*/
    a=(e-1)%K; /*new allele type-1*/
    i=(e-a-1)/K; /*locus*/
    
    /*update m*/
    m[0]=*(type+(L+1)*k+L)* P[a][*(type+(L+1)*k+i)-1];
    m[1]=L+i; 
    
    /*calculate new type*/
    for(j=0;j<L;j++)  new[j]=*(type+(L+1)*k+j);
    new[i]=a+1;

    add_new(type,nd,new);
    
    

    /*remove old type-update nd*/
    remov(type,nd,k);
    
    
  }

  if(e>K*L && e <=(K*L+L-1)){ /*recombination*/
    i=e-K*L; /*locus*/
   /*update m*/
    m[0]= *(type+(L+1)*k+L);
    m[1]=i;
    /*calculate new type*/

    max=-1;
    for(j=0;j<i;j++){ 
      new[j]=*(type+(L+1)*k+j);
      if(*(type+(L+1)*k+j)>0) min=j;}
    for(j=(i);j<L;j++){
      new[j]=0;
      if(max==-1 && *(type+(L+1)*k+j)>0 ) max=j;
    }
    /*update rec*/
    for(j=min;j<max;j++) *(rec+j)-=1;
    /*add new type*/
    add_new(type,nd,new);

    for(j=0;j<i;j++) new[j]=0;
    for(j=(i);j<L;j++) new[j]=*(type+(L+1)*k+j);

    /*add new type*/
    add_new(type,nd,new);

    /*remove old  type*/
    remov(type,nd,k);

    *nt+=1;
   
    
  }
  
  if(e>(K*L+L-1)){ /*coalesence*/
    i=e-(K*L+L); /*individual*/

    /*update rec*/
    min=-1;
    omin=-1;
    for(j=0;j<L;j++){
      if(*(type+(L+1)*i+j)>0){
	omax=j;
	if(omin==-1) omin=j;
      }
      if(*(type+(L+1)*k+j)>0){
	max=j;
	if(min==-1) min=j;
      }
    }
    min= min+(omin-min)*(omin>min);
    max= max+(omax-max)*(omax<max);
    if(min<=max) for(j=min;j<max;j++) *(rec+j)-=1;
    else for(j=max;j<min;j++) *(rec+j)+=1;

    
    for(j=0;j<L;j++){
      if(*(type+(L+1)*k+j)==0) new[j]= *(type+(L+1)*i+j);
      else{
	new[j]= *(type+(L+1)*k+j);
	if(  *(type+(L+1)*i+j)>0) n[j]-=1;
      }
    }
    /*calc m*/
    m[0]= *(type+(L+1)*k+L) *( *(type+(L+1)*i+L) - (i==k));
    m[1]=0.0;
    /*add new type*/
    add_new(type,nd,new);
    /*remove two types*/
    if ( *(type+(L+1)*k+L)==1 && k<i) i--;
      remov(type,nd,k);
    remov(type,nd,i);
    
    *nt-=1;
    
  }
  free(new);
  free(rec2);
  free(pr);
  free(rec1);
  free(old);
}
  
  
/*prob of nth member new given n-1 members (type - individual k).
  nt=n, nd = number of distinct individuals-NEW METHOD*/
double pical(int *type,int k, int *new, int nt, int nd,double rho,int L,double ****Q, double *****Qb)
{
  int i,j,l,*index,index_len=0,a,e,d;
  double temp=0.0;
  double int_temp;
  double w[]={0.6031541043416, 0.3574186924378, 0.0388879085150, 0.0005392947056}; /*weights for num int*/
  double rec;
  double p[4][LENTYPE]; 


  index=(int *)calloc(L,sizeof(int));

  /*remove type k*/
  nt--;
  *(type+k*(L+1)+L)-=1;


  /*set up index-loci ancetral in new*/
  for(i=0;i<L;i++){
    if(*(new+i)>0){ *(index+index_len)=i;
    index_len++;}
  }

  if(index_len==0){
    /*add type k back*/
  *(type+k*(L+1)+L)+=1;
  free(index);
  return(1.0);
  }
  if(index_len==1){
    /*add type k back*/
  *(type+k*(L+1)+L)+=1;
  temp=0;
  for(l=0;l<4;l++) temp+=w[l]*pi[*(new+*index)-1][*index][l];
  free(index);  
  return(temp);
  }
  
  for(i=0;i<4;i++) for(j=0;j<nd;j++) p[i][j]=1.0;
  temp=1.0;

  for(i=0;i<index_len;i++){
    /*calc recombination rate*/
    if(i>0){
      rec=0.0;  
      for(j=index[i-1];j<*(index+i);j++) rec+=rho;
      rec=rec/(rec+nt);
    }
    else rec=0.0;

    if((nt+(int)(rec+.5))>QB_MAX) a=QB_MAX-1;
    else
    a=nt+(int)(rec+.5)-1;
    d=new[index[i]];
    for(j=0;j<nd;j++){
      if(*(type+j*(L+1)+L)>0){
	e=*(type+j*(L+1)+index[i]);

	if(e>0){
	  for(l=0;l<4;l++) p[l][j]=((1-rec)*p[l][j]+rec*temp)*Qb[e-1][d-1][index[i]][l][a];
	}
	else{
	  for(l=0;l<4;l++) p[l][j]=((1-rec)*p[l][j]+rec*temp)*pi[d-1][index[i]][l];
	}
      }
    }
    temp=0;
    for(j=0;j<nd;j++) for(l=0;l<4;l++) temp+=w[l]* *(type+(L+1)*j+L)*p[l][j];
    temp/=nt;
  }

  /*add type k back*/
  *(type+k*(L+1)+L)+=1;
  free(index);

  return(temp);

  
}


void pairs(int nd, int *data, int K, double **P, double *mu, double theta[2], double rho_max, int n_pts, int NRUN, double *log_lik);

extern int kall_l(int *type, int *nd,double *rho,double theta[2], double *logq, double *logp,double *rem_p,int k, int n_pts,int K,int L,int M,double *mu,double **P, double ****Q, double *****Qb);
/***********************************************************
 *
 * Routine for calculating pairwise likelihood curve.
 * Input:
 *
 *  nd_store = number of distinct haplotypes
 *  data = vector of length 3*nd; type at each locus
 *         and multiplicity (type = {1,..,K})
 *  K = number of alleles
 *  P = mutation matrix; K by K
 * mu = mutation stationary distribution 
 * theta = mutation rate at each locus
 * rho_max = maximum value of rho
 * n_pts = number of points on likelihood surface;
 *         even,y spaced between 0.00001 and rho_max;
 *         each point is used to generate trees (estimate
 *         calculated using bridge sampling)
 * NRUN = number of trees to be generated
 * log_lik = output; NRUN length vector of log-likelihoods
 *
 *********************************************************/

#define ESS_MIN 1000
#define CF  1


void pairs(int nd_store, int *data, int K, double **P, double *mu, double theta[2], double rho_max, int n_pts, int NRUN, double *log_lik)
{
  FILE *seed_file,*seed_file_w;
  double *logp,*rem_p,logq;
  double *rho;
  double ****Q;
  double *****Qb;
  double norm;

  double *res_sq;
  double ess_min=0.0;

  int FLAG;
  int M; /*number of haplotypes*/
  int L=2; /*number of loci =2 for pairwise*/
  int nd;
  
  int i,j,k,l;
  int *type; /*copy of data*/

  /*calculate M*/
  M=0;
  for(i=0;i<nd_store;i++)
    M+=data[i*(L+1)+L];
  /*storage allocation*/

  logp = (double *)calloc(n_pts,sizeof(double));
  rem_p = (double *)calloc(n_pts,sizeof(double));
  rho = (double *)calloc(n_pts,sizeof(double));
  res_sq = (double *)calloc(n_pts,sizeof(double));

  Q = (double ****) calloc(K,sizeof(double ***));
  for(i=0;i<K;i++){
    Q[i] = (double ***) calloc(K,sizeof(double **));
    for(j=0;j<K;j++){
      Q[i][j] = (double **) calloc(L,sizeof(double *));
      for(l=0;l<L;l++)
        Q[i][j][l] = (double *)calloc(M+1,sizeof(double));
    }
  }
  
  Qb = (double *****) calloc(K,sizeof(double ****));
  for (i=0;i<K;i++){
    Qb[i] = (double ****) calloc(K,sizeof(double ***));
    for(j=0;j<K;j++){
      Qb[i][j] = (double ***) calloc(L,sizeof(double **));
      for(l=0;l<L;l++){
        Qb[i][j][l] = (double **) calloc(4,sizeof(double *));
        for(k=0;k<4;k++) Qb[i][j][l][k] = (double *)calloc(QB_MAX,sizeof(double));
      }
    }
  }
     
  type=(int *)calloc(3*LENTYPE,sizeof(int));
  /*initiate vectors*/

  for(i=0;i<n_pts;i++) {
    rho[i]= ((double)i*rho_max)/(double)(n_pts-1.0);
    log_lik[i]=0.0; /*initially store as likelihood; then normalise*/
    res_sq[i]=0.0;
  }

  
  /* initialise seed */
  seed_file = fopen(".seed","r");
  if (seed_file==NULL){
    (void) fprintf(stderr,"Can not open .seed\n Using default values.");
    *seed=1984795325;
    *(seed+1)=298376915;
  }else  {
    seed_set(seed_file);
    (void)fclose(seed_file);
  }
  
  
  /*calculate the Q matrices*/
  Qmat(Q,Qb,theta,K,L,M,mu,P);

  /*loop*/
  k=0; /*value to use for generating trees*/
  for(l=0;(l<NRUN) && (exp(ess_min)<ESS_MIN);l++){
/*    printf("\nrun %i: ess_min = ", l);*/
/*  k++*/
    if(k==n_pts) k=0;
    
    nd=nd_store;
    for(i=0;i<3*nd;i++)
      type[i]=data[i];

    /*generate tree*/
    FLAG=kall_l(type,&nd,rho,theta,&logq,logp,rem_p,k,n_pts,K,L,M,mu,P,Q,Qb);

    if(FLAG==0){
      /*store results*/
      norm=0.0;
      for(i=0;i<n_pts;i++) if (rem_p[i] > 0.0)
	norm += exp(logp[i]-logp[k])*rem_p[i];
      norm *= 1.0/((double)n_pts*rem_p[k] );
      for(i=0;i<n_pts;i++) {
/*		printf("\ni=%i k=%i rem=%e logp=%e logq=%e norm=%e", i, k, rem_p[i], logp[i], logq, norm);*/
		if (rem_p[i] > 0.00) {
			log_lik[i]+= rem_p[i]*exp(logp[i]-logq)/norm;
			res_sq[i] += CF*rem_p[i]*rem_p[i]*exp(2.0*logp[i]-2.0*logq)/(norm*norm);
/*			printf(" newlk=%e rsq=%e",log_lik[i],res_sq[i]);*/
		}
      }
      
      ess_min=-100;
      for(i=0;i<n_pts;i++)
	if ((log_lik[i] > 0.0) && (res_sq[i]> 0.0))
		if(CF*log_lik[i]*log_lik[i]/res_sq[i] > ess_min)
	  		ess_min=log(CF)+2*log(log_lik[i])-log(res_sq[i]);
/*      printf(" %f",exp(ess_min));*/
      if((l+1)%1000==0){
	(void)fprintf(stderr,"End of Run %i\n",l+1);
	(void)fprintf(stderr,"min ESS = %f \n",exp(ess_min));
      }
    }
    else{
      (void)fprintf(stderr,"Number of distinct individuals too large, ignoring tree \n");
      l--;
    }
    k++;
  }
  for(i=0;i<n_pts;i++)
    log_lik[i]=log(log_lik[i]/l); 
  
  free(logp);
  free(rem_p);
  free(rho);
  free(type);

  for(i=0;i<K;i++){
    for(j=0;j<K;j++){
      for(k=0;k<L;k++){
        free(Q[i][j][k]);
        for(l=0;l<4;l++){
          free(Qb[i][j][k][l]);
        }
        free(Qb[i][j][k]);
      }
      free(Q[i][j]);
      free(Qb[i][j]);
    }
    free(Q[i]);
    free(Qb[i]);
  }
  free(Q);
  free(Qb);
  free(res_sq);
  /*change .seed file*/
  seed_file_w=fopen(".seed","w");
  seed_put(seed_file_w);
  (void)fclose(seed_file_w);

}






