/*
*
* File main.c
* Copyright 2019 J. Guo, P. Singh, K. E. Bassler
* This program implements the Reduced Network Extremal Ensemble Learning (RenEEL) scheme for maximizing modularity. 
* Instructions for usage are contained in the file readme.txt
*
* If you use it for your work, please cite the following paper:
* Reduced network extremal ensemble learning (RenEEL) scheme for community detection in complex networks,
* Jiahao Guo, Pramesh Singh & Kevin E. Bassler 
* Scientific Reports volume 9, Article number: 14234 (2019) 
* For comments/questions and reporting any bugs please contact Kevin E. Bassler (bassler@uh.edu)
*
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"rg.h"
#include"help.h"

#define Tsmall 1e-6


int *cols;


void inputGraph(struct graph *G);
int cleansort(struct part *ensem,int kmax,int N);
int getscore(int *score,int *edgecc, int N,int know,int *ref);
void renorm(struct graph G, int *ref);
int mergroup(struct graph *G, int x, int y);
void remolast(struct part *ensemble, int ens, int *edgecc, int N)
{
	long i;
	long j;
	for (i=0;i<N;i++)		
		cols[i]=0;
	
	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (ensemble[0].pa[i] == ensemble[0].pa[j])
				edgecc[i*N + j]--;
				if (edgecc[i*N+j]==ens-1)
				cols[j]=1;
			}
		
		}
	//delete in ensemble
	i = 0;
	free(ensemble[0].pa);
	for (i = 0;i<ens - 1;i++)
	{
		ensemble[i] = ensemble[i + 1];
	}
	ensemble[ens - 1].pa = NULL;
}

void replaceone(struct part* ensemble,int know, struct part* inensem, int max, int *edgecc, int N)
{
	long i,j;
	for (i=0;i<N;i++)	cols[i]=0;
	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (inensem[max].pa[i] == inensem[max].pa[j])
				edgecc[i*N + j]++;
				if (ensemble[0].pa[i]==ensemble[0].pa[j])
				edgecc[i*N+j]--;
				if (edgecc[i*N+j]==know)
				cols[j]=1;
			}
		}

	j = -1;
	for (i = 1;i<know;i++)
		if (ensemble[i].Q>inensem[max].Q)
		{
			j = i;    break;
		}

	if (j<0)
	{
		free(ensemble[0].pa);
		for (i=1;i<know;i++)
			ensemble[i-1]=ensemble[i];
		ensemble[know-1]=inensem[max];

	}
	else
	{
		free(ensemble[0].pa);
		for (i = 1;i < j;i++)
			ensemble[i - 1] = ensemble[i];
		ensemble[j-1]=inensem[max];


	}
}


void addone(struct part *ensemble, int know, struct part *inensem, int max, int *edgecc, int N)
{
	long i, j;
	for (i=0;i<N;i++)		cols[i]=0;

	for (i = 0;i<N-1;i++)
		if (!cols[i])
		{
			cols[i]=1;
			for (j = i+1;j<N;j++)
			{
				if (inensem[max].pa[i] == inensem[max].pa[j])
				edgecc[i*N + j]++;
				if (edgecc[i*N+j]==know+1)
				cols[j]=1;
			}
		}

	j = -1;
	for (i = 1;i<know;i++)
		if (ensemble[i].Q>inensem[max].Q)
		{
			j = i;    break;
		}

	if (j<0)
	{
		ensemble[know]=inensem[max];
	}
	else
	{
		for (i = know - 1;i >= j;i--)
		{
			ensemble[i + 1] = ensemble[i];
		}
		ensemble[j]=inensem[max];

	}

}
clock_t TIME;
void puttime(FILE *fo,char *mess)
{
	fprintf(fo,mess);
	fprintf(fo,": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
	printf(mess);
	printf(": %lf\n",(double)(clock()-TIME)/CLOCKS_PER_SEC);
	TIME=clock();
}



int main(int argc, char *argv[])
{
    clock_t T1;
	time_t sec;
	sec=time(NULL);
    double ti;
	TIME=clock();

    long i,j,k;
    int size;
	FILE *fout,*frc;
	fout=fopen("result.txt","w");
	frc=fopen("records.txt","a");

    size=1;
	fprintf(fout,"number of processors: %d\nclock per second: %d\n",size,CLOCKS_PER_SEC);
	printf("number of processors: %d\nclock per second: %d\n",size,CLOCKS_PER_SEC);

	unsigned int seed, *seedlist;
    int krg,know,copy1,copy2;
    int kmax,kp;

    krg=atoi(argv[1]);   
    copy1=atoi(argv[2]);
    kmax=copy1*size;
    copy2=atoi(argv[3]);
    kp=copy2*size;

	seed=(unsigned int)time(NULL);

	fprintf(fout,"seed=%ld\n",seed);
	fprintf(frc,"seed=%ld ",seed);
	fprintf(fout,"krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);
	fprintf(frc,"krg=%d kmax=%d kp=%d ",krg,kmax,kp);
	printf("seed=%ld\n",seed);
	printf("krg=%d\nkmax=%d\nkp=%d\n",krg,kmax,kp);

    seedlist=(unsigned int*)malloc(sizeof(unsigned int)*size);
    for (i=0;i<size;i++)
	{
		rand1(&seed);	
		seedlist[i]=seed+i;
	}
	struct graph G;
	inputGraph(&G);
	fprintf(fout,"links=%d nodes=%d\n",G.n,G.N);
	printf("links=%d nodes=%d\n",G.n,G.N);
	fprintf(frc,"links=%d nodes=%d ",G.n,G.N);
   	puttime(fout,"reading time");
 //input network

    struct part *ensemble,*iterensem;
    ensemble=(struct part*)malloc(sizeof(struct part)*kmax);
    iterensem=(struct part*)malloc(sizeof(struct part)*kp);

 printf("working on initialization\n");
 fprintf(fout,"working on initialization\n");
   
    int *edgecc, *score;
    score=(int*)malloc(sizeof(int)*G.N);
    edgecc=(int*)malloc(sizeof(int)*G.N*G.N);


	{
		

	    for (j=0;j<copy1;j++)
		{

	        	ensemble[j*size]=RG(G,krg,seedlist);
			
		//get set of partitions
		}
	}

	

	puttime(fout,"initialization time");
    know=cleansort(ensemble,kmax,G.N);
	puttime(fout,"sort time");

	cols=(int*)malloc(sizeof(int)*G.N);
	i=0;
	while (i<G.N) {cols[i]=0; i++;} 

    for (i=0;i<G.N-1;i++)
       if (!cols[i])
	{
		cols[i]=1;
            	for (j=i+1;j<G.N;j++)
		{
			edgecc[i*G.N+j]=0;
			for (k=0;k<know;k++)
                	if (ensemble[k].pa[i]==ensemble[k].pa[j])
                	edgecc[i*G.N+j]++;
    			if (edgecc[i*G.N+j]==know)
			cols[j]=1;
		}
	}	


    //initialization
    puttime(fout,"count edgecc time");

    int iteration=0;
    int *ref, bp;
    // ref: new group label for each old group (group=community)
    // score: current group label for each original node (not reduced node) 
    ref=(int*)malloc(sizeof(int)*G.N);
	for (i = 0;i < G.N;i++)
		score[i] = i + 1;

   
    T1=clock();
    while (know>1)
    {

        iteration++;
	
   	printf("%d:size=%d,ensem=%d,%lf,%lf\n",iteration,G.com,know,ensemble[0].Q/(2*G.n),ensemble[know-1].Q/(2*G.n));
   	fprintf(fout,"%d:size=%d,ensem=%d,%lf,%lf\n",iteration,G.com,know,ensemble[0].Q/(2*G.n),ensemble[know-1].Q/(2*G.n));

      int sizeofRN;
        sizeofRN=getscore(score,edgecc,G.N,know,ref);
        //update G, edgecc

	renorm(G,ref);
        G.com=sizeofRN;






	puttime(fout,"reducing time");

        //partitions of newG
        for (j=0;j<copy2;j++)
	{

		{
			
            		iterensem[j*size]=RG(G,krg,seedlist);

		}
       	}

	puttime(fout,"RG time");
		//pick the best
		bp = 0;	j = 1;
		for (i=1;i<kp;i++)
			if (abso(iterensem[i].Q - iterensem[bp].Q) < Tsmall)
			{
				j++;
				if (rand1(seedlist)*j < 1)
				{
					free(iterensem[bp].pa);	
					bp = i;
				}
				else
					free(iterensem[i].pa);
			}
			else if (iterensem[i].Q > iterensem[bp].Q)
			{
				j = 1;
				free(iterensem[bp].pa);
				bp=i;
			}
			else
			{
				free(iterensem[i].pa);
			}
    
        //update ensemble and edgecc
		i = -1;
		for (j = 0;j<know;j++)
			if (abso(iterensem[bp].Q -ensemble[j].Q )<Tsmall)
			{
				if (iterensem[bp].com == ensemble[j].com)
					if (comps(iterensem[bp].pa, ensemble[j].pa, G.N))
					{
						i = j;    break;
					}
			}

		if (i >= 0)
		{
	
			remolast(ensemble, know, edgecc, G.N);
			know--;
		}
		else
		{
			if (iterensem[bp].Q<ensemble[0].Q)
			{
	
				remolast(ensemble, know, edgecc, G.N); know--;
			}
			else
			{
				if (know<kmax)
				{

	
					addone(ensemble, know, iterensem, bp, edgecc, G.N);
					know++;
				}
				else
				{

					
					replaceone(ensemble, know, iterensem, bp, edgecc, G.N);
					
				}
			}
		}

	puttime(fout,"update ensemble time");



        //free iterensem
    }
 
	fprintf(fout,"rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
	fprintf(fout,"real walltime=%d\n",time(NULL)-sec);
	fprintf(frc,"real time=%d ",time(NULL)-sec);
	printf("rest time=%lf\n",(double)(clock()-T1)/CLOCKS_PER_SEC);
	printf("real walltime=%d\n",time(NULL)-sec);

	//get score of reduced network
	   int sizeofRN;
        sizeofRN=getscore(score,edgecc,G.N,know,ref);
        //update G, edgecc
        renorm(G,ref);
        G.com=sizeofRN;
 //   for (i=0;i<know;i++)
	outpart(ensemble[0],G.n, G.N);
	double Qfinal;
	Qfinal=compQ(G);
	printf("Qfinal=%lf\n",Qfinal/(2*G.n));
	fprintf(fout,"Qfinal=%lf\n",Qfinal/(2*G.n));
	fprintf(frc,"Qf=%lf\n",Qfinal/(2*G.n));

	fclose(frc);
	fclose(fout);
	return 0;
}


int mergroup(struct graph *G, int x, int y)
{
    int i;
    int *newlist,*newl2;
    
    
    newlist=(int*)malloc(sizeof(int)*(G[0].ml[x][0]+G[0].ml[y][0]+1));
  
    for (i=1;i<=G[0].ml[x][0];i++)	newlist[i]=G[0].ml[x][i];
    for (i=1;i<=G[0].ml[y][0];i++)	newlist[i+G[0].ml[x][0]]=G[0].ml[y][i];

    newlist[0]=G[0].ml[x][0]+G[0].ml[y][0];

    free(G[0].ml[x]);
    free(G[0].ml[y]);

    G->ml[x]=newlist;
    G->ml[y]=NULL;
    
    //printf("update1\n");
    newlist=(int*)malloc(sizeof(int)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    newl2  =(int*)malloc(sizeof(int)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    newl2[0]=G->wl[x][0]+G->wl[y][0];
    newlist[0]=0;
    for (i=1;i<=G[0].nl[x][0];i++)
    {
        if (G[0].nl[x][i]==y){	newl2[0]+=G[0].wl[x][i];    }
        else{
            newlist[0]++;
            newlist[newlist[0]]=G[0].nl[x][i];
            newl2[newlist[0]]=G[0].wl[x][i];
        }
    }

    for (i=1;i<=G[0].nl[y][0];i++)
        if (G[0].nl[y][i]!=x)
        {
            int sg,k1,k2;
            sg=G[0].nl[y][i];
            k1=findxiny(x,G[0].nl[sg]);
            k2=findxiny(y,G[0].nl[sg]);
            if (k1>0)
            {
                G[0].wl[sg][k1]+=G[0].wl[sg][k2];
                int test;
                test=findxiny(sg,newlist);
                if (test>0)
                    newl2[test]=G[0].wl[sg][k1];
                else
                {
                    printf("Error2!\n");
                    exit(0);
                }
                
                if (k2!=G[0].nl[sg][0])
                {
                    G[0].nl[sg][k2]=G[0].nl[sg][G[0].nl[sg][0]];
                    G[0].wl[sg][k2]=G[0].wl[sg][G[0].nl[sg][0]];
                }
                G[0].nl[sg][0]--;
                
            }
            else
            {
                if (k2>0)
                    G[0].nl[sg][k2]=x;
                else
                {
                    printf("Error3!\n");    exit(0);
                }
                newlist[0]++;
                newlist[newlist[0]]=sg;
                newl2[newlist[0]]=G[0].wl[sg][k2];
            }
        }
    free(G[0].nl[x]);	free(G[0].wl[x]);
    G->nl[x]=newlist;
    G->wl[x]=newl2;
    free(G[0].nl[y]);	free(G[0].wl[y]);
    G->nl[y]=NULL;
    G->wl[y]=NULL;

    
    G[0].dl[x]+=G[0].dl[y];

    return 0;
}

void movegroup(struct graph *G, int f, int t)
{
    int i,j;
    G->ml[t]=G->ml[f];
	G->ml[f] = NULL;
   
    for (i=1;i<=G->nl[f][0];i++)
    {
        j=findxiny(f,G->nl[G->nl[f][i]]);
        G->nl[G->nl[f][i]][j]=t;
    }
    
    G->nl[t]=G->nl[f];
    G->wl[t]=G->wl[f];
    G->dl[t]=G->dl[f];
	G->nl[f] = NULL;
	G->wl[f] = NULL;
}


void renorm(struct graph G,int *ref)
{
    int i,j;
    
    //renormlize G
    for (i=0;i<G.com;i++)
    {
        if (ref[i]!=i)
        {

			if (G.ml[ref[i]] == NULL)
				movegroup(&G,i,ref[i]);
			else
			mergroup(&G,ref[i],i);
        }
    }



    
}


int comps(int *p, int *q, int N)
{
    int i;
	for (i = 0;i < N;i++)
		if (p[i] != q[i])
			return 0;
    
    return 1;
}

int getscore(int *score,int *edgecc, int N,int know,int *ref)
{
    long i,j,k;
    int groupsnumber;
	int *newsc;
	newsc = (int*)malloc(sizeof(int)*N);
    for (i=0;i<N;i++)
        newsc[i]=0;
    groupsnumber=0;

    for (i=0;i<N;i++)
    {
        if (!newsc[i])
        {
            ref[score[i]-1]=groupsnumber;
            groupsnumber++;
            newsc[i]=groupsnumber;
        
//below may include in above section
        	for (j=i+1;j<N;j++)
        	{
 
           	 if (edgecc[i*N+j]==know)
            	{
			ref[score[j] - 1] = newsc[i]-1;
        		newsc[j]=newsc[i];
                
            	}
        	}
	}
    }
	for (i = 0;i < N;i++)
		score[i] = newsc[i];
	free(newsc);
    return (groupsnumber);
}

int cleansort(struct part *ensem, int kmax,int N)
{
    int i,j;
    struct part temp;
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (ensem[i].Q>ensem[j].Q)
            {
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
        }
    int *del,len;
    del=(int*)malloc(sizeof(int)*kmax);
    for (i=0;i<kmax;i++)
        del[i]=0;
    for (i=0;i<kmax-1;i++)
        if (del[i]==0)
        {
            j=i+1;
            while (j<kmax&&(ensem[j].Q-ensem[i].Q)<Tsmall)
            {
                if (ensem[j].com==ensem[i].com)
                    if (comps(ensem[i].pa,ensem[j].pa,N))
                        del[j]=1;
                j++;
            }
        }
    for (i=0;i<kmax-1;i++)
        for (j=i+1;j<kmax;j++)
        {
            if (del[i]>del[j])
            {
                del[i]=0;   del[j]=1;
                temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
            }
            else if (del[i]+del[j]==0)
            {
                if (ensem[i].Q>ensem[j].Q)
                {
                    temp=ensem[i];  ensem[i]=ensem[j]; ensem[j]=temp;
                }
            }
        }
    len=0;
    while (len<kmax&&del[len]==0) len++;
    i=len;
    while (i<kmax)
    {
        free(ensem[i].pa);
        i++;
    }
    free(del);
    
    return len;
}

void inputGraph(struct graph *G)
{

	
	int i,n1,n2;
	FILE *fi;
	fi=fopen("info.txt","r");
	fscanf(fi,"%d%d",&G->N,&G->n);
	fclose(fi);

	G->com=G->N;
	G->ml=(int**)malloc(sizeof(int*)*G->N);
	G->nl=(int**)malloc(sizeof(int*)*G->N);
	G->wl=(int**)malloc(sizeof(int*)*G->N);
	G->dl=(int*)malloc(sizeof(int)*G->N);

	fi=fopen("degree.txt","r");
	for (i=0;i<G->N;i++)
		fscanf(fi,"%d",G->dl+i);
	fclose(fi);

	fi=fopen("clean.txt","r");
	for (i=0;i<G->N;i++)
	{
		G->ml[i]=(int*)malloc(sizeof(int)*2);
		G->ml[i][0]=1;
		G->ml[i][1]=i+1;
		G->nl[i]=(int*)malloc(sizeof(int)*(G->dl[i]+1));
		G->nl[i][0]=0;
		G->wl[i]=(int*)malloc(sizeof(int)*(G->dl[i]+1));
		G->wl[i][0]=0;

	}


	for (i=0;i<G->n;i++)
	{
		fscanf(fi,"%d%d",&n1,&n2);
		n1--;	n2--;


		G->nl[n1][0]++;	G->nl[n2][0]++;
		G->nl[n1][G->nl[n1][0]]=n2;	
		G->wl[n1][G->nl[n1][0]]=1;
		G->nl[n2][G->nl[n2][0]]=n1;	
		G->wl[n2][G->nl[n2][0]]=1;
		
	}
	fclose(fi);

}


