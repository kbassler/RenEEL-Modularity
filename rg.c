#include<stdio.h>
#include<stdlib.h>
#include"body.h"
#include"help.h"

int findxiny(int x, int *y)
{
    int j;
    for (j=1;j<=y[0];j++)
    
    if (y[j]==x) return (j);
    return (-1);
}


int outpart(struct part ans,int n, int N)
{
    int i;
    printf("%d\t%lf\n",ans.com,ans.Q/(2*n));
	FILE *fo;
	fo=fopen("partition.txt","w");
    	for (i=0;i<N;i++)
    		fprintf(fo,"%d\n",ans.pa[i]);
	fclose(fo);
    return 0;
}

double wex(struct graph G, int x, int y)
{
    int i,r;

    if (G.nl[x][0]>G.nl[y][0]) r=y;	else r=x;
    for (i=1;i<=G.nl[r][0];i++)
    if (G.nl[r][i]==x+y-r)
    return (G.wl[r][i]);


	return 0;

    
    
}


int updateG(struct graph *G, int x, int y)
{
    int i,k;
    int *newlist,*newl2;

    newlist=(int*)malloc(sizeof(int)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    newl2  =(int*)malloc(sizeof(int)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    k=0;
    newl2[0]=G->wl[x][0]+G->wl[y][0];
    for (i=1;i<=G[0].nl[x][0];i++)
    {
        if (G[0].nl[x][i]==y){	k=-1;	newl2[0]+=G[0].wl[x][i];}
        else{
            newlist[i+k]=G[0].nl[x][i];
            newl2[i+k]=G[0].wl[x][i];
        }
    }
    newlist[0]=G[0].nl[x][0]-1;
 
    //copy g1, g1-g2 becomes self-loop
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

                newl2[test]=G[0].wl[sg][k1];
   
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
    G->nl[y]=NULL;	G->wl[y]=NULL;
 
    if (y!=G[0].com-1)
    {
        for (i=1;i<=G[0].nl[G[0].com-1][0];i++)
        {
            int k1,k2;
            k1=G[0].nl[G[0].com-1][i];
            k2=findxiny(G[0].com-1,G[0].nl[k1]);

            G[0].nl[k1][k2]=y;
        }
        G->nl[y]=G->nl[G->com-1];
        G->wl[y]=G->wl[G->com-1];
        G->nl[G->com-1]=NULL;
        G->wl[G->com-1]=NULL;

    }
    
    G[0].dl[x]+=G[0].dl[y];
    if (y!=G[0].com-1)
    {
        G[0].dl[y]=G[0].dl[G[0].com-1];
    }
 
    return 0;
}

void extrafG(struct part ans,struct graph G,int *s)
{
    int i,j;

    for (i=0;i<G.com;i++)
    {
	if (G.ml[i]!=NULL)
		{
        	for (j=1;j<=G.ml[i][0];j++)
		if (G.ml[i][j]<=G.N&&G.ml[i][j]>=1)
        		ans.pa[G.ml[i][j]-1]=s[i]+1;
		else
			{printf("extrafG wrong when memberadsadasdaddasd\n");
			exit(1);}
		}
		else
		{
		printf("extrafG wrong because null member listasdasdas\n");
		exit(1);
		}
    }

}

void gcopy(struct graph G,struct graph *Ga)
{
    Ga->N=G.N;	Ga->com=G.com;	Ga->n=G.n;
    int i,j;


    Ga->nl=(int**)malloc(sizeof(int*)*G.com);
    Ga->wl=(int**)malloc(sizeof(int*)*G.com);
    Ga->dl=(int*)malloc(sizeof(int)*G.com);
    for (i=0;i<G.com;i++)
    {

        
        Ga->nl[i]=(int*)malloc(sizeof(int)*(G.nl[i][0]+1));
        for (j=0;j<=G.nl[i][0];j++)
        Ga->nl[i][j]=G.nl[i][j];
        
        Ga[0].wl[i]=(int*)malloc(sizeof(int)*(G.nl[i][0]+1));
        for (j=0;j<=G.nl[i][0];j++)
        Ga[0].wl[i][j]=G.wl[i][j];
  
        Ga[0].dl[i]=G.dl[i];
    }
    
    
}

void freeG(struct graph G)
{
    int i;
    for (i=0;i<G.com;i++)
    {

        free(G.nl[i]);
        free(G.wl[i]);
    }
    free(G.dl);

    free(G.nl);
    free(G.wl);
}

double movedQ(struct graph G, int child, int home, int school, int *s,int *deg)
{
	int i,j,mem;
	double ans=0;

	for (i=1;i<=G.nl[child][0];i++)
	{
		mem=s[G.nl[child][i]];
		if (mem==home)
			ans-=G.wl[child][i];
        else if (mem==school)
			ans+=G.wl[child][i];
    }
	ans=ans+G.dl[child]*1.0*(deg[home]-G.dl[child]-deg[school])/(2*G.n);

	return  (2*ans);
}

void move(struct graph *Ga,int child, int home, int school, int *s, struct graph G)
{
//first check if child is alone
	int i,j;	
	printf("in move, child=%d,home=%d,school=%d\n",child,home,school);

	if (G.ml[child][0]==Ga->ml[home][0])
	{
	printf("case1\n");
		updateG(Ga,school,home);
		s[child]=school;
		Ga->com--;
		if (home!=Ga->com)
		for (i=0;i<G.com;i++)
			if (s[i]==Ga->com)
			s[i]=home;
	}
	else
	{
	printf("case2\n");
//when home has two or more reduced nodes
		int *nl1,*nl2;
	//update ml
		nl1=(int*)malloc(sizeof(int)*(Ga->ml[home][0]-G.ml[child][0]+1));
		nl2=(int*)malloc(sizeof(int)*(Ga->ml[school][0]+G.ml[child][0]+1));
		i=1;
		nl1[0]=Ga->ml[home][0]-G.ml[child][0];
		nl2[0]=Ga->ml[school][0]+G.ml[child][0];
		while (Ga->ml[home][i]!=G.ml[child][1])
		{
			nl1[i]=Ga->ml[home][i];	i++;
		}
		while (i<=nl1[0])
		{
			nl1[i]=Ga->ml[home][i+G.ml[child][0]];	i++;
		}
		for (j=1;j<=G.ml[child][0];j++)
			nl2[j]=G.ml[child][j];
		for (i=1;i<=Ga->ml[school][0];i++)
			nl2[j+i]=Ga->ml[school][i];
		free(Ga->ml[home]);	free(Ga->ml[school]);
		Ga->ml[home]=nl1;	Ga->ml[school]=nl2;

	//update nl and wl , dl 
	Ga->dl[home]-=G.dl[child];
	Ga->dl[school]+=G.dl[child];
		int vn,sp,hp,k1,k2,temp;
		//child's neighbor, school position of home, home position of school, temp position
		sp=findxiny(school,Ga->nl[home]);
		hp=findxiny(home,Ga->nl[school]);
		
		Ga->wl[home][0]-=G.wl[child][0];
		Ga->wl[school][0]+=G.wl[child][0];
		for (i=1;i<=G.nl[child][0];i++)
		{
			vn=s[G.nl[child][i]];
			if (vn==home)
			{
				Ga->wl[home][0]-=G.wl[child][i];
				Ga->wl[home][sp]+=G.wl[child][i];
				Ga->wl[school][hp]+=G.wl[child][i];
			}
			else if (vn==school)
			{
				Ga->wl[school][0]+=G.wl[child][i];
				Ga->wl[home][sp]-=G.wl[child][i];
				Ga->wl[school][hp]-=G.wl[child][i];
			}
			else
			{
				k1=findxiny(vn,Ga->nl[home]);
				k2=findxiny(vn,Ga->nl[school]);
				if (k2<0)
				{
					nl1=(int*)malloc(sizeof(int)*(Ga->nl[school][0]+2));
					nl2=(int*)malloc(sizeof(int)*(Ga->nl[school][0]+2));
					for (j=0;j<=Ga->nl[school][0];j++)
					{
						nl1[j]=Ga->nl[school][j];	nl2[j]=Ga->wl[school][j];
					}
					nl1[j]=vn;
					nl2[j]=G.wl[child][i];
					nl1[0]++;
					free(Ga->nl[school]);	free(Ga->wl[school]);
					Ga->nl[school]=nl1;	Ga->wl[school]=nl2;
					nl1=(int*)malloc(sizeof(int)*(Ga->nl[vn][0]+2));
					nl2=(int*)malloc(sizeof(int)*(Ga->nl[vn][0]+2));
					for (j=0;j<=Ga->nl[vn][0];j++)
					{
						nl1[j]=Ga->nl[vn][j];	nl2[j]=Ga->wl[vn][j];
					}
					nl1[j]=school;
					nl2[j]=G.wl[child][i];
					nl1[0]++;
					free(Ga->nl[vn]); free(Ga->wl[vn]);
					Ga->nl[vn]=nl1;	Ga->wl[vn]=nl2;
				}
				else
				{
					Ga->wl[school][k2]+=G.wl[child][i];
					temp=findxiny(school,Ga->nl[vn]);
					Ga->wl[vn][temp]+=G.wl[child][i];
				}
				if (Ga->wl[home][k1]==G.wl[child][i])
				{
					if (k1==Ga->nl[home][0])
						Ga->nl[home][0]--;
					else
					{
						Ga->nl[home][k1]=Ga->nl[home][Ga->nl[home][0]];
						Ga->wl[home][k1]=Ga->wl[home][Ga->nl[home][0]];
						Ga->nl[home][0]--;
					}
					temp=findxiny(home,Ga->nl[vn]);
					Ga->nl[vn][temp]=Ga->nl[vn][Ga->nl[vn][0]];
					Ga->wl[vn][temp]=Ga->wl[vn][Ga->nl[vn][0]];
					Ga->nl[vn][0]--;
				}
				else
				{
					Ga->wl[home][k1]-=G.wl[child][i];
					temp=findxiny(home,Ga->nl[vn]);
					Ga->wl[vn][temp]-=G.wl[child][i];
				}

			}
		}
		
		s[child]=school;
		Ga->com--;
	}

	
	
}


struct part RG(struct graph G, int ke, unsigned int *seed)
{
    int i;
    struct part ans;
    struct graph Ga;
	
    gcopy(G,&Ga);

    ans.pa=(int*)malloc(sizeof(int)*G.N);;

    int *s,*ns,kee,*deg;
	s=(int*)malloc(sizeof(int)*G.com);
	ns=(int*)malloc(sizeof(int)*G.com);
	deg=(int*)malloc(sizeof(int)*G.com);
    double Q;
    Q=compQ(Ga);

    double dQmax,dQ;

    int j,k,g1,g2,lr,randg;

    for (i=0;i<G.com;i++)
	{
		s[i]=i; 
		ns[i]=i;
		deg[i]=G.dl[i];   
    	}

    ans.com=G.com;
    ans.Q=Q;
    int maxiter=G.com-1;
    
    for (i=1;i<=maxiter;i++)
    {
        dQmax=-2*G.n-1; lr=0;
	if (Ga.com<ke)
		kee=Ga.com;
	else
		kee=ke;
        for (j=0;j<kee;j++)
        {
   
            randg=randint(Ga.com,seed);

            for (k=1;k<=Ga.nl[randg][0];k++)
            {
                dQ=2*(Ga.wl[randg][k]-Ga.dl[randg]*1.0*Ga.dl[Ga.nl[randg][k]]/(2*Ga.n));
         
                if (abso(dQ-dQmax)<Inf_Sma)
                {
                    lr++;
                    if (rand1(seed)*lr<1)
                    {
                        g1=randg;	g2=Ga.nl[randg][k];
                    }
                }
                else if (dQ>dQmax)
                {
                    lr=1;
                    dQmax=dQ;
                    g1=randg;	g2=Ga.nl[randg][k];
                }
            }
        }

        if (lr==0) break;
        else	{
		if (Ga.nl[g1][0]<Ga.nl[g2][0])
        	{
			g1=g1+g2;	g2=g1-g2;	g1=g1-g2;
		}    	

		updateG(&Ga,g1,g2);

            Ga.com--;

		for (j=0;j<G.com;j++)
			if (s[j]==g2)
				s[j]=g1;
		if (g2!=Ga.com)
		for (j=0;j<G.com;j++)
			if (s[j]==Ga.com)
				s[j]=g2;

        }
        Q+=dQmax;

        if (Q>ans.Q)
        {

            ans.com=Ga.com;

            ans.Q=Q;
            for (j=0;j<G.com;j++)
			ns[j]=s[j];
		for (j=0;j<Ga.com;j++)
			deg[j]=Ga.dl[j];

        }

 
    }



//refine section
	int change=1,*visit;
	visit=(int*)malloc(sizeof(int)*ans.com);

	while (change)
	{
		change=0;
		for (i=0;i<G.com;i++)
		{
			g1=ns[i];
			dQmax=-1;
			for (j=0;j<ans.com;j++)
			visit[j]=0;
			for (j=1;j<=G.nl[i][0];j++)
				if (ns[G.nl[i][j]]!=g1&&visit[ns[G.nl[i][j]]]==0)
				{
					visit[ns[G.nl[i][j]]]=1;
					dQ=movedQ(G,i,g1,ns[G.nl[i][j]],ns,deg);

					if (dQ>dQmax)
					{
						dQmax=dQ;
						g2=ns[G.nl[i][j]];
					}
				}
			if (dQmax>Inf_Sma)
			{
				ans.Q+=dQmax;
				ns[i]=g2;
				deg[g1]-=G.dl[i];
				deg[g2]+=G.dl[i];

				change=1;
			}
		}
	}


	free(visit);

	extrafG(ans,G,ns);

//below is to reorder ans.partition
	j = 0;
	int *lab;
	lab=(int*)malloc(sizeof(int)*G.N);
	for (i=0;i<G.N;i++)
		lab[i]=0;
	for (i = 0;i < G.N;i++)
	{
		if (!lab[i])
		{
			j++;
			lab[i]=j;
		
		for (k = i + 1;k < G.N;k++)
			if (ans.pa[k] == ans.pa[i])
				lab[k] = lab[i];
		}
	}
	for (i=0;i<G.N;i++)
		ans.pa[i]=lab[i];
	ans.com=j;
	free(lab);
	free(s);
	free(ns);
	free(deg);
    freeG(Ga);
    
    return (ans);
}
