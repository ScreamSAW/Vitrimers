#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<string.h>
#include<cmath>
#include<algorithm>
#define N 1650
#define S 1100
using namespace std;
string s;
int guide[N+20];
float map[S+20][N+20][5];

void read()
{
	freopen("xyz.xyz","r",stdin);
	getline(cin,s);
	getline(cin,s);
	for(int i=1;i<=N;i++)
	{
		char p;
		while(1)
		{
			scanf("%c",&p);
			if(p==' ')
			continue;
			else
			break;
		}
		//printf("%d %c\n",i,p);
		switch(p)
		{
			case 'A':
				guide[i]=2;
				break;
			case 'B':
				guide[i]=3;
				break;
			case 'C':
				guide[i]=4;
				break;
			case 'N':
				guide[i]=1;
				break;
			case 'X':
				guide[i]=5;
				break;
			default:
				break;
		}
		scanf("%f%f%f",&map[0][i][0],&map[0][i][1],&map[0][i][2]);
		getline(cin,s);
	}
	for(int i=1;i<=S;i++)
	{
		getline(cin,s);
		getline(cin,s);
		for(int j=1;j<=N;j++)
		{
			char p;
			while(1)
			{
				scanf("%c",&p);
				if(p==' ')
				continue;
				else
				break;
			}
			scanf("%f%f%f",&map[i][j][0],&map[i][j][1],&map[i][j][2]);
			getline(cin,s);
		}
	}
	//printf("%d %f %f %f",guide[N],map[S][N][0],map[S][N][1],map[S][N][2]);
	fclose(stdin);
	return;
}

float dist(int id1,int id2,int step)
{
	float ans;
	ans=pow((map[step][id1][0]-map[step][id2][0]),2)+pow((map[step][id1][1]-map[step][id2][1]),2)+pow((map[step][id1][2]-map[step][id2][2]),2);
	ans=pow(ans,0.5);
	return ans;
}
//C-B-rdf
int listC[N+20];
int listB[N+20];
void xlinkdegree()
{
	freopen("xld.out","w",stdout);
	int numC=0,numB=0;
	for(int i=1;i<=N;i++)
	{
		if(guide[i]==4)
		{
			numC++;
			listC[numC]=i;
		}
		else if(guide[i]==3)
		{
			numB++;
			listB[numB]=i;
		}
		else 
			continue;
	}
	

	fclose(stdout);
	return;
}
int main()
{
	read();
	
	return 0;
}
