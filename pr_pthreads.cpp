

#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include "iostream"
#include <vector>
#include<sys/time.h>
#include <math.h>
#include<string.h>
#define DAMPING 0.85
#define epi 0.00001
using namespace std;



int numnodes, edge;
double n = 1.0;
vector< vector <int> > matvec;
vector<int> row;
//Array for storing outlinkss
vector<int> outlinks;
vector<double> pr;
vector<double> pr_old;
pthread_barrier_t barr;
pthread_barrier_t barrier;	
pthread_mutex_t mutex;


//#define threads 16
int* rank;
int threads = 16;
//int iterations = 0;
//double sum;
//int rowcnt=0, colcnt=0;
//int x;
double get_time()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	double d = t.tv_sec + (double) t.tv_usec/1000000;
	return d;
}
void* pagerank(void* arg)
{	
	int x=0;
	double sum;
	int rank = *((int *)arg);
	int i=0;
	vector<int>::iterator j;

	int start_idx = numnodes%threads>=rank ? ((numnodes/threads)+1)*rank : (numnodes/threads)*rank+numnodes%threads;
	int rows = numnodes%threads>rank ? (numnodes/threads)+1 : numnodes/threads;
	int end_idx = start_idx+rows;
	double tranv[rows+1];
	if(rank==0){
		for(i=0;i<=numnodes;i++){
		}}
	while(n > epi){

	for(i=start_idx+1;i<end_idx+1;i++){
		sum = 0.0;
		vector<int> k= matvec.at(i);
		for(j=(k).begin();j!=(k).end();j++)
		{
			sum += DAMPING*(pr_old[*j]/outlinks[*j]);
		}
		tranv[i-start_idx]=sum + (1-DAMPING)/(double) numnodes;
	}
	int rc = pthread_barrier_wait(&barr);
//	printf("Thread %d Called Barrier\n", rank);
	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
	{
		printf("Could not wait on barrier\n");
		exit(-1);
	}
//	printf("Thread %d Passed Barrier\n", rank);
	for(i=0;i<rows;i++){
		pr[start_idx+i+1]=tranv[i+1];
	}
   pthread_mutex_lock(&mutex);

//	if(rank==0){
		n = 0;
		for(x=1;x<=numnodes;x++)
		{
			n += (pr[x] - pr_old[x])*(pr[x] - pr_old[x]);
		}
		n = sqrt(n);
	//	printf("Thread %d Norm value: %lf\n",rank, n);
      pthread_mutex_unlock(&mutex);
//	}
	rc = pthread_barrier_wait(&barrier);
//	printf("Thread %d Called Barrier\n", rank);
//	if(rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD)
//	{
//		printf("Could not wait on barrier\n");
//		exit(-1);
//	}
//	printf("Thread %d Passed Barrier\n", rank);
	for(i=0;i<rows;i++){
		pr_old[start_idx+i+1]=tranv[i+1];
	}
	for(x=0;x<=numnodes;x++)
		printf("%.16f\n", pr[x]);
	}
	pthread_exit(NULL);
}

int main(int argc, char *argv[])
{
	FILE *file;	
//	printf("\nHello");
	if(argc != 3)
	{
		printf("\nUsage: Enter the name of the file containing the Page-Ranks and number of threads. \n");
		exit(EXIT_FAILURE);
	}
	//printf("\nThe name of the file entered is %s", (char *) argv[2]);
	file = fopen(argv[1], "r");
	if(file==NULL)
	{
		printf("Unable to open file\n");
		exit(EXIT_FAILURE);
	}
	char line[512];
	
	threads = atoi(argv[2]);
	
//	printf("%s\n", argv[2]);
	fgets(line, sizeof(line), file);
	fgets(line, sizeof(line), file);
	//printf("%s\n", line);
	fgets(line, sizeof(line), file);
	//getting the number of nodes and edges
	sscanf(line, "# Nodes: %d Edges: %d", &numnodes, &edge);
	int x;
	for(x=0;x<=numnodes;x++)
	{
		pr.push_back(1/(double) numnodes);
		pr_old.push_back(1/(double) numnodes);
		outlinks.push_back(0);
	}

	for(int it=0;it<numnodes+1;it++)
		matvec.push_back(row);
	int from, to;
	fgets(line, sizeof(line), file);
	while(1)
	{
		fgets(line, sizeof(line), file);		
		//printf("\nLine read is %s", line);
		sscanf(line, "%d	%d", &from, &to);
		//printf("From %d To %d\n", from, to);
		if(feof(file))
			break;
		outlinks[from] = outlinks[from] + 1;
		matvec[to].push_back(from);
		if(feof(file))
			break;
	}

	printf("\nNumber of nodes is %d and edges is %d", numnodes, edge);
	printf("\n");
	//Done reading the file into the necessary data structures

	//Now, partition the data into various sets. 

	pthread_t thr[threads];
	if(pthread_barrier_init(&barr, NULL, threads))
	{
		printf("Could not create a barrier\n");
		return -1;
	}
	if(pthread_barrier_init(&barrier, NULL, threads))
	{
		printf("Could not create a barrier\n");
		return -1;
	}
	
    if(pthread_mutex_init(&mutex, NULL))
    {
        printf("Unable to initialize a mutex\n");
        return -1;
    }

	double time_start = get_time();
	rank = (int*) malloc(sizeof(int)*threads);
	for(int i = 0; i < threads; i++)
	{
		rank[i]=0;
//		printf("Creating thread:%d\n", i);
		rank[i] = i;
		if(pthread_create(&thr[i], NULL, &pagerank, (void*)(&rank[i])))
		{
			printf("Could not create thread %d\n", i);
			return -1;
		}
	}

	for(int i = 0; i < threads; ++i)
	{
		if(pthread_join(thr[i], NULL))
		{
			printf("Could not join thread %d\n", i);
			return -1;
		}
	}

//	pagerank();
	double time_end = get_time();

	printf("The final pageranks are as follows\n");
	for(x=1;x<=numnodes;x++)
		printf("%.16f\n", pr[x]);
//	printf("\nNumber of iterations is %d", iterations);
	printf("\nThe time taken for the execution is %lf", time_end - time_start);
	printf("\nEnd of pagerank!\n");
}



