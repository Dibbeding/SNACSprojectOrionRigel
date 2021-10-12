
#include "rigel.h"
#include "node.h"
#include <time.h>
#include <fstream>
#include <iostream>
#include <math.h>
using std::ifstream;

#define MAXTIME 100000

#define TINY 1.0e-10
#define NMAX 500000
int DIM=7; 
int curve=-15;

double ftol = 0.00001;

float **medians = NULL;
int BETA;
int LANDMARK;
int *idarray;
int DIMENSIONS = 7;


int nodeCount = 0;
char *outputFormat = NULL;
Node *node1 = NULL;
Node *node2 = NULL;
int *first_beta=NULL;	
int *larray;
int *basic;
int firstin;
int clonei=0;
map<int,int> first_map;

FILE *runLogFP = NULL;


FILE* openOutputFile (char *prefix, char *suffix);
void openOutputFile (FILE *&fp, char *prefix, char *suffix, int stamp);
void closeOutputFile (FILE *fp);
void perSecondOutput (int stamp);


void printUsage () {
	
	cout<<"./rigel \n"
	<<"required parameters:\n"
	<<"-e:curvature number\n"
	<<"-b:neighbor number\n"
	<<"-o:outputprefix, three output:.coord & .time & .para for both landmark embedding and regular nodes embedding; .land file is the file storing landmark coordinates\n"
	<<"-l:measured distance file contain landmark# row * node number matrix\n"
	<<"-x:dimesion\n"
	<<"-t:input file prefix, .num, .ord\n"
	<<"-r:landmark id file name\n"
    <<"-u: total number of nodes in the graph\n"
    <<"-y: the start file id; if -1, mean embedding landmarks\n"
	<<"-------optional parameters\n"
    <<"-L:the total number of landmarks; By default: 100\n"
    <<"-i: the total number of the set of basic landmarks; By default: 16\n";
    exit (-1);
}



double linear_dist(double *a, double *b) {
		 int i;
		 double dist;
	
		 int norm = 2;
		 dist = 0;
	
		for (i=0; i<DIM; i++)
		{
			if (norm == 2) {
				dist += (a[i] - b[i])*(a[i] - b[i]);
				} else {
					dist += pow(fabs(a[i] - b[i]), norm);
				}
		}
	
		if (norm == 1) {
			return dist;
		} else
			return pow(dist, 1.0/norm);
	}


double arccosh(double x)
{
	double t,r;
	t=x+sqrt(x*x-1);
	r=log(t);
	return r;
}
double pb_dist(double *a, double *b)
{
	int i;
	double dist,tx=1.0,ty=1.0,tt=0.0;
	for(i=0;i<DIM;i++)
	{
		tx = tx + a[i]*a[i];
		ty = ty + b[i]*b[i];
		tt = tt + a[i]*b[i];
	}
	double t=sqrt(tx*ty)-tt;
	t=arccosh(t)*abs(curve);
	return t;
}

//Global simplex downhill algorithm
void simplex_downhill1(double **simplex, double *values, int d, double ftol,
                       double (*obj)(double *, int, int), int *num_eval,
                       int stuff) {
    
    int i, j, low, high, second_high, ssize;
    double rtol, sum, *simplex_sum, *test_point, test_value, mult;
    
  
    *num_eval = 0;
    ssize = d + 1; 
    
       simplex_sum = new double [d+1]; 
    for (i=1; i<=d; i++) {
        sum = 0.0;
        for(j=1; j<=ssize; j++) {
            sum = sum + simplex[j][i];
        }
        simplex_sum[i] = sum;
    }
    
   
    test_point = new double [d+1];
    while (1) {

        if (values[1] > values[2]) {
            low = 2;
            high = 1;
            second_high = 2;
        } else {
            low = 1;
            high = 2;
            second_high = 1;
        }
        for (i=1; i<=ssize; i++) {
            if (values[i] > values[high]) {
                second_high = high;
                high = i;
            } else if (values[i] > values[second_high] && i != high) {
                second_high = i;
            } else if (values[i] <= values[low]) {
                low = i;
            }
        }
        
        rtol=2.0*fabs(values[high]-values[low])/(fabs(values[high])+fabs(values[low])+TINY);
        if (*num_eval >= NMAX || rtol < ftol || values[low] < 1.0e-6) {
            values[1] = values[low];
            for (i=1; i<=d; i++) {
                simplex[1][i] = simplex[low][i];
            }
            break;
        }
        
        
        mult = 2.0/d;
        for (i=1; i<=d; i++) {
            test_point[i] = (simplex_sum[i] - simplex[high][i])*mult - simplex[high][i];
        }
        
        test_value = (*obj)(&test_point[1], d, stuff);
        
        (*num_eval)++;
        if (test_value < values[high]) {
         
            values[high] = test_value;
            for (i=1; i<=d; i++) {
                simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
                simplex[high][i] = test_point[i];
            }
        }
        
        if (test_value <= values[low]) {
            mult = -1.0/d;
            for (i=1; i<=d; i++) {
                test_point[i] = (simplex_sum[i] - simplex[high][i])*mult + 2.0*simplex[high][i];
            }
			
            test_value = (*obj)(&test_point[1], d, stuff);
			
            (*num_eval)++;
            if (test_value < values[high]) {
                values[high] = test_value;
                for (i=1; i<=d; i++) {
                    simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
                    simplex[high][i] = test_point[i];
                }
            }
			
        } else if (test_value >= values[second_high]) {
            mult = 0.5/d;
            for (i=1; i<=d; i++) {
                test_point[i] = (simplex_sum[i] - simplex[high][i])*mult + 0.5*simplex[high][i];
            }
			
            test_value = (*obj)(&test_point[1], d, stuff);
			
            (*num_eval)++;
            if (test_value < values[high]) {
                values[high] = test_value;
                for (i=1; i<=d; i++) {
                    simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
                    simplex[high][i] = test_point[i];
                }
            } else {
                for (i=1; i<=ssize; i++) {
					if (i != low) {
						for (j=1; j<=d; j++) {
                            simplex[i][j]=test_point[j]=0.5*(simplex[i][j]+simplex[low][j]);
						}
						values[i]=(*obj)(&test_point[1], d, stuff);
					}
				}
				*num_eval = *num_eval + d;
            
				for (i=1; i<=d; i++) { 
                    sum = 0.0;
                    for(j=1; j<=ssize; j++) {
                        sum = sum + simplex[j][i];
                    }
                    simplex_sum[i] = sum;
				}
            }
        }
    }
    
    delete [] simplex_sum;
    delete [] test_point;
    
}
//regular node embedding
void simplex_downhill2(double **simplex, double *values, int d, double ftol,
                       double (*obj)(double *, int, int), int *num_eval,
                       int stuff) {
    
    
    int i, j, low, high, second_high, ssize;
    double rtol, sum, *simplex_sum, *test_point, test_value, mult;
    
    
    *num_eval = 0;
    ssize = d + 1; 
    
   
    simplex_sum = (double*) malloc(sizeof(double)*(d+1));
    for (i=1; i<=d; i++) {
        sum = 0.0;
        for(j=1; j<=ssize; j++) {
            sum = sum + simplex[j][i];
        }
        simplex_sum[i] = sum;
    }

    test_point = (double*) malloc(sizeof(double)*(d+1)); 
    
    while (1) {
        if (values[1] > values[2]) {
            low = 2;
            high = 1;
            second_high = 2;
        } else {
            low = 1;
            high = 2;
            second_high = 1;
        }
        for (i=1; i<=ssize; i++) {
            if (values[i] > values[high]) {
				second_high = high;
				high = i;
            } else if (values[i] > values[second_high] && i != high) {
				second_high = i;
            } else if (values[i] <= values[low]) {
				low = i;
            }
        }
		
    
        rtol=2.0*fabs(values[high]-values[low])/(fabs(values[high])+fabs(values[low])+TINY);
        if (*num_eval >= NMAX || rtol < ftol || values[low] < 1.0e-6) {
            values[1] = values[low];
            for (i=1; i<=d; i++) {
				simplex[1][i] = simplex[low][i];
            }
            break;
        }
		
		
        mult = 2.0/d;
        for (i=1; i<=d; i++) {
            test_point[i] = (simplex_sum[i] - simplex[high][i])*mult - simplex[high][i];
        }
		
        test_value = (*obj)(&test_point[1], d, stuff);
		
        (*num_eval)++;
        if (test_value < values[high]) {
            values[high] = test_value;
            for (i=1; i<=d; i++) {
				simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
				simplex[high][i] = test_point[i];
            }
        }
		
        if (test_value <= values[low]) {
            mult = -1.0/d;
            for (i=1; i<=d; i++) {
				test_point[i] = (simplex_sum[i] - simplex[high][i])*mult + 2.0*simplex[high][i];
            }
            
            test_value = (*obj)(&test_point[1], d, stuff);
            
           
            (*num_eval)++;
            if (test_value < values[high]) {
				values[high] = test_value;
				for (i=1; i<=d; i++) {
					simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
					simplex[high][i] = test_point[i];
				}
            }
            
        } else if (test_value >= values[second_high]) {
            mult = 0.5/d;
            for (i=1; i<=d; i++) {
				test_point[i] = (simplex_sum[i] - simplex[high][i])*mult + 0.5*simplex[high][i];
            }
            
            test_value = (*obj)(&test_point[1], d, stuff);
            
            (*num_eval)++;
            if (test_value < values[high]) {
				values[high] = test_value;
				for (i=1; i<=d; i++) {
					simplex_sum[i] = simplex_sum[i] - simplex[high][i] + test_point[i];
					simplex[high][i] = test_point[i];
				}
            } else {
				for (i=1; i<=ssize; i++) {
                    if (i != low) {
                        for (j=1; j<=d; j++) {
                            simplex[i][j]=test_point[j]=0.5*(simplex[i][j]+simplex[low][j]);
                        }
                        values[i]=(*obj)(&test_point[1], d, stuff);
                    }
                }
                *num_eval = *num_eval + d;
                
                for (i=1; i<=d; i++) { 
                    sum = 0.0;
                    for(j=1; j<=ssize; j++) {
                        sum = sum + simplex[j][i];
                    }
                    simplex_sum[i] = sum;
                }
            }
        }
    }
    delete [] simplex_sum;
    delete [] test_point;
}


double fit_p1(double *array, int num, int backup)
{
    double dist, sum, error, actual;
    int i, j, temp;
	
    sum = 0;
    for(i=0;i<firstin;i++)
        for(j=0;j<firstin;j++)
        {
            if (i == j) continue;
            actual = medians[basic[i]][first_beta[basic[j]]];
            if (actual < 0) {
                printf("Fit_p1: No route to the host!\n");
                continue;
            }
            dist=pb_dist(&array[i*DIM],&array[j*DIM]);
            error = ((dist - actual)/actual)*((dist - actual)/actual);	 // one bug? directly sqrt...
            sum += error;
        }
    return sqrt(sum);
}



double fit_p5(double *array, int num, int test_node)
{
    double dist, sum, error, actual;
    int i, j, temp;
	
    sum = 0;
    for(i=0;i<BETA;i++)
    {
        /*if(i<clonei)
        {
            actual=1;
            temp=larray[i];
            int id=idarray[temp];
            dist = pb_dist(&array[0], &node2[id].vec[0]);
            error = ((dist - actual)/actual)*((dist - actual)/actual);
        }*/
        //else
        //{
            temp = larray[i];
            int mid=first_map[temp],nid=idarray[temp];
            actual = medians[mid][test_node];
            if (actual<0) {printf("actual is minus\n");continue;}
            dist = pb_dist(&array[0], &node2[nid].vec[0]);
            error = ((dist - actual)/actual)*((dist - actual)/actual);	
        //}
        sum += error;
    }
    return sqrt(sum);
}

FILE* openInputFile (char *prefix, char *suffix) {
	FILE *fp;
	char *filename;
	int fileLength = strlen(prefix)+strlen(suffix)+1;
	
	filename = new char [fileLength];
	memset (filename,0,fileLength);
	sprintf (filename, "%s%s", prefix, suffix);
	if ((fp = fopen(filename, "r")) == NULL) {
		printf("%s: file open error.\n", filename);
		exit (-1);
	}
	delete [] filename;
	return fp;
}

int main (int argc, char **argv) {
    char c =  0;
    char *sampleFile = NULL;
    char *trueLatencyMatrixFile = NULL;
    int seed = getpid();
    int ROUNDS = 0;
    int marks = 0;
    int txms = 0; 
    FILE *fp;
    char *landmarks= NULL;
    char *outputPrefix = NULL;
    int totalnode = 0;
    char *landcoord=NULL;
    char *para=NULL;
    char *inpre=NULL;
    DIM=0;
    BETA=0;
    //int landcom=0;
    int highest,closest,nclose,hdegree;
    int landcom = 1;
    LANDMARK=100;
    highest =100;
    firstin=16;
    hdegree=16;
    closest=16;
    clonei = 0;
    nclose = 0;
    
    nodeCount=0;
    totalnode=0;
    
    char *neims=NULL;
    int start=-1;
    int omit=0,end=0,treenum=1;
    
    while ((c = getopt(argc,argv,"b:e:s:o:x:g:m:t:w:i:r:l:u:p:y:v:n:L:")) >= 0) {
	cout<<c<<endl;
        switch (c) {
            case 'o':
                outputPrefix = optarg;
                break;
            case 'l':
                trueLatencyMatrixFile = optarg;//distance file
                break;
            case 'b':
                BETA = atoi(optarg);//neighbor number
                break;
            case 'x':
                DIM = atoi(optarg);//dimension
                break;
            case 't':
                inpre=optarg;
                break;
            case 'r':
                landmarks = optarg;//landmarks file 
                break;
            case 'u':
                nodeCount= atoi(optarg);
                totalnode = nodeCount;
                break;
            case 'e':
                curve=atoi(optarg);
                break;
            case 'y':
                start=atoi(optarg);//-1: only landmark;tree start number:using others
                treenum = start+1;
		break;
            case 'L':
		cout<<"I am here"<<endl;
                txms = atoi(optarg);
		cout<<"here3"<<endl;
                LANDMARK=txms;
		cout<<"here2"<<endl;
                highest=LANDMARK;
		cout<<LANDMARK<<endl;
		break;
            case 'i':
                firstin=atoi(optarg);
                closest = firstin;
                hdegree = firstin;
		break;
            default:
		cout<<c<<endl;
                printUsage ();
                break;
        }
	cout<<c<<endl;
    }
    cout<<"done"<<endl;
    
    if (outputPrefix == NULL) {
        printf("Output file prefix is required.\n");
        printUsage ();
    }
    if(inpre==NULL)
    {
        printf("input file prefix is required\n");
        printUsage();
    }
    if(DIM==0 || BETA==0)
    {
        printf("DIM and BETA are required\n");
        printUsage();
    }
    if(nodeCount==0)
    {
        printf("Graph node number is needed\n");
        printUsage();
    }
   
    if(landmarks == NULL || trueLatencyMatrixFile==NULL ){
		printf ("Landmark file or matrix file misses\n");
		printUsage ();
    }
    DIMENSIONS=DIM;
    if(start>treenum)
    {
        cout<<"start:"<<start<<" larger than end:"<<treenum<<endl;
        exit(1);
    }
    cout<<start<<endl;    
   	FILE *fpp=openOutputFile(outputPrefix,".para");
	
	fprintf(fpp,"landmark:%d\n",LANDMARK);
	fprintf(fpp,"highest degree node:%d\n",highest);
	fprintf(fpp,"map node:%d\n",nodeCount);
	fprintf(fpp,"basic landmark number:%d\n",firstin);
	fprintf(fpp,"highest nodes of basic:%d\n",closest);
	//fprintf(fpp,"closest number for normal:%d\n",nclose);
	fprintf(fpp,"highest nodes for normal:%d\n\n",hdegree);
	fprintf(fpp,"output prefix:%s\n",outputPrefix);
	fprintf(fpp,"matrix file:%s\n",trueLatencyMatrixFile);
	fprintf(fpp,"landmark file:%s\n",landmarks);
	fprintf(fpp,"input file prefix:%s\n",inpre);
	fprintf(fpp,"parameter name:%s\n",para);
	fprintf(fpp,"neighbour number:%d\n",BETA);
	fprintf(fpp,"demision:%d\n",DIM);
	fprintf(fpp,"landmark comptutation method:%d\n",landcom);
	fprintf(fpp,"sys node:%d\n",totalnode);
	fprintf(fpp,"curvature:%d\n",curve);
	fprintf(fpp,"start tree ID:%d\n",start);
	fprintf(fpp,"end tree ID:%d\n",treenum);
	fclose(fpp);
    int myId, yourId;
    int stamp = 0;
    float rawLatencySample;
    FILE *land;
    char filex[100];
    sprintf(filex,"%s%d.time",outputPrefix, start);
    FILE *TimeFile=fopen(filex,"w");
    time_t readstart=time(NULL);
    
     cout<<start<<" "<<omit<<" "<<end<<endl;
    //cout<<"map build\n";
    
    if (trueLatencyMatrixFile != NULL && landmarks != NULL) {
        
        //if ((land = openInputFile(landmarks,".mark"))== NULL){
        if((land = fopen (landmarks, "r")) == NULL){
            printf ("Cannot open file %s", landmarks);
            exit (-1);
        }
        
        first_beta = new int [LANDMARK];
        
        int num = 0,temp1,temp2;
        int c=0;
        
       /* while (fscanf (land, "%d %d %d\n",&temp1,&temp2,&num) >= 0) {
            if (num >= 0)
            {
                first_beta[c]=num;
                first_map[num]=c;
                c++;
            }
        }*/
        
        while (fscanf (land, "%d\n",&num) >= 0) {
            if (num >= 0)
            {
                first_beta[c]=num;
                first_map[num]=c;
                c++;
            }
        }
        fclose (land);
        printf("read landmarks\n");
        
        if ((fp = fopen (trueLatencyMatrixFile, "r")) == NULL) {
            printf ("Cannot open file %s", trueLatencyMatrixFile);
            exit (-1);
        }
        
        medians = new float* [LANDMARK];
        
        idarray = new int[nodeCount];
        for (int myIdT = 0; myIdT < LANDMARK; myIdT++) {
            medians[myIdT] = new float[nodeCount];
            for (int yourIdT = 0; yourIdT < nodeCount; yourIdT++) {
                medians[myIdT][yourIdT] = -1; /////////////// 0 -> -1
                idarray[yourIdT]=-1;
                //	 larray[myIdT][yourIdT]=-1;
            }
        }
        int idin = 0;
        while (idin<LANDMARK) {
            for(int i=0;i<nodeCount;i++)
            {fscanf(fp,"%f ",&medians[idin][i]);//fprintf(test,"%d ",int(medians[idin][i]));
            }
            fscanf(fp,"\n");
            idin++;
        }
        fclose (fp);
    }
    
    //treenum = start;
	
	FILE *NUMfile = openInputFile(inpre,".num");
	map<int,int> numMap;
	int tempnum,temptreeid;
	while(fscanf(NUMfile,"%d %d\n",&temptreeid,&tempnum)>=0)
	{
		numMap[temptreeid]=tempnum;
	}
	fclose(NUMfile);
    // FILE *rtime = openOutputFile (outputPrefix, ".time");
	
    time_t readend=time(NULL);
    /////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////
    fprintf(TimeFile,"%ld\n",(int)(readend-readstart));
    
   	{
        int ii, jj, kk;
        int count=0;
        node1 = new Node[LANDMARK];	// Basic
        //for(ii=0;ii<totalnode;ii++) {//node1_1[ii].Refresh();
        for(ii=0; ii<LANDMARK; ii++){
            node1[ii].Refresh();
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        int joinnum=0;
        if(start==-1)
        {
            
		    int a, b, i, j;
			int total_try1, total_try2;
			int num = DIM * firstin;
			double localftol = ftol;
			int num_eval;
			int stuff;
            //select the first join in basic landmark which are used into Global Simplex Downhill
			basic = new int[firstin];
			int jo=0;
			//double landtime=0, nortime=0;
            node2 = new Node[LANDMARK];
            totalnode=LANDMARK;
            cout <<"total node:"<<totalnode<<endl;;
            for(int im=0;im<LANDMARK;im++)
            {
                node2[im].Refresh();
            }
	        
			/*FILE *basf=openInputFile(landmarks,".bas");
			while(fscanf(basf,"%d\n",&basic[jo])>=0)
			{
				printf("%d\n",basic[jo]);
				jo++;
			}
			fclose(basf);*/
            for(jo = 0; jo<firstin; jo++)
            {
                basic[jo]=jo;
            }
			double min_fit = 100000000;
			double **p0;
			printf("start compute basic landmarks...\n");
			p0 = new double*[DIM*firstin+1+1];
			for(int tr=0; tr< DIM*firstin+1+1; tr++)
			{
				p0[tr] = new double [DIM*firstin+1];
			}
			double *y0 = new double[DIM*firstin+1+1];
			int innum=firstin;
			for(total_try1 = 0; total_try1 < 3; total_try1++)
			{
				//printf("try1:%d\n",total_try1);
				for(total_try2 = 0; total_try2 < 3; total_try2++)
				{
					//printf("try2:%d\n",total_try2);
					if (total_try2 == 0)
					{
						a = 1;
						for(b=1;b<=DIM*firstin;b++)
						{
							p0[a][b] = (double)(rand() % 500 - 250);
						}
						y0[a] = fit_p1(&p0[a][1], DIM*firstin, 0); //printf("%.2f ", y0[a]);
                        
						for(a=2;a<=DIM*firstin+1;a++)
						{
							for(b=1;b<=DIM*firstin;b++)
							{
								p0[a][b] = p0[1][b];
                                //				printf("a:%d;b:%d\n",a,b);
							}
							p0[a][a-1] += 2000;
                            
							y0[a] = fit_p1(&p0[a][1], DIM*firstin, 0); //printf("%.2f ", y0[a]);
						}
					}
					else
					{
						y0[1] = fit_p1(&p0[1][1], DIM*firstin, 0); //printf("%.2f ", y0[1]);
						for(a=2;a<=DIM*firstin+1;a++)
						{
							for(b=1;b<=DIM*firstin;b++)
							{
								p0[a][b] = p0[1][b];
							}
                            
							p0[a][a-1] += 2000;
							y0[a] = fit_p1(&p0[a][1], DIM*firstin, 0); //printf("%.2f ", y0[a]);
						}
					}
					simplex_downhill1(p0,y0,num,localftol,fit_p1,&j,i);
                    
					if (fit_p1(&p0[1][1], DIM*firstin, 0) < min_fit)
					{
						//printf("fit_p1 used here\n");
						min_fit = fit_p1(&p0[1][1], DIM*firstin, 0);
						//printf("%lf\n",min_fit);
						//printf("fit_p1 used here two\n");
						for(i=0;i<firstin;i++)
						{
							for(j=0;j<DIM;j++)
							{
								node1[basic[i]].vec[j]= p0[1][i*DIM+j+1];
								node2[basic[i]].vec[j]= p0[1][i*DIM+j+1];
								//printf("landmark in\n");
							}
						}
					}
					//printf("fit_p1 used \n");
				}
			}
            
            FILE *lands=openOutputFile(outputPrefix, ".land");
            for(int ki = 0; ki<firstin;ki++)
            {
                idarray[first_beta[basic[ki]]]=joinnum;
               // cout<<"Join basic:"<<basic[ki]<<" "<<first_beta[basic[ki]]<<endl;
                joinnum++;
                fprintf(lands,"%d ",first_beta[basic[ki]]);
                for(int kj =0; kj<DIM;kj++)
                    fprintf(lands,"%f ",node1[basic[ki]].vec[kj]);
                fprintf(lands,"\n");
            }
            
            for(int tr=0; tr< DIM*firstin+1+1; tr++)
            {
                delete [] p0[tr];
            }
            delete [] y0;
            delete [] p0;
            
            
			// Layer1 GNP (Pure GNP, Node4[])
            double **p;
            p = new double* [DIM+1+1];
            for(int trr=0;trr<DIM+1+1;trr++)
            {
                p[trr]=new double[DIM+1];
            }
            double *y= new double[DIM+1+1];
            larray=new int[BETA];
            //FILE *ordf=openInputFile(landmarks,".sec");
            
            map<int,int> mar;
            map<int,int>::iterator mIter;
            double *tempv;
            tempv = new double[DIM];
            
            for(int ii=0;ii<LANDMARK-firstin;ii++)
			{
				int join=0;
				int ok=0;
				//fscanf(ordf,"%d\n",&join);
                join = first_beta[firstin+ii];
				mar.clear();
				if(landcom)
				{
					int count = 0,nein=0;
					//cout<<"join:"<<join<<" "<<first_map[join]<<" "<<idarray[join]<<endl;
					if(idarray[join]>=0) {cout<<"already join in the system:"<<join<<endl;exit(1);}
					idarray[join]=joinnum;
					joinnum++;
					int neix=0;
					int rein=0;
					int jtemp=0;
					int neicount=0;
					int hig=0;
					clonei=0;
					
					while(hig<hdegree && count<BETA){
						int x=rand()%highest;
						int id=first_beta[x];
						mIter=mar.find(id);
						if(mIter==mar.end() && join!=id && idarray[id]>=0)
                        {
                            larray[count]=id;
                            mar[id]=1;
                            count++;
                            hig++;
                            //cout<<hig<<" "<<x<<endl;
                        }
						continue;
                        
					}
					//cout<<"high degree nodes\n";
					while(count<BETA)
					{
						int x=rand()%(LANDMARK-highest)+highest;
						int id=first_beta[x];
						mIter=mar.find(id);
						
						if(mIter==mar.end() && join!=id && idarray[id]>=0)
						{
							larray[count]=id;
							mar[id]=1;
							count++;
						}
						continue;
					}
					//cout<<"neigh select\n"<<endl;
				}
                
				//printf("JOIN NODE:%d\n",ii);
				{
                    
                    int num = DIM;
                    double localftol = ftol;
                    
                    //double min_fit = 10000;
                    double min_fit = 100000000;
                    map<int,int>::iterator it;
                    it=first_map.find(join);
                    int a, b;
                    int tmp_j;
                    
                    for(int total_try1 = 0; total_try1 < 3; total_try1++)
                    {
                        for(int total_try2 = 0; total_try2 < 3; total_try2++)
                        {
                            if (total_try2 == 0)
                            {
                                for(a=1;a<=DIM+1;a++)
                                {
                                    for(b=1;b<=DIM;b++)
                                    {
                                        p[a][b] = (double)(rand() % 500 - 250);
                                    }
                                    if (a>1) p[a][a-1] += 2000;
                                    if(!ok) y[a] = fit_p5(&p[a][1], DIM, join); //printf("%.2f ", y[a]);
                                    //else y[a]=fit_p4(&p[a][1],DIM,join);
                                }
                            }
                            else
                            {
                                if(!ok) y[1] = fit_p5(&p[1][1], DIM, join);
                               // else y[1]=fit_p4(&p[1][1],DIM,join);
                                for(a=2;a<=DIM+1;a++)
                                {
                                    
                                    for(b=1;b<=DIM;b++)
                                    {
                                        p[a][b] = p[1][b];
                                    }
                                    p[a][a-1] += 2000;
                                    
                                    if(!ok) y[a] = fit_p5(&p[a][1], DIM, join);
                                    //else y[a] = fit_p4(&p[a][1], DIM, join);
                                }
                            }
                            if(!ok) simplex_downhill2(p,y,num,localftol,fit_p5,&tmp_j,join);
                           // else simplex_downhill2(p,y,num,localftol,fit_p4,&tmp_j,join);
                            
                            if (!ok){
                                if(fit_p5(&p[1][1], DIM, join) < min_fit)
                                {
                                    min_fit = fit_p5(&p[1][1], DIM, join);
                                    for(int jj=0;jj<DIM;jj++)
                                    {
                                        node1[idarray[join]].vec[jj]=p[1][jj+1];
                                        node2[idarray[join]].vec[jj]=p[1][jj+1];
                                    }
									//tempv[jj]=p[1][jj+1];
                                }
                            }
                            /*else {
                             if(fit_p4(&p[1][1], DIM, join) < min_fit)
                             {
                             min_fit = fit_p4(&p[1][1], DIM, join);
                             for(int jj=0;jj<DIM;jj++)
                             node1[join].vec->v[jj] = p[1][jj+1];
                             node1[join].lastUpdateTo = 1; // !!!!!!!!!!!!
                             }
                             }*/
                            
                            //printf("fi1 end\n");
                            //cout<<"coordinate write in\n"<<endl;
                        }
                        
                    }
                    
                    fprintf(lands,"%d ",join);
                    for(int dl = 0; dl<DIM; dl++)
                    {
                        fprintf(lands,"%f ",node1[idarray[join]].vec[dl]);
                    }
                    fprintf(lands,"\n");
                    //cout<<"write into files\n";
				}
			}
            delete [] node1;
            return 1;
        }
        else{
            /*FILE *lands=fopen(landcoord,"r");
             joinnum=0;
             
             for(int i=0;i<LANDMARK;i++)
             {
             
             int id;
             fscanf(lands,"%d ",&id);
             idarray[id]=joinnum;
             joinnum++;
             if(id==3072440) cout<<"Landmark contains:3072440"<<endl;
             //int x=first_map[id];
             for(int j=0;j<DIM;j++)
             {
             double y;
             fscanf(lands,"%lf ",&y);
             node1[idarray[id]].vec->v[j]=y;
             }
             fscanf(lands,"\n");
             }
             fclose(lands);*/
            //for(int treeid=start;treeid<treenum && treeid<LANDMARK;treeid++)
            for(int treeid = start; treeid<=treenum; treeid++)
            {
                
                char outx[100];
                joinnum=0;
                sprintf(outx,"%s%d",outputPrefix,treeid);
                runLogFP=openOutputFile(outx,".coord");
                totalnode=numMap[treeid];
                cout<<"In tree "<<treeid<<" there are "<<totalnode<<endl;
                node2 = new Node[totalnode+LANDMARK];
                /*for(int im=0;im<totalnode+LANDMARK;im++)
                 {
                 if(im<LANDMARK)
                 {
                 for(int xo=0;xo<DIM;xo++)
                 {
                 node2[im].vec->v[xo]=node1[im].vec->v[xo];
                 }
                 cout<<"in landmark:"<<im<<"; "<<xo<<endl;
                 }
                 else {
                 node2[im].Refresh();
                 }
                 //cout<<im<<endl;
                 }
                 joinnum=LANDMARK;*/
               // FILE *lands=fopen(landcoord,"r");
                FILE *lands = openInputFile(outputPrefix,".land");
                for(int i=0;i<nodeCount;i++)
                {
                    idarray[i]=-1;
                }
                joinnum=0;
                
                for(int i=0;i<LANDMARK;i++)
                {
                    
                    int id;
                    fscanf(lands,"%d ",&id);
                    idarray[id]=joinnum;
                    joinnum++;
                    //if(id==3072440) cout<<"Landmark contains:3072440"<<endl;
                    //int x=first_map[id];
                    int j=0;
                    for(j=0;j<DIM;j++)
                    {
                        double y;
                        fscanf(lands,"%lf ",&y);
                        node2[idarray[id]].vec[j]=y;
                    }
                    fscanf(lands,"\n");
                    //cout<<"landmark "<<i<<" joining:"<<j<<endl;
                }
                fclose(lands);
                
                for(int i=LANDMARK;i<totalnode+LANDMARK;i++)
                {
                    node2[i].Refresh();
                }
                
                cout<<"Landmark joined into node2:"<<joinnum<<endl;
                double **p;
                p = new double* [DIM+1+1];
                for(int trr=0;trr<DIM+1+1;trr++)
                {
                    p[trr]=new double[DIM+1];
                }
                double *y= new double[DIM+1+1];
                larray=new int[BETA];
                char inpx[100];
                sprintf(inpx,"%s%d",inpre,treeid);
                FILE *ordf=openInputFile(inpx,".ord");
                
                map<int,int> mar;
                map<int,int>::iterator mIter;
                
                for(int ii=0;ii<totalnode;ii++)
                {
                    int join=0;
                    int ok=0;
                    clonei=0;
                    mar.clear();
                    {
                        int count = 0,nein=0;
                        //fscanf(ordf,"%d %d\n",&join,&nein);
                        fscanf(ordf,"%d\n",&join);
                        nein = -1;
                        if(idarray[join]>=0) {cout<<"already join in the system:"<<join<<endl;exit(1);}
                        idarray[join]=joinnum;
                        joinnum++;
                        
                        /*if(nclose==1)
						{
							mar[nein]=1;
							larray[count]=nein;
							count++;
							clonei=1;
						}*/
                        
                        
                        //cout<<"join:"<<join<<" close nei:"<<clonei<<" "<<count<<endl;
                        
                        srand(time(NULL)+join);
                        int hig=0;
                        
                        while(hig<hdegree && count<BETA){
                            int x=rand()%highest;
                            int id=first_beta[x];
                            mIter=mar.find(id);
                            
                            if(mIter==mar.end() && join!=id && idarray[id]>=0)
                            {
                                larray[count]=id;
                                mar[id]=1;
                                count++;
                                hig++;
                                //	cout<<"highest id:"<<id<<endl;;
                            }
                            
                        }
                        //cout<<"high degree nodes "<<hig<<" "<<count<<endl;;
                        while(count<BETA)
                        {
                            int x=rand()%(LANDMARK-highest)+highest;
                            int id=first_beta[x];
                            mIter=mar.find(id);
                            
                            if(mIter==mar.end() && join!=id && idarray[id]>=0)
                            {
                                larray[count]=id;
                                mar[id]=1;
                                count++;
                                cout<<"id "<<id<<endl;
                            }
                            
                        }	
                        //cout<<"neigh select "<<count<<endl;	
                        
                    }
                    
                    //printf("JOIN NODE:%d\n",ii);
                    {
                        //cout<<"start to compute "<<ii<<endl;
                        
                        int num = DIM;		
                        double localftol = ftol;
                        double min_fit = 100000000;
                        map<int,int>::iterator it;
                        it=first_map.find(join);
                        int a, b;
                        int tmp_j;
                        
                        for(int total_try1 = 0; total_try1 < 3; total_try1++)
                        {
                            for(int total_try2 = 0; total_try2 < 3; total_try2++)
                            {
                                if (total_try2 == 0)
                                {
                                    for(a=1;a<=DIM+1;a++)
                                    {
                                        for(b=1;b<=DIM;b++)
                                        {	
                                            p[a][b] = (double)(rand() % 500 - 250); 				
                                        }				
                                        if (a>1) p[a][a-1] += 2000;
                                        if(!ok) y[a] = fit_p5(&p[a][1], DIM, join); //printf("%.2f ", y[a]);
                                        //else y[a]=fit_p4(&p[a][1],DIM,join);
                                    }
                                }
                                else
                                {
                                    //cout<<total_try1<<" "<<total_try2<<"fit_p5"<<endl;
                                    if(!ok) y[1] = fit_p5(&p[1][1], DIM, join); //printf("%.2f ", y[1]);
                                    //else y[1]=fit_p4(&p[1][1],DIM,join);
                                    for(a=2;a<=DIM+1;a++)
                                    {
                                        
                                        for(b=1;b<=DIM;b++)
                                        {	
                                            p[a][b] = p[1][b];					
                                        }
                                        p[a][a-1] += 2000;
                                        
                                        if(!ok) y[a] = fit_p5(&p[a][1], DIM, join); //printf("%.2f ", y[a]);
                                       // else y[a] = fit_p4(&p[a][1], DIM, join);
                                    }
                                    //cout<<total_try1<<" "<<total_try2<<"fit_p5 ends"<<endl;
                                }
                                if(!ok) simplex_downhill2(p,y,num,localftol,fit_p5,&tmp_j,join);
                                //else simplex_downhill2(p,y,num,localftol,fit_p4,&tmp_j,join);
                                //	printf("simplex_downhill2 \n");
                                //				printf("%.2f J%d  ", fit_p4(&p[1][1], DIM, ii), tmp_j);
                                //cout<<"simlex downhill"<<total_try1<<' '<<total_try2<<endl;
                                if (!ok){
                                    if(fit_p5(&p[1][1], DIM, join) < min_fit)
                                    {
                                        min_fit = fit_p5(&p[1][1], DIM, join);
                                        for(int jj=0;jj<DIM;jj++)
                                        {
                                            node2[idarray[join]].vec[jj] = p[1][jj+1];
                                            //tempv[jj]=p[1][jj+1];
                                        }
                                    }
                                }
                            }
                            
                        }
                        fprintf(runLogFP,"%d ",join);
                        int tr=idarray[join];
                        for(int dl = 0; dl<DIM; dl++)
                        {
                            fprintf(runLogFP,"%f ",node2[tr].vec[dl]);
                        }
                        fprintf(runLogFP,"\n");
                        
                    }
                }
                fclose(runLogFP);
                delete [] node2;
                delete [] larray;
                
                for(int trr=0;trr<DIM+1+1;trr++)
                {
                    delete [] p[trr];
                }
                delete [] p;
                delete [] y;
                //delete [] larray;
                mar.clear();
            }
        }
    }
    /*time_t endcom=time(NULL);
    fprintf(TimeFile,"%ld\n",(int)(endcom-readend));
    
    fclose(TimeFile);
    if (medians != NULL) {
        for (int myId = 0; myId < LANDMARK; myId++) {
            delete [] medians[myId];
        }
        delete [] medians;
        //delete [] larray;
    }*/
    return 0;
}


FILE* openOutputFile (char *prefix, char *suffix) {
  FILE *fp;
  char *filename;
  int fileLength = strlen(prefix)+strlen(suffix)+1;

  filename = new char [fileLength];
  memset (filename,0,fileLength);
  sprintf (filename, "%s%s", prefix, suffix);
  if ((fp = fopen(filename, "w")) == NULL) {
    printf("%s: file open error.\n", filename);
    exit (-1);
  }
  delete [] filename;
  return fp;
}

void openOutputFile (FILE *&fp, char *prefix, char *suffix, int round) {
  char *filename;
  int fileLength = strlen(prefix)+strlen(suffix)+32;

  filename = new char [fileLength];
  memset (filename,0,fileLength);
  sprintf (filename, "%s%05d%s", prefix, round, suffix);
  if ((fp = fopen(filename, "w")) == NULL) {
    printf("%s: file open error.\n", filename);
    exit (-1);
  }
  delete [] filename;
}

void closeOutputFile (FILE *fp) {
  fclose (fp);
  fp = NULL;
}




