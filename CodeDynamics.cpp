#include <string>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector>
#include <list>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>     
#include <math.h>      
#include <utility> 
#include <stdexcept> 
#include <cmath>
#include <random>

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include <algorithm>
#include <stdlib.h>
#include<time.h>
#include <stdio.h>      
#include <math.h>       
#include <string>
#include <fstream>
#include <vector>
#include <utility> 
#include <stdexcept> 
#include <sstream>





int N=1000; 
int intial_k=5000; 
int stabilization= N/2;
double mu= 0.5; 
double nu=5; 
double nk=10; 
double epsilon=0.0001; 
int ui=0;  
int T=  pow (10,7), Tmin=pow(10,4);
int nNet= 2; 
int dist=1; 
int nd=20; 

using namespace std;





////////////////////

std:: string Net(int i){
   string s ;
   if (i==2){
   	s="ER";
   }
   return s;
}



int source(int i, vector< vector<string> > array) {
    string s;
    int x;
    s=array[i][0] ;
    stringstream geek(s);
    geek >> x; 
   return x;
}



  
int target(int i, vector< vector<string> > array) {
    string s;
    int x;
    s=array[i][1] ;
    stringstream geek(s);
    geek >> x; 
   return x;
}




double PowerLawDistribution(double xmin, double xmax, double alpha) {
        double r = (rand() / (double)RAND_MAX);
        double l = (pow(xmax, alpha + 1.0) - pow(xmin, alpha + 1.0)) * r + pow(xmin, alpha + 1.0);
        return pow(l, 1.0 / (alpha + 1.0));
    }
    
    

    double ExponentialCutoff(double mean) {
        double lambda = 1.0 / mean; 
        double r = (rand() / (double)RAND_MAX);
        return -1.0 / lambda * log(1.0 - r);
    }

double Distribution(int i){
   double number;
   
   
   
   int max=10000000;
   
   if (i==1){ 
    
    number=-1;
    while (number <0 or number>max){
    	std::random_device rd;
        std::mt19937 gen(rd());
        std::exponential_distribution<> d(1./(38));
         number=d(gen);
    
	}
    
   
   }
   
   if (i==2){ 
    number=-1;
    while (number <0 or number>max){
    	std::random_device rd;
   
    std::mt19937 gen(rd());
    std::lognormal_distribution <double> d(0,2.7);
    number=d(gen);
   
	}
   }
   
   if (i==3){ 
     number=-1;
  

    
    double xmin = 0.00045, xmax = 100000, alpha = -1.5,  desired_mean = 1000.0, upper_bound =1000.0;
   
    
    double random_number = PowerLawDistribution(xmin, xmax, alpha);
    random_number *= ExponentialCutoff(desired_mean);
        if (random_number > upper_bound) {
            random_number = upper_bound;
        }
    number=random_number;
   
   }
   
   
   return number;
}



double normDist() { 
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double number = dis(gen); 
    return number;
}





string DistributionName(int i){
   string s;
   
   if (i==1){ 
     s="(exponential)";	
   }
   
   if (i==2){ 
    s="(lognormal)";
   }
   
    if (i==3){ 
    s="(Powerlaw)";
   }
   
   
   return s;
}




struct Qline{
    int edge;  
    double active; 
   };




 bool compare(Qline q1, Qline q2)
   {
    return (q1.active < q2.active);
   }
  


int main(void)
{	
time_t start, end; 
time(&start); 


 string sk3= Net(nNet)+DistributionName(dist)+"time.csv";
 ofstream file3(sk3);
 
 
 //string sk4= 'ER'+DistributionName(dist)+"time.csv";
 //ofstream file4(sk4);
 
 
 string sk2= Net(nNet)+DistributionName(dist)+"mu"+to_string(mu).substr(0, 4)+ "dRange.csv";
 ofstream file2(sk2);
 
 file2<<"d"<<"," <<"stabilizationTime"<<","<<"finalNumber"<<"\n";

 for ( int z=4; z<=nd; z++){
   	float  d=(0.5*z)/nd;
   	cout<<"d="<<"\t"<<d<<"\n";


	////////////////////////////////////////
	
	

int stabilize=0;
int clusterNumber=0;

   
   int nNet=2;
   
   for ( int u=1+ui; u<=nu+ui; u++){



   	cout<<"\n"<< u<<"\n";
   	string su= Net(nNet)+"edges"+to_string(u)+".csv";
 
   	ifstream in(su);
   	
   	string line, field;
   	
    vector< vector<string> > array;  
    vector<string> v;           
         
    
    int E=0;
    while ( getline(in,line) )    
    {
        v.clear();
        stringstream ss(line);

        while (getline(ss,field,','))  
        {
            v.push_back(field);  
        }

        array.push_back(v);  
        E= E+1;
    }
   	
   	
   	
  
   
   for( int k=1; k<=nk; k++){





    srand(time(0));



    struct Network{
	double Opinion;
   };	
   
   
   


   Network G[N];
  
    for( int i=0; i<N; i++){
    	double r = normDist() ;  
        G[i].Opinion =r;
	  
    }
    
    
    
   std::vector<Qline> Q; 

   for (int i=0; i<E; ++i) {
   	double number = Distribution(dist);
   	Q.push_back({i,number});
   	
   }

    
    sort(Q.begin(), Q.end(), compare);
    
   int nNotFixed=0, nFixed=0 ; 
    




    for( int t=0; t<intial_k; t++){ 
         double number = Distribution(dist);
		 Q[0].active=Q[0].active+number;

		Qline qj=Q[0],qi=Q[1];
		int q=0;
		while(compare(qj, qi)==0 & q<E-1){
			 Q[q]=qi;
	          q=q+1;
	          qi=Q[q+1];
	
			}
		Q[q]=qj;
     }
     
  
  
  
  
   
    
   int change=0;
   int j; 
   for( int t=0; t<T; t++){ 
   
     
     
     
		j= Q[0].edge;  
		 
		double number = Distribution(dist);
		Q[0].active=Q[0].active+number; 
		
	
		Qline qj=Q[0],qi=Q[1];
		int q=0;
		while(compare(qj, qi)==0 & q<E-1 ){
			 Q[q]=qi;
	          q=q+1;
	          qi=Q[q+1];
	
			}
			
			
		Q[q]=qj;
		
	//	
		
	    
	    
	    

   	    int n1 = source(j, array);
   	    int n2 = target(j, array);
   	  
   	  
   	    double x1=G[n1].Opinion;
   	    double x2=G[n2].Opinion; 

   	  
   	  
   	  
		if (abs(x1-x2)<d & abs(x1-x2)>epsilon){
			double dx1= mu*(x2-x1);
			double dx2= mu*(x1-x2);
			G[n1].Opinion=x1+dx1;
			G[n2].Opinion=x2+dx2;
			nNotFixed=nNotFixed+1;
		}
		
		else{
			if(nNotFixed==0){  
				nFixed=nFixed+1;
			}
			else{
				nFixed=1;
				nNotFixed=0;	 
			}					
	    }
	
	
	
		
	

	
	

            if(t>Tmin){
        	if (nFixed== stabilization){ 
        		if (t>stabilize){ 
        		   stabilize=t;
				}
				      		
			break;
	    	}
	    }
			
    } 
    
    
    
    
    
    list <double> final;
    for( int i=0; i<N; i++){ 
    	final.push_back(G[i].Opinion);
	}
	final.sort();
	
	
	
	
	double min= *final.begin();


	int counterCluster=0;
	
  
    int sum=0;  
		for (auto it = final.begin(); it != final.end(); ++it){
			if((*it-min)<=epsilon){
				sum=sum+1;
			}
			else{
			
				if(sum>(N/100)){ 
				   counterCluster=counterCluster+1; 
				}
				sum=0;
			}
			min=*it;	
	      }
	      
	      if(sum==N){
	      	counterCluster=1; 
		  }
    
    
    

    if(counterCluster>clusterNumber){
		clusterNumber=counterCluster;
	}
 
 
 

	 file2<<d<<"," <<stabilize<<"," <<clusterNumber<<"\n";
	
	}   
}   


    
  	
}
    ////////////////
    time(&end);
    double time_taken = double(end - start);
	file3<<"runeTime"<<"," <<time_taken<<"\n";
	



	
    return 0;

}

