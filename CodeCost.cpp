#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include <algorithm>
#include <stdlib.h>
#include<time.h>
#include <stdio.h>      /* printf */
#include <math.h>       /* exp */
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

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


double alpha=0.1  ;
double beta= 1-alpha;
double U0=0.5;

int N=1000;
int Edge=0;
int c1=30, c2=100, c, d;
double a=1.2 ;
int T= 11*pow (10,4) ; 
int lagT= N ; 
int m=5;  // number of nejborhs that we check. 
int is=4, js=1, i, j; 
//string s= Net +"ListOfNeighbors"+to_string(i).substr(0, 1)+".csv";

 
double epsilon=0.0001; // minium  probability to colaborate
using namespace std;
//////////////////// calling data

int listSize;
int limitSize=20; ///Memory
string Net= "net1";
string s= "ListOfNeighbors"+to_string(is).substr(0, 1)+".csv";

ifstream in(s);

string ss= "dm"+to_string(is).substr(0, 1)+".csv";
ifstream inn(ss);

 /////////////////saving datastring s;
string s1= Net +"mean-"+to_string(is).substr(0, 1)+"-"+to_string(js).substr(0, 1)+"-"+to_string(U0).substr(0, 3)+"a"+to_string(alpha).substr(0, 4)+"b"+to_string(::beta).substr(0, 4)+".csv";
ofstream file1(s1);
string s2= Net +"final-"+to_string(is).substr(0, 1)+"-"+to_string(js).substr(0, 1)+"-"+to_string(U0).substr(0, 3)+"a"+to_string(alpha).substr(0, 4)+"b"+to_string(::beta).substr(0, 4)+".csv";
ofstream file2(s2);
string s3= Net +"Time-"+to_string(is).substr(0, 1)+"-"+to_string(js).substr(0, 1)+"-"+to_string(U0).substr(0, 3)+"a"+to_string(alpha).substr(0, 4)+"b"+to_string(::beta).substr(0, 4)+".csv";
ofstream file3(s3);
string s5= Net +"Counter-"+to_string(is).substr(0, 1)+"-"+to_string(js).substr(0, 1)+"-"+to_string(U0).substr(0, 3)+"a"+to_string(alpha).substr(0, 4)+"b"+to_string(::beta).substr(0, 4)+".csv";
ofstream file5(s5);
string s55= Net +"Similarity-"+to_string(is).substr(0, 1)+"-"+to_string(js).substr(0, 1)+"-"+to_string(U0).substr(0, 3)+"a"+to_string(alpha).substr(0, 4)+"b"+to_string(::beta).substr(0, 4)+".csv";
ofstream file55(s55);


 int count_bigger(const std::vector<double>& elems) {
 	
    return std::count_if(elems.begin(), elems.end(), [](double c){return c > epsilon;});
      } 
	
int Deg(int i, vector< vector<string> > array) { // Degree of nod i
    string s;
    int x;
    s=array[i][0] ;
    stringstream geek(s);
    geek >> x; 
   return x;
}
int neib(int i, vector< vector<string> > array,int k) { // K th neibor of of nod i
    int x;
    string s;
	s=array[i][k+1] ;
    stringstream geek(s);
    geek >> x;    
  return x;
}

double dis(int i, vector< vector<string> > array2,int k) { // dis of i and k
    double x;
    string s;
	s=array2[i][k] ;
    stringstream geek(s);
    geek >> x;    
  return x;
}



double cost(double d) { // dis of i and k
    double x;
    x= alpha *3.0 *sqrt(d)/4.0;  
  return x;
}


//edit later
double benefit(double s) { // dis of i and k
    double x;
    double ps;
    if(s==0 or s==1){
    	ps=0;
	}
	else{
		ps= -s * log(s) - (1.0 - s) * log(1.0 - s);
	}
    
    x= ::beta* ps;  
  return x;
}

double Utility(double d, double s){ // dis of i and k
    double x;
    x= benefit(s) - cost(d);
  return x;
}


bool Find(std::list<std::string> ilist, std::string x) {
    bool status = false;
    
    std::list<std::string>::iterator it;
    it = std::find(ilist.begin(), ilist.end(), x);
    if (it != ilist.end()) {
        status = true;
    }
    
    return status;
}


list<string> uncommon(list<string> a, list<string> b){   //They are only in a not b
	list<string> uc;
	int L = a.size();
	// Initialize iterator to list 
    list<string>::iterator it = a.begin();
	string x;
	for(int i=0; i<L; ++i){
		if(Find( b, *it)== 0){
            uc.push_back(*it);
		}		
		// Move the iterator by 1 elements 
        advance(it, 1); 
	}
    return uc;
} 

list<string> common(list<string> a, list<string> b){   //They are only in a not b
	list<string> uc;
	int L = a.size();
	// Initialize iterator to list 
    list<string>::iterator it = a.begin();
	string x;
	for(int i=0; i<L; ++i){
		if(Find( b, *it)== 1){
            uc.push_back(*it);
		}		
		// Move the iterator by 1 elements 
        advance(it, 1); 
	}
    return uc;
} 

string element(list<string> a, int n){ 

	// Initialize iterator to list 
    list<string>::iterator it = a.begin();
    // Move the iterator by 1 elements 
    advance(it, n); 
return(*it)	  ;	
}


list<int> FindMaxElements(std::vector<int> v)
{
	
	std::list<int> indices;
    int current_max = v[0];

    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (v[i] > current_max)
        {
            current_max = v[i];
            indices.clear();
        }

        if (v[i] == current_max)
        {
            indices.push_back(i);
        }
    }

    return  indices;
}




int FindMaxElementLenght(std::vector<int> v)
{
	

    int current_max = v[0];

    for (std::size_t i = 0; i < v.size(); ++i)
    {
        if (v[i] > current_max)
        {
            current_max = v[i];
        }
        
    }

    return  current_max;
}

   
   
string Word(list<string> a, list<string> b){ //find word in "b" which is not in node a
   
    std::vector<int> Ca;
    int L=a.size();
    list<string>::iterator it=a.begin();
    string x;
    for(int i=0; i<L;i++){
    	x=(*it);
  	    Ca.push_back(x.length());
  	    advance(it, 1);
  	}
  	//chose the ones which has complexity less that Ca
    int  MaxCa = FindMaxElementLenght(Ca);
    
    
    it=b.begin();
  

    list<string>  b_List;
    for(it = b.begin(); it != b.end(); ++it){
    	x=(*it);
  	    if(x.length()<=MaxCa){
  	    	b_List.push_back(x);
		  }
  	}
    
   
    list<string> b_choise= uncommon(b_List,a);
    int L2= b_choise.size();
    string s;
    if (L2==0){
    	s="";
	}
	else{
		int r = (int)(L2 * (double)(rand())/RAND_MAX );
		s=element(b_choise,r);
	}
    
    return(s);
} 

  
   
string mostComplexWord(list<string> a, list<string> b){ //find most complwx in "b" which is not in node a

   std::vector<int> Ca;
    int L=a.size();
    list<string>::iterator it=a.begin();
    string x;
    for(int i=0; i<L;i++){
    	x=(*it);
  	    Ca.push_back(x.length());
  	    advance(it, 1);
  	}
  	//chose the ones which has complexity less that Ca
    int  MaxCa = FindMaxElementLenght(Ca);
    
    
    it=b.begin();
  

    list<string>  b_List;
    for(it = b.begin(); it != b.end(); ++it){
    	x=(*it);
  	    if(x.length()<=MaxCa){
  	    	b_List.push_back(x);
		  }
  	}
    
   
    list<string> b_choise= uncommon(b_List,a);
    
   //////////////
    std::vector<int> Cb;
    L=b_choise.size();
    it=b_choise.begin();
    
    for(int i=0; i<L;i++){
    	x=(*it);
  	    Cb.push_back(x.length());
  	    advance(it, 1);
  	}
  	
  	
    list<int>  argMax = FindMaxElements(Cb);
    
    list<string>  b_Complexes;
    list <int> :: iterator it1; 
    for(it1 = argMax.begin(); it1 != argMax.end(); ++it1){
  	    b_Complexes.push_back(element(b_choise,*it1));
	} 
   
    int L2= b_Complexes.size();
    string s;
    if (L2==0){
    	s="";
	}
	else{
		int r = (int)(L2 * (double)(rand())/RAND_MAX );
		s=element(b_Complexes,r);
	}
    
    return(s);
} 


string Innovation(list<string> a){
	
	std::vector<int> Ca;
	list<string>::iterator it=a.begin();
	string x;
	for(int i=0; i<a.size(); i++){
        x= *it;
		Ca.push_back(x.length());
		advance(it, 1);
	}
	
	
    list<int> argMax = FindMaxElements (Ca);
    int L= argMax.size();
    int r = (int)(L * (double)(rand())/RAND_MAX );
    list<int>::iterator it1 = argMax.begin(); 
    advance(it1, r);
    int n=*it1;
    
    it=a.begin();
    advance(it, n);
    string s = *it;
    string A= "A";
    string B= "B";
    
    
    r = (int)(100 * (double)(rand())/RAND_MAX );
    
    if(r<50){
    	s=s+A; 
	}
	else{
		s=s+B; 
	}
    
   
    
    return(s);
    
}
    
double similarity(list<string> a, list<string> b){
	
	list<string> c =common(a,b);
	a.merge(b);
   
	int s1= c.size();
	int s2= a.size();
	
	return(s1*1./(s2-s1));
	
}

string AB(){
	string s;
	
	
	int r = (int)(100 * (double)(rand())/RAND_MAX );
    if(r<50){
    	s="A"; 
	}
	else{
		s="B"; 
	}
	return(s);
	
}
	
struct Network{
	std::list<string> listOfWords;
    //std::vector<int> Neighbors;
	int Complexity;
	int Diversity;
	
	int locX;
	int locY;

};	

double distance(int xi,int xj, int yi, int yj){
	double x= xi - xj;
	double y= yi - yj;
	double z= sqrt(x*x+y*y);
	return(z);
}

   
double prob(double u){
	
	  double z =  1 / ( 1 + exp(c2*(U0-u)));
	  return z;
}



int Distribution(std::vector<double> Psd){//"exponential"
   double number;
   std::random_device rd;
   std::mt19937 gen(rd());
   std::discrete_distribution<> d(Psd.begin(), Psd.end());
   
   number=d(gen);


   return number;
}
    
    
int main(void)
{
	
	
	
	time_t start, end; 
	time(&start); 
	
	/* Intializes random number generator */
    srand(time(0));
    

    
    	
//////////////////// calling data (list of neibor)
    
    string line, field;
    vector< vector<string> > array;  // the 2D array
    vector<string> v;                // array of values for one line only

    while ( getline(in,line) )    // get next line in file
    {
        v.clear();
        stringstream ss(line);

        while (getline(ss,field,','))  // break line into comma delimitted fields
        {
            v.push_back(field);  // add each field to the 1D array
        }

        array.push_back(v);  // add the 1D array to the 2D array
    }
  //////////////////////////////////////////////////////////////////
  
  
  
  //////////////////// calling data (dm)
    
    string line2, field2;
    vector< vector<string> > array2;  // the 2D array
    vector<string> v2;                // array of values for one line only

    while ( getline(inn,line2) )    // get next line in file
    {
        v2.clear();
        stringstream ss(line2);

        while (getline(ss,field2,','))  // break line into comma delimitted fields
        {
            v2.push_back(field2);  // add each field to the 1D array
        }

        array2.push_back(v2);  // add the 1D array to the 2D array
    }
  //////////////////////////////////////////////////////////////////
    
  
  
   /////////////////saving data
   
    file1 << "meanC" <<","<<"meanD"<< ","<<"meanC2" <<","<<"meanD2"<< ","<<"meanC3" <<","<<"meanD3"<<","<<"meanC4" <<","<<"meanD4"<<endl ;
	file2 << "finalC" <<","<<"finalD"<< endl ;


	
  //////////////////////////////////////////////////////////////////
 
    
	
struct Network{
	std::list<string> listOfWords;
    std::vector<int> Neighbors;
    std::vector<double> DisNeighbors;  //distance 
    std::vector<double> PrNeighbors;  //probabiliy
	int Complexity;
	int Diversity;

};	


Network G[N];
list<string>::iterator it;



	
int r;
// The loop:


for( int i=0; i<N; i++){
	
	//listOfWords
   string s= AB();
   (G[i].listOfWords).push_back(s);
   	G[i].Complexity=1;
	G[i].Diversity=1;
	
	 //Neighbors
   
      int len = Deg(i, array);
    
      for( int j=0; j<len; j++){
      	int ni= neib(i, array, j);
        (G[i].Neighbors).push_back(ni);
        
        (G[i].DisNeighbors).push_back(dis(i, array2, ni));
        (G[i].PrNeighbors).push_back(0);
        
           
	  }
   
  
   
	
}



int Ac =1000; //Accuracy of random generator

list<string> Li , Lj;
string technology;
vector<int> Ni;





	

file55 <<"t"<<","<< "Similarity "<<"," <<"Utility"<< endl;





int SumC=N;
int SumD=N;

int SumC2=N, SumD2=N;  // sigma N(1^2)=N
int SumC3=N, SumD3=N;  // sigma N(1^2)=N
int SumC4=N, SumD4=N;  // sigma N(1^2)=N

int SumI=0, SumL=0, SumSuccess=0;
double sumSimilarity = 0.0;
double avgSimilarity=0.0;
double sumUtility = 0.0;
double avgUtility =0;

for( int t=0; t<T; t++){
	
	if(t% (T/10)==0){
		cout<<t/(T/10)<<"\n";
	}
	
	
	
	//update similariry
	            

                
				if(t% lagT==0){
                         cout<< t<< end;
					
						
						//probability  (not fixed)
						Edge=0;
						sumSimilarity = 0.0;
						sumUtility=0.0;
						for( int i=0; i<N; i++){
							Ni=G[i].Neighbors;
							Li=G[i].listOfWords;
							double sumMp=0;
							int Deg=Ni.size();
							Edge= Edge+ Deg;
							
							for( int jj=0; jj<Deg; jj++){
                            	j=Ni[jj];
                            	Lj=G[j].listOfWords;
                            	double s= similarity(Li,Lj);
                            	sumSimilarity += s;
								double dstnc=  G[i].DisNeighbors[jj];
								double u= Utility(dstnc, s);
								sumUtility += u;
								double prb = prob(u);
								G[i].PrNeighbors[jj] = prb;
								sumMp=sumMp+prb;
								//std::cout<<prb<<"\t";	
							}
                            
                            //avgSimilarity = sumSimilarity / Deg;	
                           //file55 << avgSimilarity << ",";
							
							
							//sumMp=sumMp-Mp[i][i];	
							//Mp[i][i]=0;
							
							/*for( int jj=0; jj<Deg; jj++){
                            	j=Ni[jj];
								G[i].PrNeighbors[jj]= G[i].PrNeighbors[jj]/sumMp;
									
									
							}	*/				
			        	}  
			        	
			        
	        	        avgSimilarity= sumSimilarity /Edge;
	        	        avgUtility= sumUtility/Edge;
	        	        
                        file55 <<t<<","<< avgSimilarity <<"," <<avgUtility<< endl; 
                      	  
		//finish updating

            	}
	
    
    if(t% (T/100)==0){
		file1 << (SumC*1./N) <<","<<(SumD*1./N)<<","<< (SumC2*1./N) <<","<<(SumD2*1./N)<<","<< (SumC3*1./N) <<","<<(SumD3*1./N)<<","<< (SumC4*1./N) <<","<<(SumD4*1./N)<< endl ;
        file5<< SumI <<","<<SumL<< endl ;
      
	}

	
	
	int i = (int)(N * (double)(rand())/RAND_MAX );
	Li=G[i].listOfWords;
	c = G[i].Complexity;
    d = G[i].Diversity;
	
	double p =  1 / ( 1 + exp(-c1*(d-a*c) ));

	double r = rand()/(double)RAND_MAX;
	
	if(r<p){//Inovation
	//	std::cout<<"yesI"<<"\n";
		technology = Innovation(Li);
		(G[i].listOfWords).push_back(technology);
		
		
		int moment1= G[i].Complexity;
		G[i].Complexity++;
		int moment2= G[i].Complexity;
		SumC++;
		
		SumC2= SumC2-pow(moment1,2)+pow(moment2,2);
		SumC3= SumC3-pow(moment1,3)+pow(moment2,3);
		SumC4= SumC4-pow(moment1,4)+pow(moment2,4);
		
		
    	G[i].Diversity++;
    	
    	
    	SumI++; //counter of innovation
    	
		//cout<<"inovation"<<"|"<< technology<<"\n";	
		//limit size
				
		listSize=G[i].Diversity;
		if (listSize>limitSize){	
			it=G[i].listOfWords.begin();
		 //because we don't want to remove the last thech 
			int r = (int)((listSize-2) * (double)(rand())/RAND_MAX );
			std::advance(it, r);
			G[i].listOfWords.remove(*it);
			G[i].Diversity=listSize-1;
	        }
			else{
			SumD++;		
			}  
		 
   	 
	}
	else{	//learnig
		//std::cout<<"yesL"<<"\n";
		std::vector<double>  Psd;
		
		
		
		Ni=G[i].Neighbors;
		Psd=G[i].PrNeighbors;
		
		if (count_bigger (Psd)>0){
			//cout<<"count>0"<<"\n";
    	int jj= Distribution(Psd);
        j=Ni[jj];
        int Deg=Ni.size();
        
        ///if(t<50)
        //std::cout<<i<<"\t"<<Deg<<"\t"<<jj<<"\t"<<j<<"\n";



	//

	
//	if(r<pl){
	
     

		Lj=G[j].listOfWords;
		

            //copymostComplexWord j to i
	        if(Word(Li,Lj).length()>0){
	        //	std::cout<<"yesL"<<"\n";
	        	
			
	         	technology =Word(Li,Lj);
	       	    
     	        (G[i].listOfWords).push_back(technology);
     	            
     	        G[i].Diversity++;
     	        
     	       	SumL++; //counter of social learning
     	        //cout<<"i"<<"|"<< technology<<"\n";	
				 //limit size
				listSize=G[i].Diversity;
		        if (listSize>limitSize){	
				   	it=G[i].listOfWords.begin();
				  //because we don't want to remove the last thech
				   	int r = (int)((listSize-2) * (double)(rand())/RAND_MAX );
				   	std::advance(it, r);
				    G[i].listOfWords.remove(*it);
				    G[i].Diversity=listSize-1;
	             }
				 else{
				 SumD++;		
				 }   	
         	}
         	
         	//copymostComplexWord i to j
	        if(Word(Lj,Li).length()>0){
	        //	std::cout<<"yesL"<<"\n";
	         	technology =Word(Lj,Li);
     	        (G[j].listOfWords).push_back(technology);
     	        G[j].Diversity++;
     	        
     	       	SumL++; //counter od social learning
     	        //	cout<<"j"<<"|"<< technology<<"\n";
     	        //limit size
				listSize=G[j].Diversity;
		        if (listSize>limitSize){	
				   	it=G[j].listOfWords.begin();
				   	int r = (int)((listSize-2) * (double)(rand())/RAND_MAX ); //because we don't want to remove the last thech
				   	std::advance(it, r);
				    G[j].listOfWords.remove(*it);
				    G[j].Diversity=listSize-1;
	             }
				 else{
				 SumD++;		
				 }  
     	        
     	        
		       	
         	}
         	
	//	}
//	}
    }
    
    //else{
    //	cout<<"coun=0"<<"\n";
//	}
	}
	
    
}
    
    
   
    
   
    /////

	    for (int i=0; i<N; i++)
	    {
	        file2<< G[i].Complexity <<","<<G[i].Diversity<< endl ;
	    }
	    
	
	    
		

		
	time(&end);
	double time_taken = double(end - start);
	file3<< time_taken;

   
	 
    return 0;
}

