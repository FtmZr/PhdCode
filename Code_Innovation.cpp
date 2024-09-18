////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fatemeh Zarei- Gent, Belgium- 03.07.2024
//
// This C++ program simulates the dynamic evolution of agents' technological complexity
// and diversity in a social network. Each agent interacts with its social contacts
// for social learning or self-creates technologies, updating its set of technologies accordingly.
// The program reads a social network, executes the simulation for a specified number of time steps, 
// and records various measures, such as mean complexity, diversity, and similarities among agents
// at different time steps. Results are saved to multiple CSV files for analysis.
////////////////////////////////////////////////////////////////////////////////////////////////////////

//// START OF THE CODE

// CALL LIBRARIES

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
#include <iostream> 
#include <cstdlib>

// SPECIFY MODEL PARAMETERS

// N= the number of agents in the network
int N=1000; 
//c1= parameter beta in sigmoid function
//c= complexity of agents
//d= diversity of agents
int c1=30, c, d;
//a= parameter alpha in the sigmoid function
float a=1.2 , th=0.1;
//T= total number of time steps
int T= 4*pow(10,6);
using namespace std;


// CALL INPUT DATA
// The input data is a network
// 1000 row. Each row has an agent number in the first column and its social contact number in the next column. 
ifstream in("ERListOfNeighbors1.csv");

// CREATE FILES FOR SAVING THE OUTPUT
 
ofstream file1("ERmean-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the average complexity of agents at each time step 
ofstream file2("ERfinal-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the complexity of every agent at the final time step
ofstream file3("ERTime-1-2(c1=100 , a=1.2, th=0.2).csv");  ///  the time steps
ofstream file4("ERsimilarity-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the similarity of agents at each time step 
ofstream file5("ERCounter-1-2(c1=100 , a=1.2, th=0.2).csv");  // the number of creation and the number of social learning at each time step
ofstream file6("ERmidTimeI-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the complecity of every agent at t=T/4
ofstream file7("ERmidTimeJ-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the complecity of every agent at t=T/2
ofstream file8("ERmidTimek-1-2(c1=100 , a=1.2, th=0.2).csv"); /// the complecity of every agent at t=3T/4

// DEFINE FUNCTIONS TO BE USED IN THE MAIN PROGRAM

// This function reads the degree of agent i	
int Deg(int i, vector< vector<string> > array) { 
    string s;
    int x;
    s=array[i][0] ;
    stringstream geek(s);
    geek >> x; 
   return x;
}

// This function reads K_th social contact of agent i
int neib(int i, vector< vector<string> > array,int k) { 
    int x;
    string s;
	s=array[i][k+1] ;
    stringstream geek(s);
    geek >> x;    
  return x;
}

// This function checks whether the specific technology x is in the set of skills "list"
bool Find(std::list<std::string> ilist, std::string x) {
    bool status = false;
    
    std::list<std::string>::iterator it;
    // Fetch the iterator of element with value 'the'
    it = std::find(ilist.begin(), ilist.end(), x);
    // Check if iterator points to end or not
    if (it != ilist.end()) {
        // It does not point to end, it means element exists in list
        status = true;
    }
    
    return status;
}

// This function finds the technologies in the set of technologies list a  but not in the set of technologies list b
list<string> uncommon(list<string> a, list<string> b){   
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

// This function finds the technologies which are in both set of technologies list a and set of technologies list b
list<string> common(list<string> a, list<string> b){   
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

// This function finds the n_th element of list a 
string element(list<string> a, int n){ 

	// Initialize iterator to list 
    list<string>::iterator it = a.begin();
    // Move the iterator by 1 elements 
    advance(it, n); 
return(*it)	  ;	
}

// This function finds the most complex technology in the set of technologies v
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

// This function finds the most complete technologies in a set of technology list
// a which is not in the set of technologies list b//find most complex in "b" which is not in node "a"  
string mostComplexWord(list<string> a, list<string> b){ 
   
    std::vector<int> Cb;
    int L=b.size();
    list<string>::iterator it=b.begin();
    string x;
    for(int i=0; i<L;i++){
    	x=(*it);
  	    Cb.push_back(x.length());
  	    advance(it, 1);
  	}
    list<int>  argMax = FindMaxElements(Cb);
    
    list<string>  b_Complexes;
    list <int> :: iterator it1; 
    for(it1 = argMax.begin(); it1 != argMax.end(); ++it1){
  	    b_Complexes.push_back(element(b,*it1));
	} 
    list<string> b_choise= uncommon(b_Complexes,a);
    int L2= b_choise.size();
    string s;
    if (L2==0){
    	s="";
	}
	else{
////////////////////////////////
/////// PREVIOUS INCORRECT CODE
/////// int r = (rand() % L2);
////////////////////////////////
/////// CORRECTED CODE
		int r = (int)(L2 * (double)(rand())/RAND_MAX );  
////////////////////////////////

		s=element(b_choise,r);
	}
    
    return(s);
} 

// This function creates a new technology "s" to be assigned to a set of technologies of an agent "a"
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
    
////////////////////////////////
/// PREVIOUS INCORRECT CODE
/// int r = (rand() % L);
////////////////////////////////
/// CORRECTED CODE
    int r = (int)(L * (double)(rand())/RAND_MAX );
////////////////////////////////
     
    list<int>::iterator it1 = argMax.begin(); 
    advance(it1, r);
    int n=*it1;
    
    it=a.begin();
    advance(it, n);
    string s = *it;
    string A= "A";
    string B= "B";
    
////////////////////////////////
/// PREVIOUS INCORRECT CODE
/// r = (rand() % 100);
////////////////////////////////
/// CORRECTED CODE    
    r = (int)(100 * (double)(rand())/RAND_MAX ); 
////////////////////////////////
    
    if(r<50){
    	s=s+A; 
	}
	else{
		s=s+B; 
	}
    
    return(s);
    
}
    
// This function finds the similarity os set of technologies a and set of technologies b    
float similarity(list<string> a, list<string> b){
	
	list<string> c =common(a,b);
	a.merge(b);
   
	int s1= c.size();
	int s2= a.size();
	
	return(s1*1./(s2-s1));
	
}

// This function provides a random technology of A and B to use for an agent at time 0
string AB(){
	string s;
	
////////////////////////////////
/// PREVIOUS INCORRECT CODE
/// int r = (rand() % 100);
////////////////////////////////
/// CORRECTED CODE
	int r = (int)(100 * (double)(rand())/RAND_MAX );  
////////////////////////////////

    if(r<50){
    	s="A"; 
	}
	
	else{
		s="B"; 
	}
	return(s);
	
}

///==== MAIN PROGRAM ====///
int main(void)
{
	time_t start, end; 
	time(&start); 
	
	// Intialize the random number generator
    srand(time(0));
	
// create a 2-dimensional array from the called data
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
 
// save the title of the file
  
    file1 << "meanC" <<","<<"meanD"<< endl ;  //  average complexity and average diversity at each time step
	file2 << "finalC" <<","<<"finalD"<< endl ; //  the complexity of all agents at final time step
	file6 << "mid1C" <<","<<"mid1D"<< endl ; // the complexity of all agents at final t= T/4
	file7 << "mid2C" <<","<<"mid2D"<< endl ; // the complexity of all agents at final t= T/2
	file8 << "mid3C" <<","<<"mid3D"<< endl ; // the complexity of all agents at final t= 3T/4

// define the set of technology, social contacts, complexity, and diversity for each agent in the network
struct Network{
	std::list<string> listOfWords;
	std::vector<int> Neighbors;
	int Complexity;
	int Diversity;

};

Network G[N];
list<string>::iterator it;

// The loop makes the initial state for all agents of the network  

for( int i=0; i<N; i++){
	
	//list Of technologies
   string s= AB();
   (G[i].listOfWords).push_back(s);
   
   //number of social contacts
	int len = Deg(i, array);
    for( int j=0; j<len; j++){
      	int ni= neib(i, array, j);
        (G[i].Neighbors).push_back(ni);
	}
	  
	G[i].Complexity=1;
	G[i].Diversity=1;
}

list<string> Li , Lj;
string technology;
vector<int> Ni;

int SumC=N;
int SumD=N;

int SumI=0, SumL=0;

//  This is the evolution of the system - update the state of the agents
for( int t=0; t<T; t++){
	
    if(t% (T/100)==0){
		file1 << (SumC*1./N) <<","<<(SumD*1./N)<< endl ;
                file5<< SumI <<","<<SumL<< endl ;
	}

////////////////////////////////
/// PREVIOUS INCORRECT CODE
/// int i = (rand() % N);
////////////////////////////////
/// CORRECTED CODE
	int i = (int)(N * (double)(rand())/RAND_MAX ); 
////////////////////////////////

	Li=G[i].listOfWords;
	c = G[i].Complexity;
    d = G[i].Diversity;

	double p =  1 / ( 1 + exp(-c1*(d-a*c) ));
	
////////////////////////////////
/// PREVIOUS INCORRECT CODE
/// int Ac =1000; 
/// float r = (rand() % Ac )*1./Ac;
////////////////////////////////
/// CORRECTED CODE
	double r = rand()/(double)RAND_MAX;  
////////////////////////////////
	
	// creation of technology	
	if(r<p){
		
		technology = Innovation(Li);
		(G[i].listOfWords).push_back(technology);
		
		G[i].Complexity++;
		SumC++;
		
    	G[i].Diversity++;
    	SumD++;
    	
    	SumI++; //Counter of creation
    		
	}
	// social learning
	else{
		
		Ni=G[i].Neighbors;
		int Deg=Ni.size();
		
////////////////////////////////
/////// PREVIOUS INCORRECT CODE
/////// int r = (rand() % Deg);
////////////////////////////////
/////// CORRECTED CODE
        int r = (int)(Deg * (double)(rand())/RAND_MAX );
////////////////////////////////

		int j= Ni[r];
	
		Lj=G[j].listOfWords;
		float s= similarity(Li,Lj);
		if(s>th){
			
            //copy most complex technology j to i
	        if(mostComplexWord(Li,Lj).length()>0){
	         	technology =mostComplexWord(Li,Lj);
	       	    if (technology.length()<=G[i].Complexity){
     	            (G[i].listOfWords).push_back(technology);
     	            
     	        	G[i].Diversity++;
     	        	SumD++;
     	        	
     	        	SumL++; //counter of social learning
	
		       }	
         	}
         	
         	//copy most complex technology i to j
	        if(mostComplexWord(Lj,Li).length()>0){
	         	technology =mostComplexWord(Lj,Li);
	       	    if (technology.length()<=G[j].Complexity){
     	        	(G[j].listOfWords).push_back(technology);
     	        	
     	        	G[j].Diversity++;
     	        	SumD++;
     	        	
     	        	SumL++; //counter of social learning
		       }	
         	}
         	
		}	
	}
	
	// SAVE THE RESULTS - OUTPUT 
   
    // saving the complexity of all agents at t=T/4
    if(t==int (T/4)){
   	for (int i=0; i<N; i++)
	    {
	        file6<< G[i].Complexity <<","<<G[i].Diversity<< endl ;
	    }
    }
    // saving the complexity of all agents at t=T/2
	else if(t==int (T/2)){
   	for (int i=0; i<N; i++)
	    {
	        file7<< G[i].Complexity <<","<<G[i].Diversity<< endl ;
	    }
	}
	else if(t==int ((3*T)/4)){
    // saving the complexity of all agents at t=3T/4
   	for (int i=0; i<N; i++)
	    {
	        file8<< G[i].Complexity <<","<<G[i].Diversity<< endl ;
	    }
	}
}

    // save the complexity of all agents at final time step

	for (int i=0; i<N; i++)
	{
        file2<< G[i].Complexity <<","<<G[i].Diversity<< endl ;
	}
	        
	// save the similarity of all agents at the final time step	
		file4 << "similarity" << endl ;
		float s;
	    for( int i=0; i<N; i++){
	    	Li=G[i].listOfWords;
	    	s=0;
	    	Ni=G[i].Neighbors;
		    int Deg=Ni.size();
		    for( int r=0; r<Deg; r++){
		    	int j= Ni[r];
		    	Lj=G[j].listOfWords;
		    	s= s+similarity(Li,Lj);
		    	
			}
			s=s/Deg;
			file4<< s <<endl ;
		}
		
	time(&end);
	double time_taken = double(end - start);
	file3<< time_taken;

    return 0;
}

//// END OF THE CODE
