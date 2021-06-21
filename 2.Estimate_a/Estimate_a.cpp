#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <set>			//std::set
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <sstream>
#include <utility>
using namespace std;
#include "Randnum_MT.h"


long int N;         //Node number 
#define pi          3.141592653589793238462643383279502884197   
double rand_num;
MTRand Rnd(1997);// open a class with automatically generate seed(1977)
//////////////////////////////////////////////////////////////////////////
vector<vector<long int> > edge_list;
vector<double> theta;
vector<double> kappa;

void Read_meank_in_each_year(string filename, vector<double> &degree) //should check the file name.
{
	//No.              Gc_N             ave_k.Gc         ave_c.Gc  
	string line;
	ifstream file(filename.c_str());
	while (getline(file, line))
	{
		if (line[0] != '#' ) //without comment
		{
			std::istringstream iss(line);
			long int tw;
			double N,K,C;
			while(iss >> tw >> N>> K>> C)
			{
				degree.push_back(K);
			}
		}
	}
	file.close();
	cout<<"read mean degree success!"<<endl;
}
vector<vector<long int> > Read_Graph(string filename) //should check the file name.
{
	long int i,j;
	//filename+="_edgelist.txt";
	ifstream infile(filename.c_str()); 
	N=-1;
	while (infile >> i)
	{
		if(i>N)
		{
			N=i;
		}
	}
	N=N+1;//max node number// id+1
	vector<vector<long int> > adj_list(N);	
	infile.clear();
    infile.seekg(0, std::ios_base::beg);
	while (infile >> i >> j)
	{
		if ( (i!=j) && (find(adj_list[i].begin(), adj_list[i].end(), j) == adj_list[i].end()) )// can't self and multi-connection!
		{
			adj_list[i].push_back(j);	
			adj_list[j].push_back(i);	
		}
	}
	infile.close();	
	cout<<"input adj file is success!"<<endl;
	cout<<"total number:"<<N<<endl;
	return adj_list;
}
void Read_parameters(string filename, long int &realN, double &beta, double & mu) //should check the file name.
{
	ifstream infile(filename.c_str());	
	infile>>realN>>beta>>mu;	
	cout<<"read Nodes Bet Mu success!"<<endl;
	cout<<"read Nodes "<<realN<<endl;
	cout<<"beta "<<beta<<endl;
	cout<<"mu "<<mu<<endl;
	infile.close();

}
void Read_kappa_theta(string filename, vector<double> &kappa, vector<double> &theta, long int Nn) //should check the file name.
{
	// should include ....#include <sstream> #include <string> 
	//initialization kappa, theta;
	vector<double> kappa0(Nn,0.0);
	vector<double> theta0(Nn,0.0);
	kappa.swap(kappa0);
	theta.swap(theta0);
	
	//filename+="_coordinates.txt";
	string line;
	ifstream file(filename.c_str());
	while (getline(file, line))
	{
		if (line[0] != '#' ) //without comment
		{
			std::istringstream iss(line);
			long int i;
			double kappai,thetai;
			while(iss >> i >> kappai>> thetai)
			{
				kappa[i]=kappai;
				theta[i]=thetai;
			}
		}
	}
	file.close();
	cout<<"read kappa theta success!"<<endl;
}
vector<double> Read_kappa_from_stable_dist(string filename, vector<long int> &node_in_which_supernode, vector<vector<long int> > &supernode_contain_which_node,long int &N, long int &Nl1) //should check the file name.
{
	long int i,j;
	filename+=".txt";
	ifstream infile(filename.c_str()); 	
	long int node;
	long int supernode;
	double kappai;
	//find the number of subnode and supernode	
	N=-1;
	Nl1=-1;	
	while(infile >> node>> supernode>> kappai)			
	{
		if(node>Nl1)
		{
			Nl1=node;
		}
		if(supernode>N)
		{
			N=supernode;
		}
	}
	N=N+1;//max node number// id+1
	Nl1=Nl1+1;	
	//cout<<"total number of supernodes:"<<N<<endl;
	//cout<<"total number of nodes in uplayer:"<<Nl1<<endl;
	
	
	//node i in wich supernode 
	//initialization node_in_which_supernode(N, -1);
	vector<long int> tempvector(Nl1, -1);
	node_in_which_supernode.swap(tempvector);
	//initialization supernode_contain_which_node(Nn)
	vector<vector<long int> > tvector(N); //supernode[i] contain sub node j1--j2--j3......
	supernode_contain_which_node.swap(tvector);
	
	vector<double> kappa_l1(Nl1,0.0);
	
	infile.clear();
    infile.seekg(0, std::ios_base::beg);
	while(infile >> node>> supernode>> kappai)	
	{		
		node_in_which_supernode[node]=supernode;
		supernode_contain_which_node[supernode].push_back(node);
		kappa_l1[node]=kappai;
	}
	infile.close();		
	return kappa_l1;
}
void sort_index(vector<vector<long int> > &edge_list, vector<double> &kappa, vector<double> &theta, long int &N) //re-asign the id for each node
{
    // Declaring vector of pairs
    vector< pair <double, long int> > vect;
 
    // Initializing 1st and 2nd element of pairs of theta
    // Entering values in vector of pairs
    for (long int i=0; i<N; i++)
	{
		if(edge_list[i].size()>0)
		{
			vect.push_back( make_pair(theta[i], i) );
		}
	}  

    // Using simple sort() function to sort
    sort(vect.begin(), vect.end());
 
     // Printing the sorted vector(after using sort())
    vector<double> new_theta;
	vector<double> new_kappa;
	vector<long int> new_id(N,0);
    for (long int i=0; i<vect.size(); i++)
    {
        // vect[i].first-->theta, vect[i].second-->old index
        // 1st and 2nd element of pair respectively
		new_theta.push_back(vect[i].first);
		new_kappa.push_back( kappa[vect[i].second] );
        new_id[vect[i].second]=i;
    } 
	new_theta.swap(theta);
	new_kappa.swap(kappa);
	vector<vector<long int> > new_edge_list(vect.size());
	for(long int i=0;i<N;i++)
	{
		for(long int j=0;j<edge_list[i].size();j++)
		{
			long int a=new_id[i];
			long int b=new_id[edge_list[i][j]];
			new_edge_list[a].push_back(b);
		}
	
	}
	new_edge_list.swap(edge_list);
	N=vect.size();
	//cout<<"N is"<<N<<endl;
   // return new_id;
}
int Real_N(vector<vector<long int> > edge_list, long int Nn)
{
	
	long int realN=0;	
	for(long int i=0;i<Nn;i++)
	{		
		if (edge_list[i].size() > 0)
		{
			realN++;
		}
		
	}
	return realN;
}
double conncted_probability(double theta1,double theta2,double kappa1,double kappa2,long int Nn, double beta, double mu)
{
	double delta_theta=pi-abs(pi-abs(theta1-theta2));			
	double d=(double)Nn/(2*pi)*delta_theta;
	double chi=1.0/(1.0+pow(d/(mu*kappa1*kappa2), beta) );	
	return chi;
}
vector<double> cal_connect_pattern(vector<double> p)
{
	long int pattern_size=(int)(p.size());
	vector<double> P(pattern_size,0.0);	
    double prb=1.0;
	for(long int i=0;i<p.size();i++)
	{
		prb*=(1-p[i]);
	}
	for(long int i=0;i<p.size();i++)
	{
		rand_num=Rnd.randExc();
		if(rand_num<p[i]/(1-prb))
		{
			P[i]=1;
			for(long int j=i+1;j<p.size();j++)
			{
				double randx;
				randx=Rnd.randExc();
				if(randx<p[j])
				{
					P[j]=1;			
				}

			}
			break;
		}
		else
		{
			prb=prb/(1-p[i]);
		}
	}
	/*for(long int i=1;i<P.size();i++)
	{
		cout<<P[i];	
	}
	cout<<endl;*/
	return P;
}
double correct_theta(double xtheta)
{
	while(xtheta>2*pi)
	{
		xtheta=xtheta-2*pi;	
	}
	while(xtheta<0)
	{
		xtheta=xtheta+2*pi;	
	}
	return xtheta;
}
double meank_maxk(vector<vector<long int> > edge_list, long int Nn)
{
	double meank=0;
	long int maxk=0;
	double sumn = 0;
	for(long int i=0;i<Nn;i++)
	{
		if(edge_list[i].size()>maxk)
			maxk=edge_list[i].size();
		if (edge_list[i].size() > 0)
			sumn++;
		meank+=edge_list[i].size();
	}
	meank /= (double)(Nn);
	//cout<<"maxk="<<maxk<<" "<<setw(10)<<"mean_degree="<<meank<<endl;
	return meank;
}

vector<vector<long int> > Growth(vector<vector<long int> > edge_list, long int Nn, long int Nl1, vector<double> &kappa, vector<double> &theta, vector<double> kappa_l1,\
		vector<long int> &node_in_which_supernode, vector<vector<long int> > &supernode_contain_which_node, double r, double beta,double mu, double &mu_new)
{	
	vector<double> theta_l1(Nl1, 0.0);	
	vector<vector<long int> > aij(Nl1);  //aij edge list in layer 1	
	//****************************************************************************************	
	//asign theta to sub nodes 	
	for(long int i=0;i<Nn;i++)
	{	
		if(supernode_contain_which_node[i].size()==1) //i do not split
		{
			long int j=supernode_contain_which_node[i][0];
			theta_l1[j]=theta[i];
		}
		else //split 2
		{
			long int j1,j2;
			//make sure j1 in the left, and j2 in the right of i
			if(supernode_contain_which_node[i][0]<supernode_contain_which_node[i][1])
			{
				j1=supernode_contain_which_node[i][0];
				j2=supernode_contain_which_node[i][1];
			}
			else
			{
				j1=supernode_contain_which_node[i][1];
				j2=supernode_contain_which_node[i][0];
			}			
			//find j1's theta			
			double dtheta;
			double dx_super;
			dtheta=2*pi/(double)(Nl1);
			if(i==0)
			{
				dx_super=pi-abs(pi-abs(theta[i]-theta[Nn-1]));
			}
			else
			{
				dx_super=pi-abs(pi-abs(theta[i]-theta[i-1]));
			}
			if(dtheta>dx_super/2.0)
			{
				dtheta=dx_super/2.0;			
			}			
			theta_l1[j1]=theta[i]-dtheta*Rnd.randDblExc();		//rand_num between (theta[i]-dtheta, theta[i])
			theta_l1[j1]=correct_theta(theta_l1[j1]);
			// find j2's theta
			dtheta=2*pi/(double)(Nl1);
			if(i==Nn-1)
			{
				dx_super=pi-abs(pi-abs(theta[0]-theta[i]));
			}
			else
			{
				dx_super=pi-abs(pi-abs(theta[i]-theta[i+1]));
			}
			if(dtheta>dx_super/2.0)
			{
				dtheta=dx_super/2.0;			
			}			
			theta_l1[j2]=theta[i]+dtheta*Rnd.randDblExc();		//rand_num between (theta[i], dtheta+theta[i])
			theta_l1[j2]=correct_theta(theta_l1[j2]);			
		}
	}
	//cout<<"finished theta part!!"<<endl;getchar();
	
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	// connectied edges between sub nodes
	//update mu
	mu=r*mu;	//mu do not return
	for(long int i=0;i<Nl1;i++)
	{
		aij[i].clear();
	}
	long int supernode1, supernode2;
	//connected sub nodes in the same supernode 
	for(supernode1=0;supernode1<Nn;supernode1++)
	{	
		if(supernode_contain_which_node[supernode1].size()>1)
		{
			long int a1, a2;
			a1=supernode_contain_which_node[supernode1][0];
			a2=supernode_contain_which_node[supernode1][1];

			double chi=conncted_probability(theta_l1[a1],theta_l1[a2],kappa_l1[a1],kappa_l1[a2], Nl1, beta, mu);	
			rand_num=Rnd.randExc();  // real number in [0,1)
			if(rand_num<chi && (find(aij[a1].begin(), aij[a1].end(), a2) == aij[a1].end()))
			{
				aij[a1].push_back(a2);
				aij[a2].push_back(a1);
			}
		}
	}	
	//connected sub nodes between two super nodes 
	for(supernode1=0;supernode1<Nn;supernode1++)
	{		
		for(long int k=0;k<edge_list[supernode1].size();k++)
		{
			supernode2=edge_list[supernode1][k];
			if(supernode2>supernode1)// up triangle
			{
				long int i,j;
				vector <double> pr;
				pr.clear();
				long int two_points_in_a_link[4][2];
				long int m=0;				
				for(long int e1=0;e1<supernode_contain_which_node[supernode1].size();e1++)
				{	
					i=supernode_contain_which_node[supernode1][e1];
					for(long int e2=0;e2<supernode_contain_which_node[supernode2].size();e2++)
					{						
						j=supernode_contain_which_node[supernode2][e2];										
						double chi=conncted_probability(theta_l1[i],theta_l1[j],kappa_l1[i],kappa_l1[j], Nl1, beta, mu);	
						pr.push_back(chi);
						two_points_in_a_link[m][0]=i;
						two_points_in_a_link[m][1]=j;
						m++;
					}
				}
				//cal the connection pattern between four nodes.
				vector<double> Pn;
				Pn=cal_connect_pattern(pr);
				for(long int n=0;n<Pn.size();n++)
				{					
					if(Pn[n]==1)
					{
						long int i=two_points_in_a_link[n][0];
						long int j=two_points_in_a_link[n][1];
						if( find(aij[i].begin(), aij[i].end(), j) == aij[i].end() )
						{
							aij[i].push_back(j);
							aij[j].push_back(i);							
						}
					}
				}


			}	
		}
	}
    
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//reture kappa, theta and node_in_which_supernode in layer 1;
	kappa.swap(kappa_l1);
	theta.swap(theta_l1);
	mu_new=mu;
	return aij;	
	
}
vector<double > cluster_coefficence(vector<vector<long int> > edge_list, long int Nn, double &meanC) //The clustering coefficient of each node
{
	long int i,j,k,u,v;
	vector<double > Cluster_coeff;
	double C,ci;	
	C=0.0;
	double sum=0; //sum is the number of degree larger than 1
	for(i=0;i<Nn;i++)
	{
		ci=0.0;
		if(edge_list[i].size()>1)//degree larger than 1
		{
			for(j=0;j<edge_list[i].size();j++)
			{
				for(k=j+1;k<edge_list[i].size();k++)
				{
					u=edge_list[i][j];
					v=edge_list[i][k];
					if(find(edge_list[u].begin(),edge_list[u].end(),v) !=edge_list[u].end())
					{
						ci+=1.0;
					}
				}
			}		
			ci=2.0*ci/(edge_list[i].size()* (edge_list[i].size()-1.0));
			sum++;
		}
		else
			ci=0.0;
		Cluster_coeff.push_back(ci);
		C+=ci;
	}
	//C/=(double)(Nn);
	C/=sum; //sum is the number of degree larger than 1
	meanC=C;
	//cout<<"the mean clustering coefficient is:"<<" "<<setw(10)<<meanC<<endl;
	return Cluster_coeff;
}
void output_edgelist(string filename, vector<vector<long int> > edge_list,long int Nn)
{
	string file1,file2;
	file1=filename+"_edgelist.txt";	
	file2=filename+"_coordinates.txt";	
	ofstream outfile10(file1.c_str(),ios::out); 	
	ofstream outfile11(file2.c_str(),ios::out); 
	
	for(long int i=0;i<Nn;i++)
	{
		for(long int k=0;k<edge_list[i].size();k++)
		{
			if(edge_list[i][k]>i)
				outfile10<<setprecision(10)<<left<<setw(10)<<i<<" "<<setw(10)<<edge_list[i][k]<<endl;		
		}
		if(edge_list[i].size()>0)
		outfile11<<setprecision(10)<<left<<setw(10)<<i<<" "<<setw(16)<<kappa[i]<<" "<<setw(16)<<theta[i]<<endl;
	}
	outfile10.close();
	outfile11.close();
}
void output_supernode(string filename, vector<long int>  node_in_which_supernode)
{
	string file1;
	file1=filename+"_supernode.txt";		 
	ofstream outfile10(file1.c_str(),ios::out); 	
	 
	for(long int i=0;i<node_in_which_supernode.size();i++)
	{		
		outfile10<<left<<" "<<setw(10)<<i<<" "<<setw(10)<<node_in_which_supernode[i]<<endl;
	}
	outfile10.close();	
}

void Cal_mean_std(vector<double> v, double &ave, double &std)
{
	ave=0;
	for(int i=0;i<v.size();i++)
	{
		ave+=v[i]/(double)(v.size());	
	}
	std=0;
	for(int i=0;i<v.size();i++)
	{
		std+=pow((v[i]-ave),2.0);	
	}
	std=sqrt(std);
}	
int main()
{
	double time_begin=time(NULL);
	//input the slope ex of k as a function of N in log-log scale, std_es is the standard deviation of fitting slope.
	double ex=0.33440312;//[50,65], ex=0.33440312, std=0.0379642==> ex_up=0.3723673126, ex_down=0.296438919764	
	double std_ex=0.0379642;
	double r;
	//output alpha	
	
	
	ofstream outfile22("data/tw_alpha_b_q_nv.txt",ios::out);	
	outfile22<<" "<<left<<setw(10)<<"# No."<<" "<<setw(10)<<"ex"<<" "<<setw(10)<<"std.ex"<<" "<<setw(10)<<"b"<<" "<<setw(10)<<"ave.a"<<" "<<setw(10)<<"std.a"<<" "<<setw(10)<<"ave.b^-nu"<<" "<<setw(10)<<"std.b^-nu";
	outfile22<<" "<<left<<setw(20)<<"ave_a_realization"<<" "<<setw(20)<<"std_a_realization"<<endl;
	for(int n=50;n<=65;n++)	
	{
		vector<double> a_set;
		vector<double> q_set;
		a_set.clear();
		q_set.clear();
		for(int ite=0;ite<100;ite++)
		{
			
			double meanC=0.0;
			double meanK=0.0;
			long int realN;			
			double beta, mu, mu_new;

			// input parameters, edgelist, coordinates from file
			ostringstream fin1;
			fin1 <<"data/TW_"<<n<<"_parameters.txt";
			string f_para= fin1.str();		
			Read_parameters(f_para, realN, beta, mu);

			ostringstream fin2;
			fin2 <<"data/TW_"<<n<<"_edgelist.txt";
			string f_edge= fin2.str();		
			edge_list=Read_Graph(f_edge);

			ostringstream fin3;
			fin3 <<"data/TW_"<<n<<"_coordinates.txt";
			string f_coord= fin3.str();			
			Read_kappa_theta(f_coord, kappa, theta, N);
			//***********************************************
			
			
			vector<double > Cc;
			Cc=cluster_coefficence(edge_list,N,meanC);
			double meanK0=meank_maxk(edge_list,N);
			mu_new=mu;

			int layer=0;
			realN=Real_N(edge_list,N);
			//IRG process begin*************

			//for(int layer=1;layer<2;layer++)
			layer=1;
			{			
				vector<long int> node_in_which_supernode;	
				vector<vector<long int> > supernode_contain_which_node;
				vector<double> kappa_l1;
				long int N;
				long int Nl1;			
				//change the file name of kappa from stable distribution
				ostringstream file1;
				file1 <<"data/TW_"<<n<<"_kappa_l_"<<layer;
				string File_name_kappa = file1.str();	
				kappa_l1=Read_kappa_from_stable_dist(File_name_kappa, node_in_which_supernode, supernode_contain_which_node, N, Nl1);
				//****return kappa_l1, node_in_which_supernode, supernode_contain_which_node, N, Nl1****
				
				//r will be change as Nl1; meank_0 in uplayer can be calculated using the fitting parameters; 
				r=(double)(Nl1)/(double)(N);
				
								
				//get the edge list in the up layer****
				edge_list=Growth(edge_list, N, Nl1, kappa, theta, kappa_l1, node_in_which_supernode, supernode_contain_which_node,r, beta, mu,mu_new);
				//return edgelist,kappa, theta, mu_new in uplayer  
				
				//update mu			
				mu=mu*r;

				//cal mean C, mean degree and real N	
				Cc=cluster_coefficence(edge_list, Nl1, meanC);
				meanK=meank_maxk(edge_list,Nl1);
				realN=Real_N(edge_list,Nl1);

				//output alpha
				double q=meanK/meanK0;
				q_set.push_back(q);
				
				double b=r;			
				double alpha=exp(ex*log(b)-log(q));
				a_set.push_back(alpha);
				
			}
		}		
		double ave_q;
		double std_q;		
		Cal_mean_std(q_set, ave_q, std_q);
		
		double ave_a_realization;
		double std_a_realization;	
		Cal_mean_std(a_set, ave_a_realization, std_a_realization);
		
		
		double b=r;			
		double ave_a=exp(ex*log(b)-log(ave_q));
		double sum_x=pow(ex,2)*pow(log(b),2)*pow((std_ex/ex),2)+pow((std_q/ave_q),2);
		double std_a=ave_a*sqrt(sum_x);
		
		
		outfile22<<" "<<left<<setw(10)<<n<<" "<<setw(10)<<ex<<" "<<setw(10)<<std_ex<<" "<<setw(10)<<r<<" "<<setw(10)<<ave_a<<" "<<setw(10)<<std_a<<" "<<setw(10)<<ave_q<<" "<<setw(10)<<std_q;				
		outfile22<<" "<<left<<setw(20)<<ave_a_realization<<" "<<setw(20)<<std_a_realization<<endl;
		cout<<" "<<left<<setw(10)<<n<<" "<<setw(10)<<ex<<" "<<setw(10)<<std_ex<<" "<<setw(10)<<r<<" "<<setw(10)<<ave_a<<" "<<setw(10)<<std_a<<" "<<setw(10)<<ave_q<<" "<<setw(10)<<std_q<<endl;				
		
	}
	outfile22.close();
	cout<<"finish calculating!"<<endl;
	double time_end=time(NULL);
    cout<<"RUN TIME: "<<time_end-time_begin<<endl;	
	return 0;
}
