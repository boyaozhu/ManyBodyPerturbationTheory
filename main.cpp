//
//  main.cpp
//  mbpt_pairing
//
//  Created by Boyao Zhu on 6/28/18.
//  Copyright © 2018 Boyao Zhu. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "armadillo"

using namespace std;
using namespace arma;

double assym(int, int, int, int, int [8][2], double&);
double eps(int, int, int, int, int [8][2]);
double e2nd(double (*assym) (int, int, int, int, int [8][2], double&), double (*eps) (int, int, int, int, int [8][2]), int [4], int [4], int [8][2], double&);

ofstream ofile;
int main(int argc, char* argv[])
{
    char *outfilename;
        // Read in output file, abort if there are too few command-line arguments
    if (argc<=1)
    {
        cout<<"Bad Usage: "<<argv[0]<<" Read also output file, number of integration points and the final x values on same line, four variables in total"<<endl;
        exit(1);
    }
    else 
    {
        outfilename=argv[1];
    }
    
    ofile.open(outfilename);
    
    int below_fermi[4] = {0,1,2,3};
    int above_fermi[4] = {4,5,6,7};
    int mystates[8][2] = {{1,1},{1,-1},{2,1},{2,-1},{3,1},{3,-1},{4,1},{4,-1}};
    double grid[21];
    grid[0]=-1;                  //initialize grid space from -1 to 1
    for (int i=1; i<21; i++)
    {
        grid[i]=grid[i-1]+0.1;
    }
    
    double s2;   //second order
    double g;    //interaction
    mat Hamiltonian(6,6);
    vec eigval;
    double val;
    ofile<<setw(15)<<"Grid Space"<<setw(15)<<"MBPT2"<<setw(15)<<"Exact"<<endl;
    for (int i=0; i<21; i++)
    {
        g=grid[i];
        s2=e2nd(&assym, &eps, below_fermi, above_fermi, mystates, g);
        
        ///////********************************************/////////////
        
        Hamiltonian<<2-g<<-g/2<<-g/2<<-g/2<<-g/2<<0<<endr
                   <<-g/2<<4-g<<-g/2<<-g/2<<0<<-g/2<<endr
                   <<-g/2<<-g/2<<6-g<<0<<-g/2<<-g/2<<endr
                   <<-g/2<<-g/2<<0<<6-g<<-g/2<<-g/2<<endr
                   <<-g/2<<0<<-g/2<<-g/2<<8-g<<-g/2<<endr
                   <<0<<-g/2<<-g/2<<-g/2<<-g/2<<10-g<<endr;
        eig_sym(eigval, Hamiltonian);
        
        val=min(eigval)-(2-g);
        //cout<<val<<endl;
        
        
    ofile<<setw(15)<<setprecision(4)<<grid[i]<<setw(15)<<setprecision(4)<<s2<<setw(15)<<setprecision(4)<<val<<endl;
        
    }
   
    ofile.close();
    
    return 0;
}







double assym(int p, int q, int r, int s, int mystates[8][2], double& pt)
{ //p,q,r,s start from 0,1,2,3,4,5,6.7
    
    int p1=mystates[p][0];
    int p2=mystates[q][0];
    int p3=mystates[r][0];
    int p4=mystates[s][0];
    int s1=mystates[p][1];
    int s2=mystates[q][1];
    int s3=mystates[r][1];
    int s4=mystates[s][1];
    if ((p1!=p2) || (p3!=p4))
        return 0;
    else if ((s1==s2) || (s3==s4))
        return 0;
    else if ((s1==s3) && (s2==s4))
        return -pt/2;
    //if ((s1==s4) && (s2==s3))
    else
        return pt/2;
}


double eps(int p, int q, int r, int s, int mystates[8][2])
{
    double E=0;
    int p1=mystates[p][0];
    int p2=mystates[q][0];
    int p3=mystates[r][0];
    int p4=mystates[s][0];
    int ptt[]={p1, p2, p3, p4};
    //cout<<ptt[0]<<endl<<ptt[1]<<endl<<ptt[2]<<endl<<ptt[3]<<endl;
    for(int i=0; i<4; i++)
    {
        if (ptt[i]<3)
            E+=((ptt[i])-1);
        else
            E-=((ptt[i])-1);
    }
    return E;
}

double e2nd(double (*assym) (int p, int q, int r, int s, int mystates[8][2], double& pt), double (*eps) (int p, int q, int r, int s, int mystates[8][2]), int below_fermi[4], int above_fermi[4], int mystates[8][2], double& g)
{
    double s1=0;
    for (int a=above_fermi[0]; a<=above_fermi[3]; a++)
    {
        for (int b=above_fermi[0]; b<=above_fermi[3]; b++)
        {
            for (int i=below_fermi[0]; i<=below_fermi[3]; i++)
            {
                for (int j=below_fermi[0]; j<=below_fermi[3]; j++)
                {
                    s1+=0.25*assym(a,b,i,j,mystates,g)*assym(i,j,a,b,mystates,g)/eps(a,b,i,j,mystates);
                    //cout<<s1<<endl;
                }
            }
        }
    }
    return s1;
}


