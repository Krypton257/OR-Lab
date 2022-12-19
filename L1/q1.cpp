#include<stdio.h>

using namespace std;

double gauss_seidel(double m1[][],double m2[][],int m,int n, float error)

int main(){

	int m,n;
	cin>>m>>n;

	double m1[m][n];
	double m2[m];

	double error = 0.001;

	for(int i=0;i<m;i++){
		for(int j=0;j<n;j++){
			cin>>m1[i][j];
		}
	}

	for(int i=0;i<m;i++){
		cin>>m2[i];
	}




}

