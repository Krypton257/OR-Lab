#include<iostream>
#include <iomanip> 
using namespace std;


void print_table(float **a,int x_b[] ,float c_b[],float cost[],float b[], float delta[], int m, int n, int iteration){
  cout<<endl;
  cout<<"*************** SIMPLEX TABLE FOR "<<iteration<<"th ITERATION ***************"<<endl;
  cout<<endl;
  cout<<" \t \t";
  for(int i=1;i<=m+n;i++)
    cout<<cost[i]<<"\t";
  cout<<endl;
  cout<<"BV\tCB\t";
  for(int i=1;i<=m+n;i++)
    cout<<"X"<<i<<"\t";
  cout<<"B"<<endl;
  for(int i=0;i<m;i++){
    cout<<x_b[i+1]<<"\t"<<c_b[i+1]<<"\t";
    for(int j=0;j<n+m;j++){
      cout<<setprecision(4)<<a[i][j]<<"\t";
    }
    cout<<b[i]<<endl;
  }
  cout<<" \t \t";
  for(int i=1;i<=m+n;i++)
    cout<<delta[i]<<"\t";
  cout<<endl;

}

void change_table(float **a, float b[], int m, int n, int q, int p){
  float temp=a[p][q];
  b[p]=b[p]/a[p][q];
  for(int i=0;i<n+m;i++){
    a[p][i]=a[p][i]/temp;
  }
  for(int i=0;i<m;i++){
    if(i!=p){
      float f= a[i][q];
      for(int j=0;j<n+m;j++){
        a[i][j]=a[i][j]-f*a[p][j];
      }
      b[i]=b[i]-f*b[p];
    }
  }
}
void simplex_table(float **a,float b[],int n,int m,float obj_func[]){
  int x_b[m+1];
  float c_b[m+1];
  float cost[n+m+1];
  float delta[n+m+1];
  for(int i=1;i<=m;i++)
    x_b[i]=n+i;
  for(int i=1;i<=m;i++)
    c_b[i]=0;
  for(int i=1;i<=n;i++){
    cost[i]=obj_func[i-1];
  }
  for(int i=n+1;i<=n+m;i++)
    cost[i]=0;
  float min=0;
  float z[n+m+1];
  int min_index;
  float b_min=10000000;
  int b_min_index;
  int interation_count=0;
  while(1){
    
    for(int i=0;i<n+m;i++){
      z[i+1]=0;
      for(int j=0;j<m;j++){
        z[i+1]+=c_b[j+1]*a[j][i];
      }
    }
    for(int i=1;i<=m+n;i++){
      delta[i]=z[i]-cost[i];
    }
    print_table(a,x_b,c_b,cost,b,delta,m,n,interation_count);
    min=0;
    for(int i=1;i<=m+n;i++){
      if(delta[i]<min){
        min=delta[i];
        min_index=i;
      }
    }
    if(min==0){
      break;
    }
    b_min=10000000;
    for(int i=0;i<m;i++){
      if(a[i][min_index-1]>0 && b[i]>=0){
        if(b_min> b[i]/a[i][min_index-1]){
          b_min=b[i]/a[i][min_index-1];
          b_min_index=i+1;
        }
      }
    }
    x_b[b_min_index]=min_index;
    c_b[b_min_index]=cost[min_index];
    change_table(a,b,m,n,min_index-1,b_min_index-1);
    interation_count++;
  }
}

int main(){
  int n,m;
  cout<<"Enter number of variables to maximize i.e n :";
  cin>>n;
  cout<<"Enter number of equations i.e value of m :" ;
  cin>>m;
  float **a = new float* [m];
  for(int i=0;i<m;i++)
    a[i]= new float[n+m];
  cout<<"Enter the coefficients on left: "<<endl;
  for(int i=0;i<m;i++)
    for(int j=0;j<n;j++)
      cin>>a[i][j];
  cout<<"Enter the values of bj's :"<<endl;
  float b[m];
  for(int i=0;i<m;i++)
    cin>>b[i];
  cout<<"Enter the coefficients of the expression you want to maximize: "<<endl;
  float obj_func[n];
  for(int i=0;i<n;i++)
    cin>>obj_func[i];
  for(int j=0;j<m;j++){
    for(int i=n;i<n+m;i++){
      if(i-j==n)
        a[j][i]=1;
      else
        a[j][i]=0;
    }
  }
  simplex_table(a,b,n,m,obj_func);
  return 0;

  
}