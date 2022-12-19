//19MA20039
//Rahul Saini

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int combinationss[100][10];
int count;
int t = 1;

void printMat2d(int m , int n, double M[m][n])
{
    for(int i=0;i<m;i++)
    {
        for(int j =0;j<n;j++)
            printf("%f ",M[i][j]);
        printf("\n");
    }
    printf("\n");

}


long long int combinations(int n, int r)
{
    int i;
    long long int result = 1;
    if (r > n / 2)
        r = n - r;
    for (i = 1; i <= r; ++i)
    {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

int DominantTransformation(int m, int r, double M[m][m], double N[m], int V[], int R[])
{
    int i, j;
    int n = m;
    if (r == m)
    {
        double T[n][n + 1];
        double P[n];
        for (i = 0; i < n; i++)
        {
            P[i] = N[R[i]];
            for (j = 0; j < n; j++)
                T[i][j] = M[R[i]][j];
        }
        for (i = 0; i < n; i++)
        {
            N[i] = P[i];
            for (j = 0; j < n; j++)
                M[i][j] = T[i][j];
        }
        return 1;
    }
    for (i = 0; i < n; i++)
    {
        if (V[i])
            continue;
        double sum = 0;
        for (j = 0; j < n; j++)
            sum += fabs(M[i][j]);
        if (2 * fabs(M[i][r]) >= sum)
        {
            V[i] = 1;
            R[r] = i;
            if (DominantTransformation(m, r + 1, M, N, V, R))
                return 1;
            V[i] = 0;
        }
    }
    return 0;
}
int generateDom(int m, double M[m][m], double N[m])
{
    int i;
    int visited[m];
    for (i = 0; i < m; ++i)
        visited[i] = 0;
    int rows[m];
    return DominantTransformation(m, 0, M, N, visited, rows);
}
double *Gauss_Siedel(int m, int n, double a[m][n], double b[m], double er, int itr)
{
    int key = 1, i, j;
    double x0[n], sum;
    double *x = (double *)malloc(n * sizeof(double));
    double M[m][m];
    double N[m];
    for (i = 0; i < n; ++i)
    {
        x0[i] = 0;
        x[i] = 0;
    }
    for (i = 0; i < m; ++i)
    {
        N[i] = b[i];
        for (j = 0; j < m; ++j)
        {
            M[i][j] = a[i][combinationss[itr][j]];
        }
    }
    generateDom(m, M, N);
    int counter = 0;
    float allErrors;
    // printMat2d(m,m,M);
    while(key == 1)
    {
        key = 0;
        counter++;
        for (i = 0; i < m; ++i)
        {
            sum = N[i];
            for (j = 0; j < m; ++j)
                if (j != i)
                    sum -= M[i][j] * x0[combinationss[itr][j]];
            // printf("Sum = %f",sum);
            x[combinationss[itr][i]] = (sum / M[i][i]);
            allErrors = (x[i] - x0[i]) / x[i];

            if (allErrors < 0)
                allErrors *= (-1);
            

            if ((allErrors > er) && fabs((x[combinationss[itr][i]] - x0[combinationss[itr][i]]) /
                     x[combinationss[itr][i]]) > er)
            {
                key = 1;
                x0[combinationss[itr][i]] = x[combinationss[itr][i]];
            }

            if(counter > 20) // for failing case
            {
                t= 0;
                break;
            }
        }
    }
    return x;
}
void concatCombination(int arr[], int data[], int start, int end, int index, int r)
{
    int i, j;
    if (index == r)
    {
        for (j = 0; j < r; ++j)
        {
            combinationss[count][j] = data[j];
        }
        ++count;
        return;
    }
    for (i = start; i <= end && end - i + 1 >= r - index; ++i)
    {
        data[index] = arr[i];
        concatCombination(arr, data, i + 1, end, index + 1, r);
    }
}

void getCombination(int arr[], int n, int r)
{
    int data[r];
    concatCombination(arr, data, 0, n - 1, 0, r);
}


// void printMat1d(int m ,double M[m])
// {
//     for(int i=0;i<m;i++)
//     {
//         printf("%f ",M[i]);
//     }
//     printf("\n");

// }
//19MA20039

int main()
{

    int i, j, n, m;
    double er;

    printf("Enter stopping criteria :\n");
    scanf("%lf", &er);

    printf("Dimention of matrix (mxn):\n");
    scanf("%d%d", &m, &n);

    long long int solCount = combinations(n, m);

    printf("Elements of Matrix A : \n");

    double a[m][n];

    for (i = 0; i < m; ++i)
        for (j = 0; j < n; ++j)
            scanf("%lf", &a[i][j]);

    //printMat2d(m,n,a);

    printf("Elements of Matrix B : \n");
    double b[m];
    for (i = 0; i < m; ++i)
        scanf("%lf", &b[i]);

    //printMat1d(m,b);

    count = 0;
    int indices[n];
    for (i = 0; i < n; ++i)
        indices[i] = i;
    getCombination(indices, n, m);

    if (count != solCount)
        printf("error !! all combinations not found\n");
    double *x[solCount];

    printf("Solutions:\n\n");
    for (j = 0; j < solCount; ++j)
    {
        printf("Solution %d:\n", j + 1);
        x[j] = Gauss_Siedel(m, n, a, b, er, j);
        for (i = 0; i < n; ++i)
            printf("x_%d = %lf\t", i + 1, x[j][i]);
        // if(!t)
        //     break;
        printf("\n\n");
    }
    printf("Total Solutions = %d\n\n", count);
    return 0;
}
