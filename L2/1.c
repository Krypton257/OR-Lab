// 19MA20039
// Rahul Saini

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void simplex(double **A, double *B, double *C, int n, int m) // Simplex method function
{
    double sum;
    int i, j, k = n + m;

    // Foramtion of Simplex Table

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + m); j++)   //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        sim[i][n + i] = 1.0; // the slack variables are added for each equation and their coefficient is set to be 1 for each of the m equations
    }
    double *cb = (double *)malloc(((n + m) * sizeof(double))); // this stores the coefficients of the basic variables which are present in the optimality condition
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i]; // all the normal variables have their same coefficients
    }
    for (j = n; j < (n + m); j++)
    {
        cb[j] = 0.0; // the slack variables have coefficients to be 0
    }
    int keycol, keyrow;
    double maxincmz = 0.0;                                            
    double mininsol = 999999.0;                                       // this helps to find the minimum ratio in each iteration
    
    int *basv = (int *)malloc(m * sizeof(int));                       // stores the Basic variables values
    double *Z = (double *)malloc(((n + m) * sizeof(double)));         // stores the values of Z = \sum (CB_i)*(a_i_j)
    double *CminusZ = (double *)malloc(((n + m) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));               // stores the solutions of each iteration
    double *ratio = (double *)malloc(m * sizeof(double));             // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + m) * sizeof(double))); // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));         // stores the key column in a separate array
    
    for (i = 0; i < m; i++)
    {
        basv[i] = i + n;    
        sol[i] = B[i];      
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;
    
    // Going through iterations

    while (check)
    {
        if (iter >= 3) // Breaks after second iteration
            break;
        if (iter != 0) // apart from the first iteration, we don't have to change the columns
        {
            basv[keyrow] = keycol;      
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; // for the pivot row
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; // for other rows
                }

                for (j = 0; j < (m + n); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; // for pivot row
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; // for other rows
                    }
                }
            }
        }

        check = 0;
        maxincmz = -90000.0;

        for (i = 0; i < (m + n); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) // Calculation of the Z values for each variable
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cb[i] - Z[i]; 
            if (CminusZ[i] > 0)        
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) // finds the key column
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; // calculates the sum for the solution
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; // finds the ratio in each iteration
            if(ratio[i] < 0)
            {
                continue;
            }
            if (ratio[i] < mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (m + n); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        // Printing Simplex table for each iteration

        printf("\n\nIteration no: %d\n", iter);
        printf("\n");
        printf("Basic_variables : ");
        for (i = 0; i < (m + n); i++)
        {
            printf(" %lf", cb[i]);
        }
        printf("\n");
        printf(" \n****************** Simplex Table ******************\n");
        for (i = 0; i < m; i++)
        {
            printf("%d %lf ", basv[i]+1, cb[basv[i]]);
            for (j = 0; j < (m + n); j++)
            {
                printf("%lf ", sim[i][j]);
            }
            printf("%lf %lf\n", sol[i], ratio[i]);
        }
        printf("***************************************************");
        printf("\n");

        printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
        printf("Value of Z is :%lf\n", zsol);
        iter++;
    }
    printf("\nThe optimal value of Z is : %lf\n", zsol);
}
int main() // menu driven program
{
    int n, m, i, j, k;
    printf("Enter number of variables : ");
    scanf("%d", &n);
    printf("Enter number of equations : ");
    scanf("%d", &m);
    double error = 0.001;
    double **A, *B;
   
    A = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m * sizeof(double));
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("Enter A[%d][%d] :", i, j);
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Enter the values of B\n");
    for (i = 0; i < m; i++)
    {
        printf("Enter B[%d] :", i);
        scanf("%lf", &B[i]);
    }
    printf("Now enter the coefficients of the function Z\n");
    double *C = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        printf("Enter coefficient of x_%d :", i + 1);
        scanf("%lf", &C[i]);
    }

    simplex(A, B, C, n, m);

    return 0;
}