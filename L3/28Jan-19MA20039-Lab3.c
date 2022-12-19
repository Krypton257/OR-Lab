// 19MA20039
// Rahul Saini

//Example
//Minimize Z=4x1+8x2+3x3, Subject to x1+x2≥2, 2x1+x3≤5

// Input :
// Enter number of variables : 3
// Enter number of equations : 2
// Enter the values of A :
// 1 1 0
// 2 0 1
// Enter the values of B:
// 2 5
// Enter (0 for <=) OR (1 for >=) OR (2 for =) in each of the m= 2 equations
// 1 0
// Enter the coefficients of the function Z
// 4 8 3
// Enter (0 for Maximization) OR (1 for Minimization) : 1

// Code is bit lengthy because Same function is repeated multiple times for the menu system view


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define bM 100000000 

void BFS(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) 
{
    double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }

        printf("\n");
        printf("Basic Variables : ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("  %lf", cb[i]);
        }
        printf("\n");

}

void iterations(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) 
{
    double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }
    int keycol, keyrow;
    double maxincmz = 0.0;                                               
    double mininsol = 999999.0;                                             // this helps to find the minimum ratio in each iteration                                 
    int *basv = (int *)malloc(m * sizeof(int));                             // stores the Basic variables values                
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));            // stores the values of Z = \sum (CB_i)*(a_i_j)        
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));                     // stores the solutions of each iteration           
    double *ratio = (double *)malloc(m * sizeof(double));                   // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double)));    // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));               // stores the key column in a separate array
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; 
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; 
            k++;
        }
        else
        {
            basv[i] = k + n; 
            k++;
        }
        sol[i] = B[i]; 
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;

// Going through iterations

    while (check)
    {
        if (iter >= 10) 
        {
            printf("There are no feasable solutions\n");
            break;
        }
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; 
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; 
                }

                for (j = 0; j < (n + msum); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; 
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; 
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; 
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) 
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) 
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; 
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; 
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        iter++;
    }
    printf("\n");
    printf("Required iterations : %d\n", iter);
    printf("\n");
}

void nbvariables(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin ,int serial) 
{
    double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }

    int keycol, keyrow;
    double maxincmz = 0.0;                                               
    double mininsol = 999999.0;                                             // this helps to find the minimum ratio in each iteration                                 
    int *basv = (int *)malloc(m * sizeof(int));                             // stores the Basic variables values                
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));            // stores the values of Z = \sum (CB_i)*(a_i_j)        
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));                     // stores the solutions of each iteration           
    double *ratio = (double *)malloc(m * sizeof(double));                   // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double)));    // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));               // stores the key column in a separate array
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; 
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; 
            k++;
        }
        else
        {
            basv[i] = k + n; 
            k++;
        }
        sol[i] = B[i]; 
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;

//Going through iterations

    while (check)
    {
        if (iter >= 10) 
        {
            printf("There are no feasable solutions\n");
            break;
        }
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; 
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; 
                }

                for (j = 0; j < (n + msum); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; 
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; 
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; 
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) 
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) 
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; 
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; 
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        if(iter==serial)
        {
            printf("\n\nIteration no: %d\n", iter);
            printf("\n");

            printf("  Z_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                printf("%lf  ", Z[i]);
            }
            printf("\nZ-i - C_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                if (maxmin == 0)
                    printf("%lf  ", -CminusZ[i]);
                else
                {
                    printf("%lf  ", CminusZ[i]);
                }
            }
            printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
            if (maxmin == 0)
                printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -maxincmz, keycol + 1);
            else
            {
                printf("Maximum Z-i - C_i is : %lf coming at pivot column: %d\n", maxincmz, keycol + 1);
            }

            printf("Value of Z is :%lf\n", zsol);
        }
            iter++;
    }
    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            break ;
        }
    }
}

void bvariables(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin ,int serial) 
{
    double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }

    int keycol, keyrow;
    double maxincmz = 0.0;                                               
    double mininsol = 999999.0;                                             // this helps to find the minimum ratio in each iteration                                 
    int *basv = (int *)malloc(m * sizeof(int));                             // stores the Basic variables values                
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));            // stores the values of Z = \sum (CB_i)*(a_i_j)        
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));                     // stores the solutions of each iteration           
    double *ratio = (double *)malloc(m * sizeof(double));                   // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double)));    // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));               // stores the key column in a separate array
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; 
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; 
            k++;
        }
        else
        {
            basv[i] = k + n; 
            k++;
        }
        sol[i] = B[i]; 
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;

//Going through iterations

    while (check)
    {
        if (iter >= 10) 
        {
            printf("There are no feasable solutions\n");
            break;
        }
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; 
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; 
                }

                for (j = 0; j < (n + msum); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; 
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; 
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; 
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) 
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) 
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; 
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; 
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        if(iter==serial)
        {
            printf("\n\nIteration no: %d\n", iter);
            printf("\n");

            printf("  Z_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                printf("%lf  ", Z[i]);
            }
            printf("\nZ-i - C_i  ");
            for (i = 0; i < (n + msum); i++)
            {
                if (maxmin == 0)
                    printf("%lf  ", -CminusZ[i]);
                else
                {
                    printf("%lf  ", CminusZ[i]);
                }
            }
            printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
            if (maxmin == 0)
                printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -maxincmz, keycol + 1);
            else
            {
                printf("Maximum Z-i - C_i is : %lf coming at pivot column: %d\n", maxincmz, keycol + 1);
            }

            printf("Value of Z is :%lf\n", zsol);
        }
            iter++;
    }
    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            break ;
        }
    }
}

void table(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin, int serial) 
{
   double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }
    
    int keycol, keyrow;
    double maxincmz = 0.0;                                               
    double mininsol = 999999.0;                                             // this helps to find the minimum ratio in each iteration                                 
    int *basv = (int *)malloc(m * sizeof(int));                             // stores the Basic variables values                
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));            // stores the values of Z = \sum (CB_i)*(a_i_j)        
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));                     // stores the solutions of each iteration           
    double *ratio = (double *)malloc(m * sizeof(double));                   // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double)));    // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));               // stores the key column in a separate array
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; 
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; 
            k++;
        }
        else
        {
            basv[i] = k + n; 
            k++;
        }
        sol[i] = B[i]; 
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;

//Going through iterations

    while (check)
    {
        if (iter >= 10) 
        {
            printf("There are no feasable solutions\n");
            break;
        }
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; 
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; 
                }

                for (j = 0; j < (n + msum); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; 
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; 
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; 
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) 
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) 
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; 
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; 
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        if(iter==serial)
        {
            printf("\n\nIteration no: %d\n", iter);

            printf("\n");

            printf(" \n****************** Simplex Table ******************\n");
            
            for (i = 0; i < m; i++)
            {
                printf("%lf  %d  ", cb[basv[i]], basv[i] + 1);
                for (j = 0; j < (n + msum); j++)
                {
                    printf("%lf  ", sim[i][j]);
                }
                printf("%lf  %lf\n", sol[i], ratio[i]);
            }
            printf("***************************************************");
            printf("\n");
        }

        iter++;
    }
    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            break ;
        }
    }
}

void optimal_sol(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) 
{
    double sum;
    int i, j, k = n + msum;

    double **sim = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
    {
        sim[i] = (double *)malloc(k * sizeof(double));
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++)    //Coefficients for the Slack variables
        {
            sim[i][j] = 0.0;
        }
        
        if (sign[i] == 1)
        {
            sim[i][n + k] = -1.0; 
            k++;
            sim[i][n + k] = 1.0; 
            k++;
        }
        else
        {
            sim[i][n + k] = 1.0;       
            k++;
        }
    }
    double *cb = (double *)malloc(((n + msum) * sizeof(double)));       
    for (i = 0; i < n; i++)
    {
        cb[i] = C[i];   
    }
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            cb[n + k] = 0.0; 
            k++;
        }
        else if (sign[i] == 1)
        {
            cb[n + k] = 0.0; 
            k++;
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }

            k++;
        }
        else
        {
            if (maxmin == 0)
                cb[n + k] = -bM; 
            else
            {
                cb[n + k] = bM; 
            }
            k++;
        }
    }
    int keycol, keyrow;
    double maxincmz = 0.0;                                               
    double mininsol = 999999.0;                                             // this helps to find the minimum ratio in each iteration                                 
    int *basv = (int *)malloc(m * sizeof(int));                             // stores the Basic variables values                
    double *Z = (double *)malloc(((n + msum) * sizeof(double)));            // stores the values of Z = \sum (CB_i)*(a_i_j)        
    double *CminusZ = (double *)malloc(((n + msum) * sizeof(double)));   
    double *sol = (double *)malloc(m * sizeof(double));                     // stores the solutions of each iteration           
    double *ratio = (double *)malloc(m * sizeof(double));                   // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((n + msum) * sizeof(double)));    // stores the key row in a separate array
    double *keycolval = (double *)malloc(m * sizeof(double));               // stores the key column in a separate array
    k = 0;
    for (i = 0; i < m; i++)
    {
        if (sign[i] == 0)
        {
            basv[i] = k + n; 
            k++;
        }
        else if (sign[i] == 1)
        {
            k++;
            basv[i] = k + n; 
            k++;
        }
        else
        {
            basv[i] = k + n; 
            k++;
        }
        sol[i] = B[i]; 
    }
    int check = 1;
    int iter = 0;
    double solkey;
    double zsol;

//Going through iterations

    while (check)
    {
        if (iter >= 10) 
        {
            printf("There are no feasable solutions\n");
            break;
        }
        if (iter != 0) 
        {
            basv[keyrow] = keycol; 
            for (i = 0; i < m; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol]; 
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; 
                }

                for (j = 0; j < (n + msum); j++)
                {
                    if (i == keyrow)
                    {
                        sim[i][j] = sim[i][j] / keyrowval[keycol]; 
                    }
                    else
                    {
                        sim[i][j] = sim[i][j] - (keycolval[i] * keyrowval[j]) / keyrowval[keycol]; 
                    }
                }
            }
        }
        check = 0;
        maxincmz = -90000.0;
        for (i = 0; i < (n + msum); i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            if (maxmin == 0)
                CminusZ[i] = cb[i] - Z[i]; 
            else
            {
                CminusZ[i] = Z[i] - cb[i];
            }

            if (CminusZ[i] > 0) 
            {
                check = 1;
            }
            if (CminusZ[i] > maxincmz) 
            {
                keycol = i;
                maxincmz = CminusZ[i];
            }
        }
        sum = 0;
        for (i = 0; i < m; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        mininsol = 999999.0;
        for (i = 0; i < m; i++)
        {
            keycolval[i] = sim[i][keycol]; 
        }
        for (i = 0; i < m; i++)
        {
            ratio[i] = sol[i] / keycolval[i]; 
            if ((sol[i] == 0.000) && (keycolval[i] < 0))
            {
                continue;
            }
            if (ratio[i] < 0.0000000)
            {
                continue;
            }
            if (ratio[i] <= mininsol)
            {
                mininsol = ratio[i];
                keyrow = i;
            }
        }
        for (i = 0; i < (n + msum); i++)
        {
            keyrowval[i] = sim[keyrow][i];
        }
        solkey = sol[keyrow];

        iter++;
    }
    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            break ;
        }
    }

    printf("\nThe optimal value of Z is : %lf\n", zsol);
}

int main() 
{
    int n, m, i, j, k, itr1;

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

    int *sign = (int *)malloc(m * sizeof(int));

    B = (double *)malloc(m * sizeof(double));

    printf("Enter the values of A : \n");

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Enter the values of B: \n");
    for (i = 0; i < m; i++)
    {
        scanf("%lf", &B[i]);
    }

    int msum = 0;
    printf("Enter (0 for <=) OR (1 for >=) OR (2 for =) in each of the m= %d equations \n", m);
    for (i = 0; i < m; i++)
    {
        scanf("%d", &sign[i]);
        if (sign[i] == 1)
        {
            msum += 2;
        }
        else
        {
            msum++;
        }
    }

    printf("Enter the coefficients of the function Z\n");
    double *C = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        scanf("%lf", &C[i]);
    }

    printf("Enter (0 for Maximization) OR (1 for Minimization) : ");
    int maxmin = 0;
    scanf("%d", &maxmin);

    while(1){

        printf("\n");
        printf("Enter the serial number for your operation :\n");
        printf("1. List of all BFS\n");
        printf("2. Number of iterations to solve the problem\n");
        printf("3. List of all Non-basic variables along with net evaluations in ith (user input) iteration\n");
        printf("4. List of Basic variables along with min ratios in ith iteration\n");
        printf("5. Simplex table of ith (user input) iteration\n");
        printf("6. Optimal solution (if exists otherwise generate report for infeasibility, unboundedness, alternative optimum etc.)\n");
        printf("---- Enter any other number to EXIT ----\n");
        printf("\n");

        int num2;
        scanf("%d",&num2);

        if(num2==1){
            BFS(A, B, C, n, m, msum, sign, maxmin);
        }
        else if(num2==2){
            iterations(A, B, C, n, m, msum, sign, maxmin);
        }
        else if(num2==3){

            int serial;
            printf("Enter the value of i for ith iteration : ");
            scanf("%d",&serial);

            nbvariables(A, B, C, n, m, msum, sign, maxmin,serial);
        }
        else if(num2==4){

            int serial;
            printf("Enter the value of i for ith iteration : ");
            scanf("%d",&serial);

            bvariables(A, B, C, n, m, msum, sign, maxmin,serial);
        }
        else if(num2==5){

            int serial;
            printf("Enter the value of i for ith iteration : ");
            scanf("%d",&serial);

            table(A, B, C, n, m, msum, sign, maxmin,serial);
        }
        else if(num2==6){
            optimal_sol(A, B, C, n, m, msum, sign, maxmin);
        }
        else{
            return 0;
        }
    }

}