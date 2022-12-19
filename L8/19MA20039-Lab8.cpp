// 19MA20039
// Rahul Saini
// Cutting Plane method only

#include <bits/stdc++.h>
using namespace std;

#define lli long long
#define bM 100000000

pair<int, double> decompose(double x)
{
    if (x == 0.0)
    {
        return {0.0, 0.0};
    }
    double sign = x / abs(x);
    x = abs(x);
    double FracPart = x - (double)((int)(x));
    double IntPart = (int)(x);
    if (abs(FracPart - 1) < 1e-8)
        IntPart += 1, FracPart = 0;
    return {sign * IntPart, sign * FracPart};
}

int cutting_plane(double **A, double *B, double *C, int n, int m, int msum, int *sign, int maxmin) 
{
    double sum;
    int i, j, k = n + msum;
    int num_equations = m;
    int num_var = n + msum;
    int doublecheck = 0;

    // Now we need to generate the cutting_plane table first and foremost
    double **sim = (double **)malloc(100 * sizeof(double *));
    for (i = 0; i < 100; i++)
    {
        sim[i] = (double *)malloc(100 * sizeof(double));
    }

    k = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            sim[i][j] = A[i][j];
        }
        for (j = n; j < (n + msum); j++) 
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

    double *cb = (double *)malloc(((100) * sizeof(double))); 
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
    double maxincmz = 0.0;                                          // this helps to find the minimum Z_i - C_i in each iteration
    double mininsol = 999999.0;                                     // this helps to find the minimum ratio in each iteration
    int *basv = (int *)malloc(100 * sizeof(int));                   // stores the basic variables values
    double *Z = (double *)malloc(((100) * sizeof(double)));         //stores the values of Z = \sum (CB_i)*(a_i_j)
    double *CminusZ = (double *)malloc(((100) * sizeof(double)));   // calculates the Z_i - C_i
    double *sol = (double *)malloc(100 * sizeof(double));           // stores the solutions of each iteration
    double *ratio = (double *)malloc(100 * sizeof(double));         // calculates the ratios of each iteration
    double *keyrowval = (double *)malloc(((100) * sizeof(double))); // stores the key row in a separate array
    double *keycolval = (double *)malloc(100 * sizeof(double));     // stores the key column in a separate array
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

        //now we print the cutting_plane table for each iteration

        printf("\n\nIteration no: %d\n", iter);
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n\t CB_i \t C_j ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("\t %lf", cb[i]);
        }
        printf("\n \t \t BV. ");
        for (i = 0; i < (msum + n); i++)
        {
            printf("\t     x_%d", i + 1);
        }
        printf("\t Solution\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        for (i = 0; i < m; i++)
        {
            printf("\t %0.2lf    x_%d ", cb[basv[i]], basv[i] + 1);
            for (j = 0; j < (n + msum); j++)
            {
                printf("\t %lf ", sim[i][j]);
            }
            printf("\t %lf  \n", sol[i]);
        }
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i < (n + msum); i++)
        {
            printf("\t %lf", Z[i]);
        }
        printf("\t %lf", zsol);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (n + msum); i++)
        {
            if (maxmin == 0)
                printf("\t %lf ", -CminusZ[i]);
            else
            {
                printf("\t %lf", CminusZ[i]);
            }
        }
        printf("\n");
        for (i = 0; i < (22 * (m + n) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        printf("\nMinimum ratio is : %lf coming at pivot row : %d\n", mininsol, keyrow + 1);
        if (maxmin == 0)
            printf("Minimum Z-i - C_i is : %lf coming at pivot column: %d\n", -maxincmz, keycol + 1);
        else
        {
            printf("Maximum Z-i - C_i is : %lf coming at pivot column: %d\n", maxincmz, keycol + 1);
        }

        printf("Value of Z is :%lf\n", zsol);
        iter++;
    }

    for (i = 0; i < m; i++)
    {
        if (cb[basv[i]] < -1000000.0)
        {
            printf("The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (infeasible).\n");
            return 0;
        }
    }
    int *integer_sol = (int *)malloc(100 * sizeof(int));
    double *fractional_part = (double *)malloc(100 * sizeof(double));
    double maxfrac = -1;
    double *varibls = (double *)malloc((100) * sizeof(double));
    int maxfrac_ind;
    double frac_indicator;
    double intPart, fractPart;
    double epsilon = 1e-9;
    int flag = 1;
    printf("\n The final optimal values are : ");
    for (i = 0; i < m; i++)
    {

        printf(" x_ %d = %lf ", basv[i] + 1, sol[i]);
    }
    printf(" And rest all are 0\n And the optimal value of Z is : %lf\n", zsol);
    iter = 0;

    while (1)
    {
        if (iter++ > 6)
        {
            break;
        }
        printf("\n\n");
        maxfrac = -1.0;
        for (i = 0; i < num_equations; i++)
        {
            integer_sol[i] = decompose(sol[i]).first;
            fractional_part[i] = decompose(sol[i]).second;
            if (maxfrac < fractional_part[i])
            {
                maxfrac = fractional_part[i];
                maxfrac_ind = i;
            }
        }

        printf("\n*************************************************************************************\n");

        for (i = 0; i < num_equations; i++)
        {
            printf("Fractional part is %lf, and the integer part is %d of x_%d\n", fractional_part[i], integer_sol[i], basv[i] + 1);
        }
        for (i = 0; i < num_equations; i++)
        {
            sim[i][num_var] = 0.0;
        }
        printf("\n");

        printf("The maximum f_i is coming at the basic variable no:%d and it is %lf\nSo, we need to add the following Gomorian Constraint :\n\n", maxfrac_ind + 1, maxfrac);
        printf("%lf = ", -maxfrac);

        // now the number of variables and the number of equations need to be incremented
        num_var++;
        num_equations++;

        k = 0;

        for (i = 0; i < num_var; i++)
        {
            if (i < (num_var - 1))
            {
                fractPart = decompose(sim[maxfrac_ind][i]).second;
                if (abs(fractPart) > 0.9999)
                {
                    fractPart = 0.0;
                    varibls[i] = fractPart;
                    printf("(%lf) * x_%d + ", varibls[i], i + 1);
                    continue;
                }
                if (fractPart >= 0.0)
                {
                    varibls[i] = -fractPart;
                    printf("(%lf) * x_%d + ", varibls[i], i + 1);
                }
                else
                {
                    varibls[i] = -1.0 - fractPart;
                    printf("(%lf) * x_%d + ", varibls[i], i + 1);
                }
            }
            else
            {
                varibls[i] = 1.0;
            }
            sim[num_equations - 1][i] = varibls[i];
        }

        k = 0;

        printf(" x_%d\n\n", (num_var)); 

        // now we are going to use the dual simplex method
        printf("The number of variables now is : %d\t And the number of equations is %d\nSo, printing the table after adding the extra variable, we get :\n", num_var, num_equations);
        printf("\n");
        for (i = 0; i < (20 * (num_var) + 20); i++)
        {
            printf("-");
        }
        sol[num_equations - 1] = -maxfrac;
        cb[num_var - 1] = 0.0;
        basv[num_equations - 1] = (num_var - 1);

        printf("\n\t CB_i \t C_j ");
        for (i = 0; i < (num_var); i++)
        {
            printf("\t %lf", cb[i]);
        }
        printf("\n \t \t BV. ");
        for (i = 0; i < (num_var); i++)
        {
            printf("\t     x_%d", i + 1);
        }
        printf("\t Solution\n");
        for (i = 0; i < (20 * (num_var) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        for (i = 0; i < num_equations; i++)
        {

            printf("\t %0.2lf    x_%d ", cb[basv[i]], basv[i] + 1);
            for (j = 0; j < (num_var); j++)
            {
                printf("\t %lf ", sim[i][j]);
            }
            printf("\t %lf  \n", sol[i]);
        }
        for (i = 0; i < (20 * (num_var) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        for (i = 0; i < num_var; i++)
        {
            sum = 0;
            for (j = 0; j < m; j++) 
            {
                sum += cb[basv[j]] * sim[j][i];
            }
            Z[i] = sum;
            CminusZ[i] = cb[i] - Z[i];
        }
        printf("\n\t Z_j \t ");
        for (i = 0; i < (num_var); i++)
        {
            printf("\t %lf", Z[i]);
        }
        sum = 0;
        for (i = 0; i < num_equations; i++)
        {
            sum += cb[basv[i]] * sol[i]; 
        }
        zsol = sum;
        printf("\t %lf", zsol);
        printf("\n \t C_j - Z-j ");
        for (i = 0; i < (num_var); i++)
        {
            printf("\t %lf", CminusZ[i]);
        }
        printf("\n");
        for (i = 0; i < (20 * (num_var) + 20); i++)
        {
            printf("-");
        }
        printf("\n");
        int baniter = 0;
        
        // finding the leaving variable with the minimum value of sol[i]
        while (1)
        {
            if (baniter++ > 5)
                break;
            mininsol = 0;
            for (i = 0; i < num_equations; i++)
            {
                if (sol[i] < mininsol)
                {
                    mininsol = sol[i];
                    keyrow = i;
                }
                if (sol[i] < 0)
                {
                    check = 1;
                }
            }
            printf("The most negative value of Solution is coming at row %d, which is %lf. So the leaving variable is : x_%d \n", keyrow + 1, sol[keyrow], basv[keyrow] + 1);
            printf("\nWe now have to find the entering variable, so we compute the following table :\n");
            
            // finding the entering variable by taking ratio of leaving variable row and C[j] - Z[j]
            for (i = 0; i < (20 * (num_var) + 5); i++)
            {
                printf("-");
            }
            printf("\n");
            printf("Variables ");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t   x_%d       ", i + 1);
            }
            printf("\n");
            for (i = 0; i < (20 * (num_var) + 5); i++)
            {
                printf("-");
            }
            printf("\n -(C_j - Z-j)");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t %lf", -1 * CminusZ[i]);
            }
            printf("\n x_%d \t", keyrow + 1);
            for (i = 0; i < (num_var); i++)
            {
                printf("\t %lf", sim[keyrow][i]);
            }
            double minratio = -100000;
            double ratio2;
            printf("\n Ratio \t");
            for (i = 0; i < (num_var); i++)
            {
                if (sim[keyrow][i] < -0.0000001)
                {
                    ratio2 = (-1 * CminusZ[i]) / sim[keyrow][i];
                    printf("\t %lf", ratio2);
                    if (ratio2 > minratio)
                    {
                        minratio = ratio2;
                        keycol = i;
                    }
                }
                else
                {
                    printf("\t  --      ");
                }
            }
            printf("\n");
            for (i = 0; i < (20 * (num_var) + 5); i++)
            {
                printf("-");
            }
            printf("\nHere, the maximum value of the Ratio is %lf. So the entering variable is x_%d \n", minratio, keycol + 1);
            for (i = 0; i < (num_var); i++)
            {
                keyrowval[i] = sim[keyrow][i];
            }
            solkey = sol[keyrow];
            for (i = 0; i < num_equations; i++)
            {
                keycolval[i] = sim[i][keycol]; // now it evaluates the key colvalues
            }
            basv[keyrow] = keycol; //entering variable is put in the basic variable column
            printf("Key row is :%d , and key column is :%d\nwith the key element being %lf\n", keyrow, keycol, keyrowval[keycol]);
            for (i = 0; i < num_equations; i++)
            {
                if (i == keyrow)
                {
                    sol[i] = sol[i] / keyrowval[keycol];
                }
                else
                {
                    sol[i] = sol[i] - (keycolval[i] * solkey) / keyrowval[keycol]; // for other rows
                }

                for (j = 0; j < (num_var); j++)
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
            printf("\nNow updating the table we get :\n");
            printf("\n");
            for (i = 0; i < (22 * (num_var) + 20); i++)
            {
                printf("-");
            }
            printf("\n\t CB_i \t C_j ");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t %lf", cb[i]);
            }
            printf("\n \t \t BV. ");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t     x_%d", i + 1);
            }
            printf("\t Solution\n");
            for (i = 0; i < (22 * (num_var) + 20); i++)
            {
                printf("-");
            }
            printf("\n");
            for (i = 0; i < num_equations; i++)
            {

                printf("\t %0.2lf    x_%d ", cb[basv[i]], basv[i] + 1);
                for (j = 0; j < (num_var); j++)
                {
                    printf("\t %lf ", sim[i][j]);
                }
                printf("\t %lf  \n", sol[i]);
            }
            for (i = 0; i < (22 * (num_var) + 20); i++)
            {
                printf("-");
            }
            printf("\n");
            for (i = 0; i < num_var; i++)
            {
                sum = 0;
                for (j = 0; j < m; j++) 
                {
                    sum += cb[basv[j]] * sim[j][i];
                }
                Z[i] = sum;
                CminusZ[i] = cb[i] - Z[i]; 
            }
            printf("\n\t Z_j \t ");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t %lf", Z[i]);
            }
            sum = 0;
            for (i = 0; i < num_equations; i++)
            {
                sum += cb[basv[i]] * sol[i];
            }
            zsol = sum;
            printf("\t %lf", zsol);
            printf("\n \t C_j - Z-j ");
            for (i = 0; i < (num_var); i++)
            {
                printf("\t %lf", CminusZ[i]);
            }
            printf("\n");
            for (i = 0; i < (22 * (num_var) + 20); i++)
            {
                printf("-");
            }
            printf("\n");
            frac_indicator = 0.0001;
            doublecheck = 1;
            for (i = 0; i < num_equations; i++)
            {
                if (sol[i] < 0)
                {
                    doublecheck = 0;
                }
            }
            if (doublecheck)
            {
                break;
            }
            else
            {
                cout << "\nWe have to update the table again\n";
            }
        }

        flag = 1;
        for (i = 0; i < num_equations; i++)
        {
            if (basv[i] >= n)
                continue;
            fractPart = decompose(sol[i]).second;
            integer_sol[i] = decompose(sol[i]).first;
            fractional_part[i] = fractPart;
            printf("So, x_%d = %lf + %d\n", basv[i] + 1, fractional_part[i], decompose(sol[i]).first);
            if (fractional_part[i] > frac_indicator)
            {
                flag = 0;
            }
        }
        if (flag)
        {
            break;
        }
        else
        {
            cout << "Here, all the solutions don't have the fractional parts 0, so we need to do the steps again.\n";
        }
    }
    if (iter > 5)
    {
        cout << "The method wasn't able to find the Integer Solutions\n";
    }
    else
    {
        printf("\n The final optimal values are : ");
        for (i = 0; i < num_equations; i++)
        {

            printf(" x_%d = %0.0lf ", basv[i] + 1, sol[i]);
        }
        printf(" and rest all are 0\n The optimal value of Z => %lf\n", zsol);
    }
    
    return 0;
}
int main()
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
    int maxmin = 0, ck;
    scanf("%d", &maxmin);

    ck = cutting_plane(A, B, C, n, m, msum, sign, maxmin);

    return 0;
}