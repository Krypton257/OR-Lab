#include <bits/stdc++.h>
using namespace std;

#define lli long long
#define bM 100000000

struct coordinates 
{
    int x, y;
};

void printTable(double *a, double *b, double **c, int m, int n, double **x)      // Prints the table
{
    cout << "\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n\t \t Depot\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n\t";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << i + 1;
    }
    cout << "\t Stock\n";
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n";
    for (int i = 0; i < m; i++)
    {
        cout << "A"<<i+1<<"\t| ";
        for (int j = 0; j < n; j++)
        {
            cout << "\t" << x[i][j]<<"("<<c[i][j]<<")";
        }

        cout << "\t" << a[i] << "\n\n";
    }
    for (int i = 0; i < n; i++)
    {
        cout << "--------------";
    }
    cout << "\n";
    cout << "Requirements :";
    for (int i = 0; i < n; i++)
    {
        cout << "\t" << b[i];
    }
    cout << "\n";
}

coordinates minInC(double **c, int m, int n)        
{
    double absolute = bM;
    coordinates min;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (c[i][j] < absolute)
            {
                min.x = i;
                min.y = j;
                absolute = c[i][j];
            }
        }
    }
    return min;
}

void MatrixMinimaMethod(double *a, double *b, double **c, int m, int n, double **x)    // Finds the BFS using Matrix Minima Method method
{
    double *aa = new double[m];
    double *bb = new double[n];
    double **cc = new double *[m];
    int row, col;
    for (int i = 0; i < m; i++)
    {
        cc[i] = new double[n];
    }
    double delivery = 0;
    for (int i = 0; i < m; i++)
    {
        aa[i] = a[i];
    }
    for (int i = 0; i < n; i++)
    {
        bb[i] = b[i];
        delivery += b[i]; 
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i][j] = 0; 
            cc[i][j] = c[i][j]; 
        }
    }
    coordinates min;
    int iter = 0;
    while (delivery > 0.0001) 
    {
        if (iter++ > 10)
            break;
        min = minInC(cc, m, n); // finds the minimum cell
        row = min.x;
        col = min.y;
        if (aa[min.x] > bb[min.y]) // and then the cell values are changed accordingly
        {
            for (int i = 0; i < m; i++)
            {
                cc[i][min.y] = bM;
            }
            x[row][col] = bb[col];
            aa[row] -= bb[col];
            bb[col] = 0.0;
            printTable(aa, bb, c, m, n, x);
        }
        else
        {
            for (int j = 0; j < n; j++)
            {
                cc[min.x][j] = bM;
            }
            x[row][col] = aa[row];
            bb[col] -= aa[row];
            aa[row] = 0.0;
            row++;
            printTable(aa, bb, c, m, n, x);
        }

        delivery = 0.0;
        for (int i = 0; i < n; i++)
        {
            delivery += bb[i];
        }
        cout << "Delivery = " << delivery << "\n"; 
    }

    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += c[i][j] * x[i][j];
        }
    }
    cout << "\nThe BFS from the Matrix Minima Method is : " << sum << "\n";
}

void northWestCorner(double *a, double *b, double **c, int m, int n, double **x)        // Finds BFS using North West Corner Method
{
    double *aa = new double[m];
    double *bb = new double[n];
    for (int i = 0; i < m; i++)
    {
        aa[i] = a[i];
    }
    for (int i = 0; i < n; i++)
    {
        bb[i] = b[i];
    }

    int col = 0, row = 0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            x[i][j] = 0;
        }
    }
    int iter = 0;
    while ((col != (n)) && (row != (m))) 
    {
        if (iter++ > (m + n))
        {
            break;
        }
        if (aa[row] > bb[col]) 
        {
            x[row][col] = bb[col];
            aa[row] -= bb[col];
            bb[col] = 0.0;
            col++;
            printTable(aa, bb, c, m, n, x);
        }
        else // other wise this is done
        {
            x[row][col] = aa[row];
            bb[col] -= aa[row];
            aa[row] = 0.0;
            row++;
            printTable(aa, bb, c, m, n, x);
        }
    }
    double sum = 0.0;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            sum += c[i][j] * x[i][j];
        }
    }
    cout << "\nThe BFS from the North West Corner method is : " << sum << "\n";
    cout<<"\n";
}

int main()
{
    int n, m;
    cout << "Enter number of Units : ";
    cin >> m;
    cout << "Enter number of Depots : ";
    cin >> n;

    double *a = new double[m];          // Stores the Stock of the Units
    double *b = new double[n];          // Stores the Requirements of the Depots
    double **x = new double *[100];     

    for (int i = 0; i < 100; i++) 
    { 
        x[i] = new double[100];
    }

    double asum = 0.0, bsum = 0.0;
    cout << "-- INPUT Stock --\n";
    for (int i = 0; i < m; i++)
    {
        cout << "Stock at Unit " << i + 1 << " :";
        cin >> a[i];
        asum += a[i];
    }
    
    cout << "-- INPUT Requirement--\n";
    for (int i = 0; i < n; i++)
    {
        cout << "Requirement at Depot " << i + 1 << " :";
        cin >> b[i];
        bsum += b[i];
    }
    if (abs(bsum - asum) < 0.001)
    {
        cout << "Transportation problem is balanced\n";
    }
    else
    {
        cout << "The problem is not balanced\n";
        return 0;
    }

    double **c = new double *[m];           // Stores the Quantity of Unit from Unit i to Depot j in Matrix Form
    for (int i = 0; i < m; i++)
    {
        c[i] = new double[n];
    }
    cout << "Enter Quantity of Unit i to be shipped to the Depot j from the Unit i in matrix form\n";
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> c[i][j];
            x[i][j] = 0;
        }
    }

    printTable(a, b, c, m, n, x);
    MatrixMinimaMethod(a, b, c, m, n, x);
    printTable(a, b, c, m, n, x);
    northWestCorner(a, b, c, m, n, x);

}