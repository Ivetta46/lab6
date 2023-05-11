#include <iostream>

void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
}

double nextIteration()
{
    return 0;
}

double bottom()
{
    return 0;
}

double top()
{
    return 0;
}

double left()
{
    return 0;
}

double right()
{
    return 0;
}

double initial()
{
    return 0;
}

int main()
{
    int L, M, N;
    std::cin >> L >> M >> N;

    double step;
    std::cin >> step;

    double** grid = new double*[L];
    for (int l = 0; l < L; l++)
    {
        grid[l] = new double[M];
    }

    double*** u = new double**[N];
    for (int i = 0; i < N; i++)
    {   
        u[i] = new double*[L];
        
        for (int l = 0; l < L; l++)
        {
            u[i][l] = new double[M];
        }
    }

    for (int i = 0; i < N-1; i++)
    {
        u[i+1][0][0] = nextIteration();
    }

    for (int l = 0; l < L; l++)
    {
        delete[] grid[l];
    }
    
    delete[] grid;

    return 0;
}