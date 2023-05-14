#include <iostream>
#include "conditions.hpp"
#include <algorithm> 

double compare(double** left, double** right, int L, int M)
{   
    double maxDif = -1;
    double max = -1;
    for (int l = 1; l < L - 1; l++)
    {
        for (int m = 1; m < M - 1; m++)
        {
            maxDif = std::max(maxDif, std::abs(left[m][l] - right[m][l]));
            max = std::max(max, std::abs(left[m][l]));
        }
    }
    std::cout << "maxdiff: " << maxDif << std::endl;
    std::cout << "max: " << max << std::endl;
    return maxDif / max;
}

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

void nextIteration(
    double* a, 
    double* b, 
    double* c, 
    int L, 
    int M, 
    double* kPlus, 
    double* kMinus, 
    double** initialLayer,
    double** currentLayer,
    double** intermediateLayer,
    double** nextLayer)
{
    while(true)
    {
        for (int m = 1; m < M - 1; m++){
            // kminus = (k1 + k0) / 2  // 1..L-1
            // kplus  = (k1 + k2) / 2  // 1..L-1
            for (int l = 1; l < L - 1; l++)
            {       
                kMinus[l] = 0.5 * (K(currentLayer[m][l - 1]) + K(currentLayer[m][l]));
                kPlus[l] = 0.5 * (K(currentLayer[m][l]) + K(currentLayer[m][l + 1]));
            }

            // Заполним матрицу СЛАУ для первого слоя
            for (int l = 1; l < L - 1; l++)
            {
                c[l - 1] = 1 + tau / (hx * hx) * (kMinus[l] + kPlus[l]); // Главная диагональ
            }
            
            b[L - 3] = 0;
            for (int l = 1; l < L - 2; l++)
            {
                b[l - 1] = - tau / (hx * hx) * kMinus[l]; // Диагональ над главной
            }

            a[0] = 0;
            for (int l = 2; l < L - 1; l++)
            {
                a[l - 1] = - tau / (hx * hx) * kPlus[l]; // Диагональ под главной
            }

            double* F = new double[L-2]; 
            double* result = new double[L-2]; 

            for (int l = 1; l < L - 1; l++)
            {
                F[l - 1] =  initialLayer[m][l];
            }

            solveMatrix(L - 2, a, c, b, F, result);        
            
            for (int l = 1; l < L - 1; l++)
            {
                intermediateLayer[m][l] = result[l - 1];
            }

            delete[] F;  
            delete[] result;
        }

        for (int l = 1; l < L - 1; l++){
            // kminus = (k1 + k0) / 2 = (u0 + u1) / 2 // 1..L-1
            // kplus  = (k1 + k2) / 2 = (u1 + u2) / 2 // 1..L-1
            for (int m = 1; m < M - 1; m++)
            {       
                kMinus[m] = 0.5 * (K(currentLayer[m - 1][l]) + K(currentLayer[m][l]));
                kPlus[m] = 0.5 * (K(currentLayer[m][l]) + K(currentLayer[m + 1][l]));
            }

            // Заполним матрицу СЛАУ для второго слоя
            for (int m = 1; m < M - 1; m++)
            {
                c[m - 1] = 1 + tau / (hy * hy) * (kMinus[m] + kPlus[m]); // Главная диагональ
            }
            
            b[M - 3] = 0;
            for (int m = 1; m < M - 2; m++)
            {
                b[m - 1] = - tau / (hy * hy) * kMinus[m]; // Диагональ над главной
            }

            a[0] = 0;
            for (int m = 2; m < M - 1; m++)
            {
                a[m - 1] = - tau / (hy * hy) * kPlus[m]; // Диагональ под главной
            }

            double* F = new double[M-2]; 
            double* result = new double[M-2]; 

            for (int m = 1; m < M - 1; m++)
            {
                F[m - 1] = intermediateLayer[m][l];
            }

            solveMatrix(M - 2, a, c, b, F, result);        
            
            for (int m = 1; m < M - 1; m++)
            {
                nextLayer[m][l] = result[m - 1];
            }

            delete[] F;
            delete[] result;
        }

        std::cout << compare(nextLayer, currentLayer, L, M) << std::endl;

        if (compare(nextLayer, currentLayer, L, M) < epsilon)
        {
            break;
        }

        for (int l = 1; l < L - 1; l++)
        {
            for (int m = 1; m < M - 1; m++)
            {
                currentLayer[m][l] = nextLayer[m][l];
            }
        }  
    }
}

int main()
{

    /////////////////////////////////////////////////////////////
    ///-----------------Allocate resources--------------------///
    /////////////////////////////////////////////////////////////

    int L, M, N;
    std::cin >> L >> M >> N;

    tau = 1.0 / (N - 1);
    hx = 1.0 / (L - 1);
    hy = 1.0 / (M - 1);

    epsilon = 0.00001;

    double*** u = new double**[N];
    for (int i = 0; i < N; i++)
    {   
        u[i] = new double*[M];
        
        for (int m = 0; m < M; m++)
        {
            u[i][m] = new double[L];
        }
    }

    /////////////////////////////////////////////////////////////
    ///-----------------Set initial values--------------------///
    /////////////////////////////////////////////////////////////

    for (int l = 0; l < L; l++)
    {
        for (int m = 0; m < M; m++)
        {
            u[0][m][l] = initial(l * hx, m * hy);
        }
    }

    for (int n = 1; n < N; n++)
    {
        for (int l = 0; l < L; l++)
        {
            u[n][0][l] = bottom(l * hx, (n) * tau);
            u[n][M - 1][l] = top(l * hx, (n) * tau);
        }

        for (int m = 0; m < M; m++)
        {
            u[n][m][0] = left(m * hy, (n) * tau);
            u[n][m][L - 1] = right(m * hy, (n) * tau);
        }
    }
    
    /////////////////////////////////////////////////////////////
    ///-----------------Solve diff equation-------------------///
    /////////////////////////////////////////////////////////////

    double* a = new double[std::max(L, M)];
    double* b = new double[std::max(L, M)];
    double* c = new double[std::max(L, M)];
    double* kPlus = new double[std::max(L, M)];
    double* kMinus = new double[std::max(L, M)];

    double** currentLayer = new double*[M];
    double** intermediateLayer = new double*[M];
    double** nextLayer = new double*[M];

    for (int m = 0; m < M; m++)
    {
        currentLayer[m] = new double[L];
        intermediateLayer[m] = new double[L];
        nextLayer[m] = new double[L];
    }

    for (int n = 1; n < N; n++)
    {
        std::cout << n << std::endl;
        for (int m = 0; m < M; m++)
        {
            for (int l = 0; l < L; l++)
            {
                currentLayer[m][l] = u[n - 1][m][l];
            }
        }
        
        nextIteration(a, b, c, L, M, kPlus, kMinus, u[n - 1], currentLayer, intermediateLayer, nextLayer);
        
        for (int m = 1; m < M - 1; m++)
        {
            for (int l = 1; l < L - 1; l++)
            {
                u[n][m][l] = nextLayer[m][l];
            }
        } 
    }


    for (int m = 0; m < M; m++)
    {
        for (int l = 0; l < L; l++)
        {
            std::cout << u[N - 1][m][l] << " ";
        }
        std::cout << std::endl;
    } 

    /////////////////////////////////////////////////////////////
    ///-----------------Deallocate resources------------------///
    /////////////////////////////////////////////////////////////

    return 0;
}