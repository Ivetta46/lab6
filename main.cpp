#include <iostream>
#include "conditions.hpp"
#include <algorithm> 

double compare(double** left, double** right, int L, int M)
{   
    double maxDif = -1;
    double max = -1;
    for (int l = 0; l < L; l++)
    {
        for (int m = 0; m < M; m++)
        {
            maxDif = std::max(maxDif, std::abs(left[m][l] - right[m][l]));
            max = std::max(max, std::abs(right[m][l]));
        }
    }

    return maxDif / max;
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
        for (int m = 0; m < M; m++){
            // kminus = (k1 + k0) / 2 = (u0 + u1) / 2 // 1..L-1
            // kplus  = (k1 + k2) / 2 = (u1 + u2) / 2 // 1..L-1
            for (int l = 0; l < L - 1; l++)
            {       
                kMinus[l] = 0.5 * (K(currentLayer[m][l]) + K(currentLayer[m][l + 1]));
                kPlus[l] = 0.5 * (K(currentLayer[m][l + 1]) + K(currentLayer[m][l + 2]));
            }

            // Заполним матрицу СЛАУ для первого слоя
            for (int l = 0; l < L; l++)
            {
                c[l] = 1 + tau / (hx * hx) * (kMinus[l] + kPlus[l]); // Главная диагональ
            }
            
            b[L - 1] = 0;
            for (int l = 0; l < L - 1; l++)
            {
                b[l] = - tau / (hx * hx) * kPlus[l]; // Диагональ над главной
            }

            a[0] = 0;
            for (int l = 1; l < L; l++)
            {
                a[l] = - tau / (hx * hx) * kMinus[l]; // Диагональ под главной
            }

            solveMatrix(L, a, c, b, initialLayer[m], intermediateLayer[m]);        
        }

        for (int l = 0; l < L; l++){
            // kminus = (k1 + k0) / 2 = (u0 + u1) / 2 // 1..L-1
            // kplus  = (k1 + k2) / 2 = (u1 + u2) / 2 // 1..L-1
            for (int m = 0; m < M - 1; m++)
            {       
                kMinus[m] = 0.5 * (K(currentLayer[m][l]) + K(currentLayer[m + 1][l]));
                kPlus[m] = 0.5 * (K(currentLayer[m + 1][l]) + K(currentLayer[m + 2][l]));
            }

            // Заполним матрицу СЛАУ для второго слоя
            for (int m = 0; m < M; m++)
            {
                c[m] = 1 + tau / (hy * hy) * (kMinus[m] + kPlus[m]); // Главная диагональ
            }
            
            b[M - 1] = 0;
            for (int m = 0; m < M - 1; m++)
            {
                b[m] = - tau / (hy * hy) * kPlus[m]; // Диагональ над главной
            }

            a[0] = 0;
            for (int m = 1; m < M; m++)
            {
                a[m] = - tau / (hy * hy) * kMinus[m]; // Диагональ под главной
            }

            solveMatrix(M, a, c, b, intermediateLayer[l], nextLayer[l]);        
        }

        if (compare(nextLayer, currentLayer, L, M) < epsilon)
        {
            break;
        }

        for (int l = 0; l < L; l++)
        {
            for (int m = 0; m < M; m++)
            {
                currentLayer[m][l] = nextLayer[m][l];
            }
        }  
    }
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
        //u[i+1][0][0] = nextIteration();
        int a = 1;
    }

    for (int l = 0; l < L; l++)
    {
        delete[] grid[l];
    }
    
    delete[] grid;

    return 0;
}