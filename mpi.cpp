#include <iostream>
#include <stdlib.h>
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[])
    {
        int numprocs, myid;
        int *a, *b, *c;
        int n, p, l;
        int N=10, M=10;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Status status;

        if (myid == 0) 
        {
            a = new int[N*M]; 
            b = new int[M];
            c = new int[N];

            for (int i=0; i<M; i++)
            {
                b[i] = 1;
            }

            for (int j=0; j<M*N; j++)
            {
                a[j] = 1;
            }

            n = int(N/numprocs);

            for (int k = 1; k<numprocs; k++) 
            {
                MPI_Send(&n, 1, MPI_INT, k, 0, MPI_COMM_WORLD); 
                MPI_Send(&a[k*n*M], n*M, MPI_INT, k, 1, MPI_COMM_WORLD);
                MPI_Send(&b[0], M, MPI_INT, k, 2, MPI_COMM_WORLD); 
            }
        }
        else 
        {
            MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            a = new int[n*M];
            b = new int[M];
            c = new int[n];

            MPI_Recv(a, n*M, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
            MPI_Recv(b, M, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier(MPI_COMM_WORLD); 
        
        for (int i=0; i<n; i++)
        {
            c[i] = 0; 
        }
	
	double t = MPI_Wtime();

        MPI_Barrier(MPI_COMM_WORLD);

        for (int i=0; i<n; i++)
        {
            for (int j=0; j<M; j++)
            {
                c[i] += a[j+i*n] * b[j];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
 
        if (myid==0) 
        {
            for (int i=1; i<numprocs; i++)
            {
                MPI_Recv(&c[i*n], n, MPI_INT, i, 3, MPI_COMM_WORLD, &status);
            } 
        }
        else 
        {
            MPI_Send(c, n, MPI_INT, 0, 3, MPI_COMM_WORLD);
        }
    
        if (myid==0)
        {
            t = MPI_Wtime() - t;
            cout << "Time: " << t << "\n";
        }

        return 0;
}

