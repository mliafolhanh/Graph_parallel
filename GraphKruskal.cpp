#include <cstdio>
#include <iostream>
#include <mpi.h>
#include <string.h>
#include <cstdlib>
#include <sstream>

using namespace std;
const int oo = 1000000007;

struct dis{
    double val;
    int index;
    int parent;
};

void Prim(int worldRank, int n, int s, dis *& distant, bool *&boo, int startRow, int numberRow,
          int ** cost, double *& finalResult, int *& resultEdgeA, int *& resultEdgeB)
{
    //printf("costt[0][1] = %d\n ", cost[0][1]);
    dis minDist;
    for (int k = 0; k < n; k++)
    {
        //printf("check dis[5]: %f %d\n", distant[5].val, distant[5].index);
        MPI_Barrier(MPI_COMM_WORLD);
        minDist.val = oo;
        int finishRow = startRow + numberRow - 1;
        for (int i = 0; i < numberRow; i++)
        {
            if (minDist.val > distant[i].val && !boo[i])
            {
                minDist = distant[i];
            }
        }
       // printf("%d min: %d %f\n", worldRank, minDist.index, minDist.val);
        dis result;
        MPI_Allreduce(&minDist, &result, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (worldRank == 0)
        {
            if (k != 0)
            {
                resultEdgeA[k - 1] = result.index;
                resultEdgeB[k - 1] = result.parent;
            }
        }
        int vertex;
        for (int i = 0; i < numberRow; i++)
        {
            vertex = distant[i].index;
            if (distant[i].val > cost[i][result.index])
            {
                distant[i].val = cost[i][result.index];
                distant[i].parent = result.index;
            }

            if (vertex == result.index)
                boo[i] = true;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void printResult(double * &distant, int n)
{
    FILE * f = fopen("outputPrim.out","w");
    for (int i = 0; i < n; i++)
        fprintf(f, "%d %f\n", i, distant[i]);
    fclose(f);

}

int main(int argc, char ** argv)8
{
    int numberFile = atoi(argv[1]);
    int worldRank, worldSize;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    int startFile = (numberFile/worldSize) * worldRank;
    int finishFile = (numberFile/worldSize) + startFile - 1;
    if (worldRank == worldSize - 1)
        finishFile = numberFile - 1;
    ostringstream convert;
    const char * fileName;
    string tam;
    FILE *f = NULL;
    int n;
    int s;
    dis * distant;
    bool * boo;
    double * finalResult;
    int * parent;
    if (worldRank == 0)
    {
        convert << worldRank;
        tam = "partFileDijkstra" + convert.str() + ".inp";
        fileName = tam.c_str();
        f = fopen(fileName, "r");
        fscanf(f, "%d %d", &n, &s);
        finalResult = new double[n];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&s, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int startRow = (n / numberFile) * startFile;
    int numberRow = (n / numberFile) * (finishFile - startFile + 1);
    if (worldRank == worldSize - 1)
        numberRow += n % numberFile;
    int ** adjNode;
    adjNode = new int*[numberRow];
    for (int i = 0; i < numberRow; i++)
        adjNode[i] = new int[n];
    distant = new dis[numberRow];
    boo = new bool[numberRow];
    int rowNow = 0;
    for (int i = startFile; i <= finishFile; i++)
    {
        if (i != 0)
        {
            ostringstream convert;
            convert << i;
            tam = "partFileDijkstra" + convert.str() + ".inp";
            fileName = tam.c_str();
            f = fopen(fileName, "r");
        }
        int rowForThisFile = n / numberFile;
        if (i == numberFile - 1)
            rowForThisFile += n % numberFile;
        int u;
        int col = 0;
        for (int k = 0; k < rowForThisFile; k++)
        {
            for (col = 0; col < n; col++)
            {
                fscanf(f, "%d", &u);
                if (u < 0)
                    u = oo;
                adjNode[rowNow][col] = u;
                //printf("%d cost[%d][%d] = %d\n",worldRank, rowNow,col, adjNode[rowNow][col]);
            }
            rowNow++;
        }
        fclose(f);
    }
    double averageTime = 0;
    int * resultEdgeA;
    int * resultEdgeB;
    for (int times = 0; times < 1; times++)
    {
        double totalTime = 0.0;
        if (worldRank == 0)
        {
            totalTime -= MPI_Wtime();
            resultEdgeA = new int[n];
            resultEdgeB = new int[n];
        }
        for (int i = 0; i < numberRow; i++)
        {
            distant[i].val = oo;
            distant[i].index = startRow + i;
        //printf("%d %f %d\n", worldRank, distant[i].val, distant[i].index);
            if (startRow + i == s)
            {
                distant[i].val = 0;
            }
        }
        memset(boo, false, numberRow * sizeof(bool));
        MPI_Barrier(MPI_COMM_WORLD);
        Dijkstra(worldRank, n, s, distant, boo, startRow, numberRow, adjNode, finalResult, resultEdgeA, resultEdgeB);
        MPI_Barrier(MPI_COMM_WORLD);
        totalTime += MPI_Wtime();
        averageTime += totalTime;
    }
    if (worldRank == 0)
    {

        printResult(finalResult, n);
        printf("time: %f\n", averageTime/1);
    }
    //printf("costs[0][1] = %d\n", adjNode[0][1]);

}

                                                                          ĵ
