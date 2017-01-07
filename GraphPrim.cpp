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
    int indexAndParent;
};

void Prim(int worldRank, int n, dis *& distant, bool *&boo, int startRow, int numberRow,int ** cost,
          int *& resultEdgeA, int *& resultEdgeB, float *& resultEdgeC, double & totalResult, double & communacationTime)
{
    //printf("costt[0][1] = %d\n ", cost[0][1]);
    for (int k = 0; k < n; k++)
    {
        if (k % 1000 == 0)
            printf("vertex: %d %d\n", worldRank, k);
        //printf("check dis[5]: %f %d\n", distant[5].val, distant[5].index);
        dis minDist;
        minDist.val = oo;
        int finishRow = startRow + numberRow - 1;
        for (int i = 0; i < numberRow; i++)
        {
            if (minDist.val > distant[i].val && !boo[i])
            {
                minDist = distant[i];
            }
        }
        dis result;
        communacationTime -= MPI_Wtime();
        MPI_Allreduce(&minDist, &result, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        communacationTime += MPI_Wtime();
        int resultIndex = result.indexAndParent / n;
        int resultParent = result.indexAndParent % n;
        totalResult = totalResult + result.val;
        //printf("%d min: %d %d %f %d\n", worldRank, k, resultIndex, result.val, resultParent);
        if (worldRank == 0)
        {
            if (k != 0)
            {
                resultEdgeA[k - 1] = resultIndex;
                resultEdgeB[k - 1] = resultParent;
                resultEdgeC[k - 1] = result.val;
            }
        }
        int vertex;
        for (int i = 0; i < numberRow; i++)
        {
            vertex = distant[i].indexAndParent / n;
            if (distant[i].val > cost[i][resultIndex])
            {
                //printf("check: %d %d %f cost[%d][%d] = %d\n",worldRank, vertex, distant[i].val, i, resultIndex, cost[i][resultIndex]);
                distant[i].val = cost[i][resultIndex];
                distant[i].indexAndParent = vertex * n + resultIndex;
            }
            //printf("check: %d %d %f\n",worldRank, vertex, distant[i].val);
            if (vertex == resultIndex)
                boo[i] = true;
        }

    }
}

void printResult(double result, int *& resultEdgeA, int *& resultEdgeB ,
                 float *& resultEdgeC, int n, int worldSize, double averageTime, double communicationTime)
{
    ostringstream convert;
    convert << worldSize;
    string tam = "resultPrimOfNumberProcess" + convert.str() + ".out";
    const char * fileName1 = tam.c_str();
    tam = "resultPrimTimeOfProcess" + convert.str() + ".out";
    const char * fileName2 = tam.c_str();
    FILE * f1 = fopen(fileName1,"w");
    FILE * f2 = fopen(fileName2,"w");
    fprintf(f1, "The sum: %f\n", result);
    fprintf(f1, "The list edges: \n");
    for (int i = 0; i < n - 1; i++)
        fprintf(f1, "%d %d %f\n", resultEdgeA[i], resultEdgeB[i], resultEdgeC[i]);
    fprintf(f2, "totaltime: %f\n", averageTime);
    fprintf(f2, "commnunacationTime: %f\n", communicationTime);
    fclose(f1);
    fclose(f2);

}

int main(int argc, char ** argv)
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
    dis * distant;
    bool * boo;
    double * finalResult;
    int * parent;
    double timeReading = -MPI_Wtime();
    if (worldRank == 0)
    {
        convert << worldRank;
        tam = "partFileDijsktraP" + convert.str() + ".inp";
        //cout << tam << endl;
        fileName = tam.c_str();
        f = fopen(fileName, "r");
        fscanf(f, "%d", &n);
        int s;
        fscanf(f, "%d", &s);
        finalResult = new double[n];
    }
   // cout << " ok " << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
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
            tam = "partFileDijsktraP" + convert.str() + ".inp";
            fileName = tam.c_str();
            f = fopen(fileName, "r");
        }
        //cout << "ok" << endl;
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
                if (u <= 0)
                    u = oo;
                adjNode[rowNow][col] = u;
                //printf("%d cost[%d][%d] = %d\n",worldRank, rowNow,col, adjNode[rowNow][col]);
            }
            rowNow++;
        }
        fclose(f);
    }
    double averageTime = 0;
    double commnucationTime = 0.0;
    int * resultEdgeA;
    int * resultEdgeB;
    float * resultEdgeC;
    double totalResult = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    if (worldRank == 0)
    {
        timeReading += MPI_Wtime();
        printf("time reading: %f\n", timeReading);
    }
    for (int times = 0; times < 1; times++)
    {
        double totalTime = 0.0;
        if (worldRank == 0)
        {
            totalTime -= MPI_Wtime();
            resultEdgeA = new int[n];
            resultEdgeB = new int[n];
            resultEdgeC = new float[n];
        }
        for (int i = 0; i < numberRow; i++)
        {
            distant[i].val = oo;
            distant[i].indexAndParent = (startRow + i) * n;
        //printf("%d %f %d\n", worldRank, distant[i].val, distant[i].index);
            if (worldRank == 0 && i == 0)
            {
                distant[i].val = 0;
            }
        }
        memset(boo, false, numberRow * sizeof(bool));
        MPI_Barrier(MPI_COMM_WORLD);
        Prim(worldRank, n, distant, boo, startRow, numberRow, adjNode,
              resultEdgeA, resultEdgeB, resultEdgeC, totalResult, commnucationTime);
        MPI_Barrier(MPI_COMM_WORLD);
        totalTime += MPI_Wtime();
        averageTime += totalTime;
    }
    if (worldRank == 0)
    {
        printResult(totalResult, resultEdgeA, resultEdgeB, resultEdgeC, n, worldSize, averageTime, commnucationTime);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(distant);
    free(boo);
    for (int i = 0; i < numberRow; i++)
        free(adjNode[i]);
    free(adjNode);
    if (worldRank == 0)
    {
        free(resultEdgeA);
        free(resultEdgeB);
        free(resultEdgeC);
    }

    MPI_Finalize();
    //printf("costs[0][1] = %d\n", adjNode[0][1]);

}

