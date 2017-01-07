#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string.h>
#include <vector>
#include <sstream>

using namespace std;

void buildNodes(vector<int> & inputNode, vector<int> * &adjNode, int numberRow, int startRow, int n,
                bool *& isVisited, vector<pair<int,int> > &result, int worldRank)
{
    for (int k = 0; k < inputNode.size(); k++)
        isVisited[inputNode[k]] = true;
    int finishRow = startRow + numberRow - 1;
    int u, v, nodes;
    for (int k =0; k < inputNode.size(); k++)
    {
        u = inputNode[k];
        if (startRow <= u && u <= finishRow)
        {
            v = u - startRow;
            for (int i = 0; i < adjNode[v].size(); i++)
            {
                nodes = adjNode[v][i];
                if (!isVisited[nodes])
                {
                    result.push_back(make_pair(nodes, u));
                }
            }
        }
    }
    /*
    for (int k = 0; k < result.size(); k++)
        printf("result %d %d %d\n", worldRank, result[k].first, result[k].second);
        */

}
void BFS(int worldRank, int worldSize, vector<int> *& adjNode, int numberRow, int startRow,
          int n, bool *&isVisited, int *& parent, int *& rankN, double & communicationTime)
{
    vector<pair<int, int> > nextGen;
    bool stop = false;
    int s;
    vector<int> nodes;
    if (worldRank == 0)
    {
        for (int k = 0; k < n; k++)
        {
            if (k % 1000 == 0)
                printf("vertex: %d\n", k);
            if (!isVisited[k])
            {
                parent[k] = -1;
                isVisited[k] = true;
                bool isFind = true;
                nodes.push_back(k);
                rankN[k] = 0;
                while (isFind)
                {
                    isFind = false;
                    s = nodes.size();
                    communicationTime -= MPI_Wtime();
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast(&stop, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast(&s,1 , MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast((void *)nodes.data(), nodes.size(), MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    communicationTime += MPI_Wtime();
                    nextGen.clear();
                    buildNodes(nodes, adjNode, numberRow, startRow, n, isVisited, nextGen, worldRank);
                    nodes.clear();
                    for (int p = 0; p < worldSize; p++)
                    {
                        if (p != 0)
                        {
                            MPI_Status status;
                            communicationTime -= MPI_Wtime();
                            MPI_Probe(p, 0, MPI_COMM_WORLD, &status);
                            communicationTime += MPI_Wtime();
                            int countB;
                            MPI_Get_count(&status, MPI_BYTE, &countB);
                            int tmpC = countB/sizeof(pair<int, int>);
                            nextGen.resize(tmpC);
                            communicationTime -= MPI_Wtime();
                            MPI_Recv((void *) nextGen.data(), countB, MPI_BYTE, p, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                            communicationTime += MPI_Wtime();
                        }
                        int t = nextGen.size();
                        //printf("check %d %d\n", p, t);
                        for (int h = 0; h < nextGen.size(); h++)
                        {
                            int u = nextGen[h].first;
                            int v = nextGen[h].second;
                            if (!isVisited[u])
                            {
                                isVisited[u] = true;
                                parent[u] = v;
                                rankN[u] = rankN[v] + 1;
                                nodes.push_back(u);
                            }
                        }
                    }
                    if (nodes.size() != 0)
                    isFind = true;
                }

            }
        }
        stop = true;
        communicationTime -= MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&stop, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        communicationTime += MPI_Wtime();
    }
    else
    {
        while (!stop)
        {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&stop, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            if (stop)
                break;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&s,1 , MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            nodes.resize(s);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast((void *)nodes.data(), s, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            /*printf("world rank %d: %d\n", worldRank, s);
            for (int k = 0; k < s; k++)
                printf("world rank %d: %d\n", worldRank, nodes[k]);*/

            nextGen.clear();
            buildNodes(nodes, adjNode, numberRow, startRow, n, isVisited, nextGen, worldRank);
            MPI_Send((void*)nextGen.data(), sizeof(pair<int,int>) * nextGen.size(),
                     MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

void printResult(int *& parent, int *& rankN, int n, int worldSize, double averageTime, double communicationTime)
{
    ostringstream convert;
    convert << worldSize;
    string tam = "resultBFSOfNumberProcess" + convert.str() + ".out";
    const char * fileName1 = tam.c_str();
    tam = "resultBFSTimeOfProcess" + convert.str() + ".out";
    const char * fileName2 = tam.c_str();
    FILE * f1 = fopen(fileName1,"w");
    FILE * f2 = fopen(fileName2,"w");
    for (int i = 0; i < n; i++)
    {
        fprintf(f1,"%d %d %d\n", i, rankN[i], parent[i]);
    };
    fprintf(f2, "totaltime: %f\n", averageTime);
    fprintf(f2, "commnunacationTime: %f\n", communicationTime);
    fclose(f1);
    fclose(f2);
}


int main(int argc, char ** argv)
{
    int numberFile = atoi(argv[1]);
    vector<int> * adjNode;
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
    double timeReading = -MPI_Wtime();
    int n;
    if (worldRank == 0)
    {
        convert << worldRank;
        tam = "partFileBFS" + convert.str() + ".inp";
        fileName = tam.c_str();
        f = fopen(fileName, "r");
        fscanf(f, "%d", &n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int startRow = (n / numberFile) * startFile;
    int numberRow = (n / numberFile) * (finishFile - startFile + 1);
    if (worldRank == worldSize - 1)
        numberRow += n % numberFile;
    adjNode = new vector<int>[numberRow];
    for (int i = 0; i < numberRow; i++)
        adjNode[i].clear();
    int rowNow = 0;
    for (int i = startFile; i <= finishFile; i++)
    {
        if (i != 0)
        {
            ostringstream convert;
            convert << i;
            tam = "partFileBFS" + convert.str() + ".inp";
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
                if (u == 1)
                {
                    adjNode[rowNow].push_back(col);
                }
            }
            rowNow++;
        }
        fclose(f);
    }
    /*
    for (int i = 0; i < numberRow; i++)
    {
        for (int j = 0; j < adjNode[i].size(); j++)
            printf("%d %d %d\n", worldRank, i + startRow, adjNode[i][j]);
    }
*/

    //BFS code

    bool * isVisited = new bool[n];
    int * parent;
    int * rankN;
    double averageTime = 0;
    double commnucationTime = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    if (worldRank == 0)
    {

        timeReading += MPI_Wtime();
        printf("time reading: %f\n", timeReading);
        parent  = new int[n];
        rankN = new int[n];
    }
    for (int times = 0; times < 1; times++)
    {
        memset(isVisited, false, sizeof(bool) * n);
        double totalTime = 0.0;
        if (worldRank == 0)
        {
            totalTime -= MPI_Wtime();
        }
        BFS(worldRank, worldSize, adjNode, numberRow, startRow, n, isVisited, parent,rankN, commnucationTime);
        if (worldRank == 0)
        {
            totalTime += MPI_Wtime();
            averageTime += totalTime;
        }
    }
    if (worldRank == 0)
        printf("ok\n");
    free(isVisited);
    MPI_Barrier(MPI_COMM_WORLD);
    if (worldRank == 0)
    {
        printResult(parent, rankN, n, worldSize, averageTime, commnucationTime);
        free(parent);
        free(rankN);
    }
    else
    {
        for (int i = 0; i < numberRow; i++)
            adjNode[i].clear();
        //free(adjNode);
        //free(isVisited);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
