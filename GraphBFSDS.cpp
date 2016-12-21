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
          int n, bool *&isVisited, int *& parent, int *& rankN)
{
    vector<pair<int, int> > nextGen;
    bool stop = false;
    int s;
    vector<int> nodes;
    if (worldRank == 0)
    {
        for (int k = 0; k < n; k++)
        {
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
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast(&stop, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast(&s,1 , MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Bcast((void *)nodes.data(), nodes.size(), MPI_INT, 0, MPI_COMM_WORLD);
                    MPI_Barrier(MPI_COMM_WORLD);
                    nextGen.clear();
                    buildNodes(nodes, adjNode, numberRow, startRow, n, isVisited, nextGen, worldRank);
                    nodes.clear();
                    for (int p = 0; p < worldSize; p++)
                    {
                        if (p != 0)
                        {
                            MPI_Status status;
                            MPI_Probe(p, 0, MPI_COMM_WORLD, &status);
                            int countB;
                            MPI_Get_count(&status, MPI_BYTE, &countB);
                            int tmpC = countB/sizeof(pair<int, int>);
                            nextGen.resize(tmpC);
                            MPI_Recv((void *) nextGen.data(), countB, MPI_BYTE, p, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
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
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&stop, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
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

void printResult(int *& parent, int *& rankN, int n)
{
    FILE * g = fopen("result.out", "w");
    for (int i = 0; i < n; i++)
    {
        fprintf(g,"%d %d\n", i, rankN[i]);
    }
    fclose(g);
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
    int n, m;
    if (worldRank == 0)
    {
        convert << worldRank;
        tam = "partFile" + convert.str() + ".inp";
        fileName = tam.c_str();
        f = fopen(fileName, "r");
        fscanf(f, "%d %d", &n, &m);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    int startVertex = (n / worldSize) * worldRank;
    int numberVertex = n / worldSize;
    if (worldRank == worldSize - 1)
        numberVertex += n % worldSize;
    adjNode = new vector<int>[numberVertex];
    int rowNow = 0;
    for (int i = startFile; i <= finishFile; i++)
    {
        if (i != 0)
        {
            ostringstream convert;
            convert << i;
            tam = "partFile" + convert.str() + ".inp";
            fileName = tam.c_str();
            f = fopen(fileName, "r");
        }
        int rowForThisFile = n / numberFile;
        if (i == numberFile - 1)
            rowForThisFile += n % numberFile;
        int u, v;
        pair<int, int> tmp;
        int col = 0;
        for (int k = 0; k < rowForThisFile; k++)
        {
            for (col = 0; col < n; col++)
            {
                fscanf(f, "%d %d", &u, &v);
                tmp.first = u;
                tmp.second = v;
                MPI_Send(&tmp, sizeof(tmp), MPI_BYTE, )

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
    if (worldRank == 0)
    {
        parent  = new int[n];
        rankN = new int[n];
    }
    for (int times = 0; times < 20; times++)
    {
        memset(isVisited, false, sizeof(bool) * n);
        double totalTime = 0.0;
        if (worldRank == 0)
        {
            totalTime -= MPI_Wtime();
        }
        BFS(worldRank, worldSize, adjNode, numberRow, startRow, n, isVisited, parent,rankN);
        if (worldRank == 0)
        {
            totalTime += MPI_Wtime();
            averageTime += totalTime;
        }
    }

    free(isVisited);
    if (worldRank == 0)
    {
        printResult(parent, rankN, n);
        cout <<"time = " << averageTime/20 << endl;
    }
    MPI_Finalize();
}
