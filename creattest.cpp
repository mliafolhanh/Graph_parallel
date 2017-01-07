#include <cstdio>
#include <iostream>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <time.h>

using namespace std;

int a[30001][30001];

int main(int argc, char ** argv)
{
    int n = 30000;
    srand(time(NULL));
    int s = rand() %n;
    for (int i = 0; i < n; i++)
    {
        a[i][i] = 0;
        for (int j = i + 1; j < n; j++)
        {
            a[i][j]  = rand() % 1000 + 50;
            a[j][i] = a[i][j];
        }
    }
    cout <<"ok"<<endl;
    int numberFile = atoi(argv[1]);
    FILE * f;
    string tam;
    const char * fileName;

    for (int k = 0; k < numberFile; k++)
    {
        ostringstream convert;
        convert << k;
        tam = "partFileDijsktraP" + convert.str() + ".inp";
        fileName = tam.c_str();
        cout << fileName << endl;
        f = fopen(fileName, "w");
        int startRow = (n / numberFile) * k;
        int finishRow = startRow + (n/numberFile) - 1;
        if (k == numberFile - 1)
            finishRow = n - 1;
        if (k == 0)
        {
            fprintf(f, "%d %d\n", n, s);
        }
        for (int i = startRow; i <= finishRow; i++)
        {
            for (int j = 0; j < n; j++)
                fprintf(f, "%d ", a[i][j]);
            fprintf(f, "\n");
        }
        fclose(f);
    }

}
