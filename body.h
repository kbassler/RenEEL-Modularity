#define Inf_Sma 1e-6

struct graph
{
    int N;
    int n;
    int com;
    int **ml;
    int **nl;
    int **wl;
    int *dl;
};

struct part
{
    double Q;
    int com;
    int *pa;
};