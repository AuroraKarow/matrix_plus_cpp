#include "matrix"

using namespace std;
using namespace mtx;
using namespace bagrt;

int main(int argc, char *argv[], char *envp[])
{
    cout << "hello, world." << endl;
    vect a = {{ 0,-1},
              { 1, 0}};
    matrix b = {{ 2, 1},
                { 8, 7}};
    matrix c = {{1, 1, 1, 0, 1, 1, 2, 0},
                {1, 1, 1, 1, 0, 1, 1, 0},
                {2, 2, 2, 1, 1, 2, 3, 1},
                {3, 3, 3, 2, 1, 3, 4, 1}};
    matrix d = {{1},
                {2},
                {3},
                {4}};
    return EXIT_SUCCESS;
}