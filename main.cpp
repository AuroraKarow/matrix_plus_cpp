// This is a c++ 17 standard file.

#include <iostream>
#include "matrix"

using std::cout;
using std::endl;

using mtx::vect;

int main(int argc, char *argv[], char *envp[])
{
    cout << "hello, world." << endl;
    vect a = {{0},
              {1},
              {2},
              {3}},
         b = {{0, 1, 2},
              {3, 4, 5},
              {6, 7, 8}},
         c = {{0, 1},
              {2, 3}},
         d = {{0, 1}},
         e = {{0, 1, 2, 3}},
         f = {{0}},
         g = {{0, 1, 2}},
         h = {{0, 1, 2}};
    vect i = {{a, {{b, {{c},
                        {d}}},
                    {e,  f}}},
              {g,        h}};
    cout << i << endl;
    return EXIT_SUCCESS;
}