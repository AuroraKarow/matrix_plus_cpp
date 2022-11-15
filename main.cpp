// This is a c++ 20 standard file.
#pragma once

#define MATRIX_UNROLL    4
#define MATRIX_BLOCKSIZE 32

#define OMP_MATRIX_MODE  false

#include <iostream>
#include "matrix"

using std::cout;
using std::endl;

using neunet::vect;
using neunet::net_sequence;

int main(int argc, char *argv[], char *envp[]) {
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
          h = {{0, 1, 2}},
          i = {{a, {{b, {{c},
                         {d}}},
                    {e,   f}}},
               {g,        h}};
     net_sequence<vect> seq = {a, b, c, d, e, f, g, h, i};
     for (auto temp : seq) {
          cout << temp << '\n';
          cout << endl;
     }
     return EXIT_SUCCESS;
}