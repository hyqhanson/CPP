#include <iostream>
#include <stdexcept>
#include "inf_precision_Ped.hpp"
#include <string>

using namespace std;
int main()
{
     try
     {
          string str_Ped = "+000065138464834800098656454517931564602317"; // 9835 -06517  +06517  a06517 06517 +00006513846483480009861202317
          inf_precision_Ped numb1(str_Ped);
          cout << "Num#1: " << numb1 << "\n";

          // inf_precision_Ped numb1 = inf_precision_Ped(str_Ped);
          string str_Ped2 = "-6236487630129"; // 9935 152468
          inf_precision_Ped numb2(str_Ped2);
          cout << "Num#2: " << numb2 << "\n";

          cout << "Num#1 + Num#2= " << numb1 + numb2 << "\n";
          cout << "Num#1 - Num#2= " << numb1 - numb2 << "\n";
          cout << "Num#1 * Num#2= " << numb1 * numb2 << "\n";

          cout << "++Num#1 = " << ++numb1 << "\n";
          cout << "Num#1++ = " << numb1++ << "\n";

          cout << "--Num#1 = " << --numb1 << "\n";
          cout << "Num#1-- = " << numb1-- << "\n";

          bool ped = numb1 > numb2;
          cout << "Num#1"
               << " > "
               << "Num#2"
               << "?\n"
               << boolalpha << ped << '\n';
          cout << "Num#1"
               << " < "
               << "Num#2"
               << "?\n"
               << boolalpha << (ped = numb1 < numb2) << '\n';
          cout << "Num#1"
               << " == "
               << "Num#2"
               << "?\n"
               << boolalpha << (ped = numb1 == numb2) << '\n';
          cout << "Num#1"
               << " >= "
               << "Num#2"
               << "?\n"
               << boolalpha << (ped = numb1 >= numb2) << '\n';
          cout << "Num#1"
               << " <= "
               << "Num#2"
               << "?\n"
               << boolalpha << (ped = numb1 <= numb2) << '\n';
          cout << "Num#1"
               << " != "
               << "Num#2"
               << "?\n"
               << boolalpha << (ped = numb1 != numb2) << '\n';

          cout << "Num#1 *= Num#2= " << (numb1 *= numb2) << "\n";
          // cout << "Num#1: " << numb1 << "\n";
          cout << "Num#1 -= Num#2= " << (numb1 -= numb2) << "\n";
          // cout << "Num#1: " << numb1 << "\n";
          cout << "Num#1 += Num#2= " << (numb1 += numb2) << "\n";
          // cout << "Num#1: " << numb1 << "\n";
          cout << "-Num#2= " << -numb2 << "\n";
     }
     catch (const exception &e)
     {
          cout << "Error: " << e.what() << '\n';
     }
}