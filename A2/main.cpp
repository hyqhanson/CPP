#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <stdexcept>
#include "infPrecision.hpp"
using namespace std;

int main()
{
    try
    {
        // Default constructor
        infPrecision _default = infPrecision();
        cout << "The default number is: " << _default << "\n";

        // Constructor with one string of integer
        string str1 = "-9223372036854775807000";
        infPrecision infP1 = infPrecision(str1);
        cout << "Using string to construct the number: " << infP1 << "\n";
        // If the string has non-numerical elements (except '-' in the beginning,
        // indicates this is a negative integer), it will throw an exception.
        // Example:
        /*
        string str_err1 = "1234d3";
        infPrecision strErr1 = infPrecision(str_err1);  // It will throw the exception
        */

        // If the string has leading zeros (if the integer is not 0 but 0 is the
        // left-most digit), it will throw an exception.
        // Example:
        /*
        string str_err2 = "013423";
        infPrecision strErr2 = infPrecision(str_err2); // It will throw the exception
        */

        // Constructor with one 64 bits fixed-width signed integer
        int64_t _64bit1 = 123400;
        infPrecision infP2 = infPrecision(_64bit1);
        cout << "Using fixed-width integer to construct the number: " << infP2 << "\n";

        // Constructor with one 64 bits fixed-width signed integer vector
        vector<int64_t> vec1 = {-1, 0, 0, 0, 4};
        infPrecision infP3 = infPrecision(vec1);
        cout << "Using fixed-width integer vector to construct the number: " << infP3 << "\n\n";
        // If the vector has a negative integer (other than the first integer
        // in the beginning), it will throw an exception.
        // Example:
        /*
        vector<int64_t> vec_err1 = {2, -3, 4};
        infPrecision vecErr1 = infPrecision(vec_err1); // It will throw the exception
        */

        // If the vector has an element that is not a single-digit integer,
        // it will throw an exception.
        // Example:
        /*
        vector<int64_t> vec_err2 = {2, 33, 4};
        infPrecision vecErr2 = infPrecision(vec_err2); // It will throw the exception
        */

        // If the vector has leading zeros (if the integer is not 0 but 0 is the
        // first item in the vector), it will throw an exception.
        // Example:
        /*
        vector<int64_t> vec_err3 = {0,2, 3, 4};
        infPrecision vecErr3 = infPrecision(vec_err3);  // It will throw the exception
        */

        // Some member functions
        // 1. get_vec(), return the vector of the infPrecision object
        infPrecision M1("12345");
        cout << "M1 = " << M1 << "\n";
        vector<int64_t> m1 = M1.get_vec();
        cout << "m1 contains the vector from M1: ";
        for (int64_t &i : m1)
        {
            cout << i;
        }
        cout << "\n";
        // 2. is_negative(), return true if the integer in the infPrecision
        // object is negative, false otherwise
        infPrecision M2("-12");
        cout << "Is M2 = " << M2 << " negative? " << boolalpha << (M2.is_negative()) << "\n";
        infPrecision M3("0");
        cout << "Is M3 = " << M3 << " negative? " << (M3.is_negative()) << "\n\n";

        // Additions
        // 1. Two positive integers
        string add1 = "1";
        string add2 = "9999";
        cout << "add1 = " << add1 << ", add2 = " << add2 << "\n";
        infPrecision Add1(add1);
        infPrecision Add2(add2);
        cout << "add1 + add2 = " << Add1 + Add2 << "\n";
        // 2. One negative integer and one positive integer
        string add3 = "15";
        string add4 = "-3015";
        cout << "add3 = " << add3 << ", add4 = " << add4 << "\n";
        infPrecision Add3(add3);
        infPrecision Add4(add4);
        cout << "add3 + add4 = " << Add3 + Add4 << "\n";
        // 3. Two negative integers
        int64_t add5 = -1111;
        int64_t add6 = -9999;
        cout << "add5 = " << add5 << ", add6 = " << add6 << "\n";
        infPrecision Add5(add5);
        infPrecision Add6(add6);
        cout << "add5 + add6 = " << Add5 + Add6 << "\n";
        // 4. Add and assign
        Add5 += Add6;
        cout << "add5 += add6, add5 = " << Add5 << "\n\n";

        // Subtractions
        // 1. Two positive integers
        string sub1 = "352331";
        string sub2 = "4520";
        cout << "sub1 = " << sub1 << ", sub2 = " << sub2 << "\n";
        infPrecision Sub1(sub1);
        infPrecision Sub2(sub2);
        cout << "sub1 - sub2 = " << Sub1 - Sub2 << "\n";
        // 2. One negative integer and one positive integer
        string sub3 = "200000000000000000000000000000000";
        string sub4 = "-1111";
        cout << "sub3 = " << sub3 << ", sub4 = " << sub4 << "\n";
        infPrecision Sub3(sub3);
        infPrecision Sub4(sub4);
        cout << "sub3 - sub4 = " << Sub3 - Sub4 << "\n";
        cout << "sub4 - sub3 = " << Sub4 - Sub3 << "\n";
        // 3. Two negative integers
        int64_t sub5 = -10000;
        int64_t sub6 = -9997;
        cout << "sub5 = " << sub5 << ", sub6 = " << sub6 << "\n";
        infPrecision Sub5(sub5);
        infPrecision Sub6(sub6);
        cout << "sub5 - sub6 = " << Sub5 - Sub6 << "\n";
        cout << "sub6 - sub5 = " << Sub6 - Sub5 << "\n";
        // 4. Subtract and assign
        Sub5 -= Sub6;
        cout << "sub5 -= sub6, sub5 = " << Sub5 << "\n\n";

        // Multiplication
        // 1. Two positive integers
        string mul1 = "123";
        string mul2 = "99";
        cout << "mul1 = " << mul1 << ", mul2 = " << mul2 << "\n";
        infPrecision Mul1(mul1);
        infPrecision Mul2(mul2);
        cout << "mul1 * mul2 = " << Mul1 * Mul2 << "\n";
        // 2. One negative integer and one positive integer
        string mul3 = "-99";
        string mul4 = "99";
        cout << "mul3 = " << mul3 << ", mul4 = " << mul4 << "\n";
        infPrecision Mul3(mul3);
        infPrecision Mul4(mul4);
        cout << "mul3 * mul4 = " << Mul3 * Mul4 << "\n";
        // 3. Two negative integers
        int64_t mul5 = -54321;
        int64_t mul6 = -1;
        cout << "mul5 = " << mul5 << ", mul6 = " << mul6 << "\n";
        infPrecision Mul5(mul5);
        infPrecision Mul6(mul6);
        cout << "mul5 * mul6 = " << Mul5 * Mul6 << "\n";
        // 4. Multiply and assign
        Mul5 *= Mul6;
        cout << "mul5 *= mul6, mul5 = " << Mul5 << "\n\n";

        // Negation
        // 1. Negative testing
        infPrecision N1("100");
        cout << "n1 = " << N1 << "\n";
        cout << "Is n1 negative? " << boolalpha << N1.is_negative() << "\n";
        cout << "-n1 = " << -N1 << "\n";
        cout << "Is -n1 negative? " << (-N1).is_negative() << "\n";
        // 2. Combined calculation using negation
        string n2 = "-2";
        string n3 = "125";
        string n4 = "-5";
        cout << "n2 = " << n2 << ", n3 = " << n3 << ", n4 = " << n4 << "\n";
        infPrecision N2(n2);
        infPrecision N3(n3);
        infPrecision N4(n4);
        cout << "((-n2) - n3) * (-n4) = " << ((-N2) - N3) * (-N4) << "\n\n";

        // Increment & decrement
        // 1. ++ (prefix)
        infPrecision I1("-1");
        cout << "i1 = " << I1 << "\n";
        cout << "++i1 = " << ++I1 << "\n";
        // 2. ++ (postfix)
        infPrecision I2("100000000");
        cout << "i2 = " << I2 << "\n";
        cout << "i2++ = " << I2++ << ", i2 = " << I2 << "\n";
        // 3. -- (prefix)
        infPrecision D1("0");
        cout << "d1 = " << D1 << "\n";
        cout << "--d1 = " << --D1 << "\n";
        // 4. -- (postfix)
        infPrecision D2("10000000");
        cout << "d2 = " << D2 << "\n";
        cout << "d2-- = " << D2-- << ", d2 = " << D2 << "\n\n";

        // Comparisons
        // 1. == & !=
        string c1 = "100000199999000";
        string c2 = "100000199999000";
        string c3 = "100000019999900";
        cout << "c1 = " << c1 << ", c2 = " << c2 << ", c3 = " << c3 << "\n";
        infPrecision C1(c1);
        infPrecision C2(c2);
        infPrecision C3(c3);
        cout << "c1 == c2 is " << (c1 == c2) << "\n";
        cout << "c1 == c3 is " << (c1 == c3) << "\n";
        cout << "c1 != c3 is " << (c1 != c3) << "\n";
        cout << "c3 != c3 is " << (c3 != c3) << "\n";

        // 2. < & <=
        infPrecision C7("0");
        infPrecision C8("1");
        infPrecision C9("1");
        cout << "c7 = " << C7 << ", c8 = " << C8 << ", c9 = " << C9 << "\n";
        cout << "c7 < c8 is " << (C7 < C8) << "\n";
        cout << "c8 < c7 is " << (C8 < C7) << "\n";
        cout << "c8 <= c9 is " << (C8 <= C9) << "\n";
        cout << "c9 <= c7 is " << (C9 <= C7) << "\n\n";

        // 3. > & >=
        infPrecision C4("-999");
        infPrecision C5("-1000");
        infPrecision C6("-999");
        cout << "c4 = " << C4 << ", c5 = " << C5 << ", c6 = " << C6 << "\n";
        cout << "c4 > c5 is " << (C4 > C5) << "\n";
        cout << "c5 > c4 is " << (C5 > C4) << "\n";
        cout << "c4 >= c6 is " << (C4 >= C6) << "\n";
        cout << "c5 >= c6 is " << (C5 >= C6) << "\n\n";

        // Assignment
        infPrecision A1("1234");
        infPrecision A2("4321");
        cout << "A1 = " << A1 << ", A2 = " << A2 << "\n";
        A2 = A1;
        cout << "After assign A2 = A1, A2 = " << A2 << "\n\n";

        // Insertion
        infPrecision In1("1111000000000000000000000000000");
        cout << "In1 = " << In1 << "\n";
    }

    catch (const infPrecision::not_numerical &e)
    {
        cout << "The string contains a character that is not numerical!\n";
    }
    catch (const infPrecision::cannot_negative &e)
    {
        cout << "There should only be negative number in the front of the original vector!\n";
    }
    catch (const infPrecision::cannot_mulDigits &e)
    {
        cout << "There should only be single digit number in the input vector!\n";
    }
    catch (const infPrecision::cannot_leadZero &e)
    {
        cout << "Non-zero Integers should not have leading zeros!\n";
    }
}