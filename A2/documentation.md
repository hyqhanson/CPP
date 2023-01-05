## Main Goal
The C++ program *main.cpp* are focusing on demonstrating the class **infPrecision** which is storing and operating infinite-precision integer objects. All declarations and implementations are in *infPrecision.hpp* header file. 
<br/>
The implementations includes:
1. [Constructor without any argument](#default-constructor)
2. [Constructor by numerical string](#constructor-with-a-string)
3. [Constructor by 64 fixed-width signed integer](#constructor-with-an-integer)
4. [Constructor by 64 fixed-width signed integer vector](#constructor-with-a-vector)
5. [Addition(+ and +=)](#addition)
6. [Subtraction (- and -=)](#subtraction)
7. [Multiplication (* and *=)](#multiplication)
8. [Negation (unary -)](#negation-unary)
9. [Increment (++) and decrement (--)](#prefix_increment)
10. [Comparison (==, !=, <, >, <=, and >=)](#equality)
11. [Assignment (=)](#assignment)
12. [Insertion (<<)](#insertion)

<br/>

## Environment
This program contains several C++ header files: 

`<iostream>`, `<stdexcept>`, `<vector>`, `<string>`, `<algorithm>`

<br/>

## class infPrecision
Class *infPrecision* has a private variable `vector<int64_t> integer;`. Every initialized infPrecision object will convert input variable in this vector, each vector element contains one single-digit integer.
### Constructors:
1. ### Default constructor
If a *infPrecision* object is initialized without any argument, it will be constructed as a integer 0.
```cpp
// Default constructor
infPrecision _default = infPrecision();
cout << "The default number is: " << _default << "\n";
```
```
The default number is: 0
```
2. ### Constructor with a string:
The argument must be a string contains only numerical characters (0~9). The heading character can contain a single negation '-'.
```cpp
// Constructor with one string of integer
string str1 = "-9223372036854775807000";
infPrecision infP1 = infPrecision(str1);
cout << "Using string to construct the number: " << infP1 << "\n";
```
```
Using string to construct the number: -9223372036854775807000
```
If the string has non-numerical elements (except '-' in the beginning, indicates this is a negative integer), it will throw an exception.
```cpp
string str_err1 = "1234d3";
infPrecision strErr1 = infPrecision(str_err1);  // It will throw the exception
```
```
The string contains a character that is not numerical!
```
If the string has leading zeros (if the integer is not 0 but 0 is the left-most digit), it will throw an exception.
```cpp
string str_err2 = "013423";
infPrecision strErr2 = infPrecision(str_err2);  // It will throw the exception
```
```
Non-zero Integers should not have leading zeros!
```
3. ### Constructor with an integer:
The argument is a 64 bit signed integer. 
```cpp
// Constructor with one 64 bits fixed-width signed integer
int64_t _64bit1 = 123400;
infPrecision infP2 = infPrecision(_64bit1);
cout << "Using fixed-width integer to construct the number: " << infP2 << "\n";
``` 
```
Using fixed-width integer to construct the number: 123400
```
4. ### Constructor with a vector:
The argument must be 64 bits signed integer vector. Each element in this vector must be a single-digit integer, and must be non-negative (other than the first element). 
```cpp
// Constructor with one 64 bits fixed-width signed integer vector
vector<int64_t> vec1 = {-1, 0, 0, 0, 4};
infPrecision infP3 = infPrecision(vec1);
cout << "Using fixed-width integer vector to construct the number: " << infP3 << "\n\n";
```
```
Using fixed-width integer vector to construct the number: -10004
```
If the vector has a negative integer (other than the first integer in the beginning), it will throw an exception:
```cpp
vector<int64_t> vec_err1 = {2, -3, 4};
infPrecision vecErr1 = infPrecision(vec_err1); // It will throw the exception
```
```
There should only be negative number in the front of the original vector!
```
If the vector has an element that is not a single-digit integer, it will throw an exception.
```cpp
vector<int64_t> vec_err2 = {2, 33, 4};
infPrecision vecErr2 = infPrecision(vec_err2); // It will throw the exception
```
```
There should only be single digit number in the input vector!
```
If the vector has leading zeros (if the integer is not 0 but 0 is the first item in the vector), it will throw an exception.
```cpp
vector<int64_t> vec_err3 = {0,2, 3, 4};
infPrecision vecErr3 = infPrecision(vec_err3);  // It will throw the exception
```
```
Non-zero Integers should not have leading zeros!
```

### Overloaded operators as member function:
1. ### Prefix Increment
Overloaded prefix increment for infinite precision integer. The integer is incremented by 1, then it returns the value.
```cpp
// ++ (prefix)
infPrecision I1("-1");
cout << "i1 = " << I1 << "\n";
cout << "++i1 = " << ++I1 << "\n";
```
```
i1 = -1
++i1 = 0
```
2. ### Postfix Increment
Overloaded postfix increment for infinite precision integer. The original integer returns firstly, then it is incremented by 1.
```cpp
++ (postfix)
infPrecision I2("100000000");
cout << "i2 = " << I2 << "\n";
cout << "i2++ = " << I2++ << ", i2 = " << I2 << "\n";
```
```
i2 = 100000000
i2++ = 100000000, i2 = 100000001
```
3. ### Prefix Decrement
Overloaded prefix decrement for infinite precision integer. The integer is decremented by 1, then it returns the value.
```cpp
-- (prefix)
infPrecision D1("0");
cout << "d1 = " << D1 << "\n";
cout << "--d1 = " << --D1 << "\n";
```
```
d1 = 0
--d1 = -1
```
4. ### Postfix Decrement
Overloaded postfix decrement for infinite precision integer. The original integer returns firstly, then it is decremented by 1.
```cpp
-- (postfix)
infPrecision D2("10000000");
cout << "d2 = " << D2 << "\n";
cout << "d2-- = " << D2-- << ", d2 = " << D2 << "\n";
```
```
d2 = 10000000
d2-- = 10000000, d2 = 9999999
```

5. ### Assignment
Overloaded binary operator = to assign one infPrecision object a new value from another infPrecision
```cpp
infPrecision A1("1234");
infPrecision A2("4321");
cout << "A1 = " << A1 << ", A2 = " << A2 << "\n";
A2 = A1;
cout << "After assign A2 = A1, A2 = " << A2 << "\n\n";
```
```
A1 = 1234, A2 = 4321
After assign A2 = A1, A2 = 1234
```

### Other member functions:
1. get_vec() <br/>
It will return the infPrecision object's stored integer vector. 
```cpp
infPrecision M1("12345");
cout << "M1 = " << M1 << "\n";
vector<int64_t> m1 = M1.get_vec();
cout << "m1 contains the vector from M1: ";
for (int64_t &i : m1)
{
    cout << i;
}
cout << "\n";
```
```
M1 = 12345
m1 contains the vector from M1: 12345
```
2. is_negative() <br/>
It will return true if the infinite-precision integer is negative, false otherwise. It is implemented by testing if the first element in the vector smaller than 0.
```cpp
infPrecision M2("-12");
cout << "Is M2 = " << M2 << " negative? " << boolalpha << (M2.is_negative()) << "\n";
infPrecision M3("0");
cout << "Is M3 = " << M3 << " negative? " << (M3.is_negative()) << "\n\n";
```
```
Is M2 = -12 negative? true
Is M3 = 0 negative? false
```
<br/>

## Overloaded operators
1. ### Addition (+)
Overloaded binary operator (+) to calculate the summation of the first infPrecision object with the second infPrecision object. It considered all situations when the two infPrecision object is negative or positive. If both integers are positive, add them digit by digit. If one digit(from right to left) is greater than 0, divide it by 10 and keep the remainder. Let the digit on its left increase by 1. If both integers are negative, they can be considered as negative of two positive integer addition (-1 + -2 = - (1 + 2)). If the first integer is positive and the second integer is negative, it is implemented as a regular subtraction. Firstly it will compare two integer's size. If they have unequal size, fill the shorter one by 0s. (if x = 1, y = -100, x will be filled to x = 001). Then again, do digit by digit subtraction. If one digit is smaller than 0, and next digit on its left will reduce by 1(If the next digit is 0, it will turn to 9 instead).
```cpp
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
```
```
add1 = 1, add2 = 9999
add1 + add2 = 10000
add3 = 15, add4 = -3015
add3 + add4 = -3000
add5 = -1111, add6 = -9999
add5 + add6 = -11110
```

2. ### Add and Assign (+=)
Overloaded binary operator (+=) to assign the new value of the left infPrecision object after it add the right infPrecision object. Directly uses (+) to implement this operator.
```cpp
int64_t add5 = -1111;
int64_t add6 = -9999;
cout << "add5 = " << add5 << ", add6 = " << add6 << "\n";
infPrecision Add5(add5);
infPrecision Add6(add6);
Add5 += Add6;
cout << "add5 += add6, add5 = " << Add5 << "\n\n";
```
```
add5 = -1111, add6 = -9999
add5 += add6, add5 = -11110
```

3. ### Subtraction (-)
Overloaded binary operator (-) to calculate the subtraction of the left infPrecision object with the right infPrecision object. It considered all situations when the two infPrecision object is negative or positive. Directly uses (+) to implement this operator.
```cpp
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
```
```
sub1 = 352331, sub2 = 4520
sub1 - sub2 = 347811
sub3 = 200000000000000000000000000000000, sub4 = -1111
sub3 - sub4 = 200000000000000000000000000001111
sub4 - sub3 = -200000000000000000000000000001111
sub5 = -10000, sub6 = -9997
sub5 - sub6 = -3
sub6 - sub5 = 3
```

4. ### Subtract and Assign (-=)
Overloaded binary operator (-=) to assign the new value of the first infPrecision object after it subtracts the second infPrecision object. Directly uses (-) to implement this operator.
```cpp
int64_t sub5 = -10000;
int64_t sub6 = -9997;
cout << "sub5 = " << sub5 << ", sub6 = " << sub6 << "\n";
infPrecision Sub5(sub5);
infPrecision Sub6(sub6);
Sub5 -= Sub6;
cout << "sub5 -= sub6, sub5 = " << Sub5 << "\n\n";
```
```
sub5 = -10000, sub6 = -9997
sub5 -= sub6, sub5 = -3
```

5. ### Multiplication (*) 
Overloaded binary operator (*) to calculate the multiplication of the first infPrecision object with the second infPrecision object. It is implemented by the idea how the multiplication is calculated by hand: each time use one of the second integer's digit(from the right to left) multiply the first integer, store the result. In the end place them in a stair shape (next summation places one index forward) and add them up.
 ```
   125
 *  15
-----
   625
+125
-----
  1875
 ```
 ```cpp
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
 ```
 ```
mul1 = 123, mul2 = 99
mul1 * mul2 = 12177
mul3 = -99, mul4 = 99
mul3 * mul4 = -9801
mul5 = -54321, mul6 = -1
mul5 * mul6 = 54321
 ```
6. ### Multiply and Assign (*=)
Overloaded binary operator (\*=) to assign the new value of the first infPrecision object after it multiplies the second infPrecision object. Directly uses (*) to implement this operator.
```cpp
int64_t mul5 = -54321;
int64_t mul6 = -1;
cout << "mul5 = " << mul5 << ", mul6 = " << mul6 << "\n";
infPrecision Mul5(mul5);
infPrecision Mul6(mul6);
Mul5 *= Mul6;
cout << "mul5 *= mul6, mul5 = " << Mul5 << "\n";
```
```
mul5 = -54321, mul6 = -1
mul5 *= mul6, mul5 = 54321
```
7. ### Negation (unary -)
Overloaded unary operator - to get the negation of infPrecision object. It is implemented by multiply -1 to the first integer in infPrecision's vector.
```cpp
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
```
```
n1 = 100
Is n1 negative? false
-n1 = -100
Is -n1 negative? true
n2 = -2, n3 = 125, n4 = -5
((-n2) - n3) * (-n4) = -615
```
8. ### equality (==)
Overloaded binary operator == to determine if the first infPrecision object is equal to the second infPrecision object. It is implemented by compare vectors in two infPrecision objects. If they are equal return true, otherwise false.
```cpp
string c1 = "100000199999000";
string c2 = "100000199999000";
string c3 = "100000019999900";
cout << "c1 = " << c1 << ", c2 = " << c2 << ", c3 = " << c3 << "\n";
infPrecision C1(c1);
infPrecision C2(c2);
infPrecision C3(c3);
cout << "c1 == c2 is " << (c1 == c2) << "\n";
cout << "c1 == c3 is " << (c1 == c3) << "\n";
```
```
c1 = 100000199999000, c2 = 100000199999000, c3 = 100000019999900
c1 == c2 is true
c1 == c3 is false
```
9. ### Inequality (!=)
Overloaded binary operator != to determine if the first infPrecision object is unequal to the second infPrecision object. It is implemented by compare vectors in two infPrecision objects. If they are unequal return true, otherwise false.
```cpp
string c1 = "100000199999000";
string c3 = "100000019999900";
cout << "c1 = " << c1 << ", c3 = " << c3 << "\n";
infPrecision C1(c1);
infPrecision C3(c3);
cout << "c1 != c3 is " << (c1 != c3) << "\n";
cout << "c3 != c3 is " << (c3 != c3) << "\n";
```
```
c1 = 100000199999000, c3 = 100000019999900
c1 != c3 is true
c3 != c3 is false
```

10. ### Smaller (<)
Overloaded binary operator < to determine if the first infPrecision object is smaller than the second infPrecision object. If both objects contains a positive integer, if the first integer has smaller size (1 < 10), return true. If they are the same size, it will compare digit by digit from right to left. If the first integer has smaller number at one digit, return true. If both objects contains a negative integer, use the reversed idea. If first integer has longer size, return true (-10 < -1). If they have the same size, it will compare digit by digit (in absolute value) from right to left. If the first integer has bigger number at one digit, return true. (-2 < -1)
```cpp
infPrecision C7("0");
infPrecision C8("1");
cout << "c7 = " << C7 << ", c8 = " << C8 << "\n";
cout << "c7 < c8 is " << (C7 < C8) << "\n";
cout << "c8 < c7 is " << (C8 < C7) << "\n\n";
```
```
c7 = 0, c8 = 1
c7 < c8 is true
c8 < c7 is false
```

11. ### Smaller than or equal to (<=)
Overloaded binary operator (<=) to determine if the first infPrecision object is smaller than or equal to the second infPrecision object. Directly uses (==) and (<) to implement this operator.
```cpp
infPrecision C7("0");
infPrecision C8("1");
infPrecision C9("1");
cout << "c7 = " << C7 << ", c8 = " << C8 << ", c9 = " << C9 << "\n";
cout << "c8 <= c9 is " << (C8 <= C9) << "\n";
cout << "c9 <= c7 is " << (C9 <= C7) << "\n\n";
```
```
c7 = 0, c8 = 1, c9 = 1
c8 <= c9 is true
c9 <= c7 is false
```

12. ### Greater (>)
Overloaded binary operator (>) to determine if the first infPrecision object is greater than the second infPrecision object. It is using exactly same idea as in overloading (<) operator.
```cpp
infPrecision C4("-999");
infPrecision C5("-1000");
cout << "c4 = " << C4 << ", c5 = " << C5 << "\n";
cout << "c4 > c5 is " << (C4 > C5) << "\n";
cout << "c5 > c4 is " << (C5 > C4) << "\n";
```
```
c4 = -999, c5 = -1000
c4 > c5 is true
c5 > c4 is false
```

13. ### Greater than or equal to (>=)
Overloaded binary operator (>=) to determine if the first infPrecision object is greater and equal to the second infPrecision object. Directly uses (==) and (>) to implement this operator.
```cpp
infPrecision C4("-999");
infPrecision C5("-1000");
infPrecision C6("-999");
cout << "c4 = " << C4 << ", c5 = " << C5 << ", c6 = " << C6 << "\n";       
cout << "c4 >= c6 is " << (C4 >= C6) << "\n";
cout << "c5 >= c6 is " << (C5 >= C6) << "\n\n";
```
```
c4 = -999, c5 = -1000, c6 = -999
c4 >= c6 is true
c5 >= c6 is false
```

14. ### Insertion (<<)
Overloaded binary operator << to print the integer stored in the infPrecision object.
```cpp
infPrecision A1("1234");
infPrecision A2("4321");
cout << "A1 = " << A1 << ", A2 = " << A2 << "\n";
```
```
A1 = 1234, A2 = 4321
```

