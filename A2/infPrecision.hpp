/**
 * @file infPrecision.hpp
 * @author Yiqi Huang (huany132@mcmaster.ca)
 * @brief Construct a class that allow to store a infinite precision integer
 * by vector, and overload several basic binary and unary operators of
 * mathmatical arithmetics and comparison.
 * @version 0.1
 * @date 2022-11-27
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

/**
 * @brief Construct the class of infinite precision integer
 *
 */
class infPrecision
{
private:
    vector<int64_t> integer;

public:
    /**
     * @brief Construct a default infPrecision object with integer 0
     *
     */
    infPrecision()
    {
        integer = {0};
    }

    /**
     * @brief Construct a new infPrecision object by numerical string
     *
     * @param _str A numerical string
     */
    infPrecision(const string &_str)
    {
        // If an integer which is not 0, has leading 0s, throw exception
        if (_str[0] == '0' && _str.size() > 1)
        {
            throw cannot_leadZero();
        }
        else if (_str[0] == '-') // If the first character is negative sign
        {
            if (isdigit(_str[1]) == 1)
            {
                integer.push_back((_str[1] - '0') * (-1)); // store the first digit as negative number
            }
            else
            {
                throw not_numerical();
            }
            for (uint64_t i = 2; i < _str.size(); i++)
            {
                // if all char in string is numerical or the first char is '-'
                if (isdigit(_str[i]) == 1)
                {
                    integer.push_back(_str[i] - '0'); // store string char by char in integer form
                }
                else
                {
                    throw not_numerical();
                }
            }
        }
        else
        {
            for (uint64_t i = 0; i < _str.size(); i++)
            {
                // If all char in string is numerical or the first char is '-'
                if (isdigit(_str[i]) == 1)
                {
                    integer.push_back(_str[i] - '0'); // store string char by char in integer form
                }
                else
                {
                    throw not_numerical();
                }
            }
        }
    }
    /**
     * @brief Construct a new infPrecision object by 64bit fixed-width signed integer
     *
     * @param _64bit A 64bit fixed-width signed integer
     */
    infPrecision(int64_t &_64bit)
    {
        int64_t absolute = abs(_64bit); // take the absolute value of the input integer
        while (absolute > 0)
        {
            integer.push_back((int)(absolute % 10));
            absolute /= 10;
        }
        reverse(integer.begin(), integer.end());
        if (_64bit < 0)
        {
            integer[0] = integer[0] * -1;
        }
    }

    /**
     * @brief Construct a new infPrecision object by a int64_t type vector
     *
     * @param _vec A int64_t type vector
     */
    infPrecision(vector<int64_t> &_vec)
    {
        // If an integer which is not 0, has leading 0s, throw exception
        if (_vec[0] == 0 && _vec.size() > 1)
        {
            throw cannot_leadZero();
        }
        for (uint64_t i = 1; i < _vec.size(); i++)
        {
            // If an integer in the vector other than the first one is
            // a negative number, throw exception
            if (_vec[i] < 0)
            {
                throw cannot_negative();
            }

            // If an integer in the vector has two digits, throw exception
            if ((abs)(_vec[i]) >= 10)
            {
                throw cannot_mulDigits();
            }
        }
        integer = _vec;
    }

    /**
     * @brief Get the vector of a infPrecision object
     *
     */
    vector<int64_t> get_vec() const
    {
        return integer;
    }

    /**
     * @brief Test if the integer in a infPrecision object is negative
     *
     */
    bool is_negative() const
    {
        if (integer[0] < 0)
        {
            return 1; // The integer is negative
        }
        return 0;
    }

    /**
     * @brief Overloaded prefix increment for infinite precision integer
     *
     * @return infPrecision& InfPrecision object with prefix increment
     */
    infPrecision &operator++()
    {
        // Reverse the vector for doing iteration easier later
        reverse(integer.begin(), integer.end());
        uint64_t s = integer.size(); // size of the vector

        if (integer[s - 1] >= 0) // if the number is non-negative
        {
            integer[0] = integer[0] + 1; // Increment the number by 1
            for (uint64_t i = 0; i < s; i++)
            {
                // If the sum caused a digit greater than 10
                if (integer[i] >= 10)
                {
                    integer[i] = integer[i] % 10; // Each element only keep one digit
                    // If this happens on the first digit,
                    // add one more digit to hold the incremented 1.
                    if (i == s - 1)
                    {
                        integer.push_back(0);
                    }
                    integer[i + 1] += 1; // Next element will also increment by 1
                }
            }
            // reverse the vector back to normal
            reverse(integer.begin(), integer.end());
        }
        // If the number is negative, treat it as a positive number with decrement
        else if (integer[s - 1] < 0)
        {
            // Turn the number to positive
            integer[s - 1] = integer[s - 1] * (-1);
            // Decrement the last digit
            integer[0] = integer[0] - 1;

            for (uint64_t i = 0; i < s; i++)
            {
                // If one digit is negative after decrement
                if (integer[i] < 0)
                {
                    integer[i] = 9;      // It will become 9 instead
                    integer[i + 1] -= 1; // Decrement its next digit
                }
            }
            reverse(integer.begin(), integer.end());

            // Eliminate leading 0s until a non-zero integer, or a single 0 left
            if (integer.size() > 1)
            {
                while (integer[0] == 0)
                {
                    integer.erase(integer.begin());
                }
            }
            // Turn the number back to negative
            integer[0] = integer[0] * (-1);
        }

        return *this;
    }

    /**
     * @brief Overloaded postfix increment for infinite precision integer
     *
     * @return infPrecision& InfPrecision object with postfix increment
     */
    infPrecision operator++(int)
    {
        infPrecision temp = *this;
        operator++();
        return temp;
    }

    /**
     * @brief Overloaded prefix decrement for infinite precision integer
     *
     * @return infPrecision& InfPrecision object with prefix decrement
     */
    infPrecision &operator--()
    {
        // Reverse the vector for doing iteration easier later
        reverse(integer.begin(), integer.end());
        uint64_t s = integer.size(); // size of the vector

        // If the number is negative, treat it as a positive number with increment
        if (integer[s - 1] < 0)
        {
            // Turn the number to positive
            integer[s - 1] = integer[s - 1] * (-1);
            // Decrement of a negative number is plus 1
            integer[0] = integer[0] + 1;
            for (uint64_t i = 0; i < s; i++)
            {
                // If the sum caused a digit greater than 10
                if (integer[i] >= 10)
                {
                    integer[i] = integer[i] % 10; // Each element only keep one digit
                    // If this happens on the first digit,
                    // add one more digit to hold the incremented 1.
                    if (i == s - 1)
                    {
                        integer.push_back(0);
                    }
                    integer[i + 1] += 1; // Next element will also increment by 1
                }
            }
            // reverse the vector back to normal
            reverse(integer.begin(), integer.end());
            // Turn the number back to negative
            integer[0] = integer[0] * (-1);
        }
        // If the number is non-negative
        else if (integer[s - 1] >= 0)
        {
            // Decrement the last digit
            integer[0] = integer[0] - 1;

            for (uint64_t i = 0; i < s; i++)
            {
                // If there is one digit is negative after decrement and the
                // number size is greater than 1
                if (integer[i] < 0 && s > 1)
                {
                    integer[i] = 9;      // It will become 9 instead
                    integer[i + 1] -= 1; // Decrement its next digit
                }
                // If the current number is just 0
                else if (integer[i] < 0 && s == 1)
                {
                    integer[i] = -1;
                }
            }
            reverse(integer.begin(), integer.end());

            // Eliminate leading 0s until a non-zero integer, or a single 0 left
            if (integer.size() > 1)
            {
                while (integer[0] == 0)
                {
                    integer.erase(integer.begin());
                }
            }
        }
        return *this;
    }

    /**
     * @brief Overloaded postfix decrement for infinite precision integer
     *
     * @return infPrecision& InfPrecision object with postfix decrement
     */
    infPrecision operator--(int)
    {
        infPrecision temp = *this;
        operator--();
        return temp;
    }

    // Copy constructor
    infPrecision(const infPrecision &_infP)
    {
        integer = _infP.get_vec();
    }

    /**
     * @brief Overloaded binary operator = to assign one infPrecision object
     * a new value from another infPrecision
     *
     * @param _infP infPrecision object we want to assign from
     * @return infPrecision& infPrecision object with the new value
     */
    infPrecision &operator=(const infPrecision &_infP)
    {
        if (this == &_infP)
            return *this;
        integer = _infP.get_vec();
        return *this;
    }

    /**
     * @brief Exception occurs when string in constructor has
     * non-numerical character
     *
     */
    class not_numerical : public invalid_argument
    {
    public:
        not_numerical() : invalid_argument("Input is not numerical!"){};
    };

    /**
     * @brief Exception occurs when vector in constructor has an item (not the
     * first item) is negative numbers
     *
     */
    class cannot_negative : public invalid_argument
    {
    public:
        cannot_negative() : invalid_argument("Cannot have negative number other than the first item!"){};
    };

    /**
     * @brief Exception occurs when vector in constructor has an item is
     * mult-digits number.
     *
     */
    class cannot_mulDigits : public invalid_argument
    {
    public:
        cannot_mulDigits() : invalid_argument("Cannot have a number more than 1 digit in vector!"){};
    };

    /**
     * @brief Exception occurs when vector in constructor or string in
     * constructor has leading zeros in the front, where the integer is not 0
     *
     */
    class cannot_leadZero : public invalid_argument
    {
    public:
        cannot_leadZero() : invalid_argument("Cannot have leading zeros!"){};
    };
};

// Forward declaration
infPrecision operator+(const infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator+=(infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator-(const infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator-=(infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator*(const infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator*=(infPrecision &_infP1, const infPrecision &_infP2);
infPrecision operator-(const infPrecision &_infP);
bool operator==(const infPrecision &_infP1, const infPrecision &_infP2);
bool operator!=(const infPrecision &_infP1, const infPrecision &_infP2);
bool operator>(const infPrecision &_infP1, const infPrecision &_infP2);
bool operator>=(const infPrecision &_infP1, const infPrecision &_infP2);
bool operator<(const infPrecision &_infP1, const infPrecision &_infP2);
bool operator<=(const infPrecision &_infP1, const infPrecision &_infP2);

/**
 * @brief Overloaded binary operator + to calculate the summation of the
 * first infPrecision object with the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision A infPrecision object as the result of addition.
 */
infPrecision operator+(const infPrecision &_infP1, const infPrecision &_infP2)
{
    uint64_t s = 0;    // size of longer integer
    uint64_t diff = 0; // size difference between two integers
    vector<int64_t> vec1 = _infP1.get_vec();
    vector<int64_t> vec2 = _infP2.get_vec();
    reverse(vec1.begin(), vec1.end());
    reverse(vec2.begin(), vec2.end());

    // If the first integer is longer
    if (vec1.size() >= vec2.size())
    {
        s = vec1.size();
        diff = vec1.size() - vec2.size();

        // the shorter integer expands with extra 0s
        // to be th same size as the longer one
        vec2.insert(vec2.end(), diff, 0);
    }
    // If the second integer is longer
    else
    {
        s = vec2.size();
        diff = vec2.size() - vec1.size();

        // the shorter integer expands with extra 0s
        // to be th same size as the longer one
        vec1.insert(vec1.end(), diff, 0);
    }

    vector<int64_t> sum; // store result

    // if x > 0, y < 0, then x + y = x +(-y'), where y' = -y
    // This is also the process of doing regular subtraction (where both integers
    // are positive)
    if (_infP1.is_negative() == 0 && _infP2.is_negative() == 1)
    {
        int64_t temp;
        // if first integer is bigger than second integer when they are all positive
        if (_infP1 >= -(_infP2))
        {
            for (uint64_t i = 0; i < vec1.size(); i++)
            {
                // if the correct digit in vec1 is smaller
                if (vec1[i] < (abs)(vec2[i]))
                {
                    // firstly add 10 in vec1[i] then subtract vec2[i]
                    temp = (vec1[i] + 10) - (abs)(vec2[i]);
                    sum.push_back(temp);
                    for (uint64_t j = i + 1; j < vec1.size(); j++)
                    {
                        // If next digit is 0, decrease to 9
                        if (vec1[j] == 0)
                        {
                            vec1[j] = 9;
                        }
                        // Otherwise, next digit decreases by 1
                        else
                        {
                            vec1[j] -= 1;
                            break; // Only one digit decreases
                        }
                    }
                    continue;
                }
                // Otherwise, do the regular subtraction
                temp = vec1[i] - (abs)(vec2[i]);
                sum.push_back(temp);
            }
            reverse(sum.begin(), sum.end());
            // If the result is not 0, remove all leading 0s
            if (sum.size() > 1)
            {
                while (sum[0] == 0)
                {
                    sum.erase(sum.begin());
                }
            }
            infPrecision summation = infPrecision(sum);
            return summation;
        }

        // if first integer is smaller than second integer when they are all positive
        else if (_infP1 < -(_infP2))
        {
            // Recall this function in a new order so we can
            // solve this by previous if statement
            // Ex: a = 9; b = -10
            // 9 + (-10) = -(10 + (-9))
            infPrecision holder = -((-_infP2) + (-_infP1));
            return holder;
        }
    }

    // if x < 0, y > 0
    else if (_infP1.is_negative() == 1 && _infP2.is_negative() == 0)
    {
        // x + y = y + x, then we can use the first if-statement to solve this
        infPrecision holder = _infP2 + _infP1;
        return holder;
    }

    // if x < 0, y < 0
    else if (_infP1.is_negative() == 1 && _infP2.is_negative() == 1)
    {
        // x + y = -(-x + -y), where -x > 0 and -y > 0.
        // Consider this as a negation of two positive integers' addition.
        // We use the up-coming statements (when x,y>0) to implement this.
        infPrecision holder = -((-_infP1) + (-_infP2));
        return holder;
    }

    // If x > 0 and y > 0, we have the regular summation
    for (uint64_t i = 0; i < s; i++)
    {
        int64_t temp;
        temp = (vec1[i] + vec2[i]); // add digit by digit
        sum.push_back(temp);
    }
    for (uint64_t i = 0; i < s; i++)
    {
        // If the previous digits' addition is above 10, get rid of the 10
        // and add 1 to the next digit
        if (sum[i] >= 10)
        {
            sum[i] = sum[i] % 10; // Each element only keep one digit
            // If this happens on the first digit,
            // add one more digit to hold the incremented 1.
            if (i == s - 1)
            {
                sum.push_back(0);
            }
            sum[i + 1] += 1;
        }
    }

    reverse(sum.begin(), sum.end());
    infPrecision summation = infPrecision(sum);
    return summation;
}

/**
 * @brief Overloaded binary operator += to assign the new value of the left
 * infPrecision object after it plus the right infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision A infPrecision object as the result of subtraction.
 */
infPrecision operator+=(infPrecision &_infP1, const infPrecision &_infP2)
{
    _infP1 = _infP1 + _infP2;
    return _infP1;
}

/**
 * @brief Overloaded binary operator - to calculate the subtraction of the
 * left infPrecision object with the right infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision Reassign _infP1 with the result of addition
 */
infPrecision operator-(const infPrecision &_infP1, const infPrecision &_infP2)
{
    infPrecision subtraction;
    // if x > 0, y < 0, then x - y = x + y' = x + (-y), where y'= -y
    if (_infP1.is_negative() == 0 && _infP2.is_negative() == 1)
    {
        subtraction = _infP1 + (-_infP2);
        return subtraction;
    }

    // if x < 0, y < 0, then x - y = (-x') - (-y') = y'- x' = (-y) + x,
    // where x'= -x, y'= -y
    else if (_infP1.is_negative() == 1 && _infP2.is_negative() == 1)
    {
        subtraction = (-_infP2) + _infP1;
        return subtraction;
    }

    // if x < 0, y > 0, then x - y = -x'- y = -(x'+ y), where x'= -x
    else if (_infP1.is_negative() == 1 && _infP2.is_negative() == 0)
    {
        subtraction = -((-_infP1) + _infP2);
        return subtraction;
    }

    // If x > 0 and y > 0, we have the regular subtraction
    // x - y = x + (-y)
    subtraction = _infP1 + (-(_infP2));
    return subtraction;
}

/**
 * @brief Overloaded binary operator -= to assign the new value of the first
 * infPrecision object after it subtracts the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision Reassign _infP1 with the result of multiplication
 */
infPrecision operator-=(infPrecision &_infP1, const infPrecision &_infP2)
{
    _infP1 = _infP1 - _infP2;
    return _infP1;
}

/**
 * @brief Overloaded binary operator * to calculate the multiplication of the
 * first infPrecision object with the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision The result of multiplication in a new infinite
 * precision integer
 */
infPrecision operator*(const infPrecision &_infP1, const infPrecision &_infP2)
{
    vector<int64_t> vec1 = _infP1.get_vec();
    vector<int64_t> vec2 = _infP2.get_vec();

    // Return true if one of the number is negative, which means we will have
    // a negative product
    bool neg_product = 0;
    if ((vec1[0] < 0 && vec2[0] > 0) || (vec2[0] < 0 && vec1[0] > 0))
    {
        neg_product = 1;
    }
    reverse(vec1.begin(), vec1.end());
    reverse(vec2.begin(), vec2.end());

    // Create an empty vector for storing the result of multiplication
    // The digits of product <= digits of first number + digits of second number
    // Ex: 99 * 999 = 98901
    vector<int64_t> mult(vec1.size() + vec2.size(), 0);

    // Fill multiplicand by 0s if its size is smaller than mult
    if (mult.size() > vec2.size())
    {
        uint64_t diff = mult.size() - vec2.size();
        vec2.insert(vec2.end(), diff, 0);
    }

    for (uint64_t i = 0; i < vec2.size(); i++)
    {
        vector<int64_t> temp(mult.size(), 0);
        for (uint64_t j = 0; j < vec1.size(); j++)
        {
            // temp stores the current number of the multiplication between
            // a single digit from second integer and entire first integer
            temp[j] += (abs)(vec1[j]) * (abs)(vec2[i]);
            // If one element in the temp is a two-digit number
            if (temp[j] >= 10)
            {
                // add extra numbers into next digit
                temp[j + 1] += temp[j] / 10;

                // Each element only keep one digit
                temp[j] = temp[j] % 10;
            }
        }

        while (temp.back() == 0)
        {
            temp.pop_back();
            if (temp.size() == 1)
            {
                break;
            }
        }

        // The starting index of the first number is with i = 0.
        // i increases in next iteration, so we can also shift the
        // index to match up next number for the addition
        uint64_t k = i;
        for (uint64_t l = 0; l < temp.size(); l++)
        {
            mult[k] += temp[l];
            if (mult[k] >= 10)
            {
                // Each element only keep one digit
                mult[k] = mult[k] % 10;
                // Increment one in next digit
                mult[k + 1] += 1;
            }
            // Increase k for storing every digit from current number in mult
            k++;
        }

        // cout << "\n";
        temp.clear();
    }

    // If the result is not 0, remove all leading 0s
    while (mult.back() == 0)
    {
        mult.pop_back();
        if (mult.size() == 1)
        {
            break;
        }
    }
    reverse(mult.begin(), mult.end());

    // If there is one negative number, then the product should be negative
    if (neg_product == 1)
    {
        mult[0] *= (-1);
    }
    infPrecision multiplication = infPrecision(mult);
    return multiplication;
}

/**
 * @brief Overloaded binary operator *= to assign the new value of the first
 * infPrecision object after it multiplies the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return infPrecision Reassign _infP1 with the result of multiplication
 */
infPrecision operator*=(infPrecision &_infP1, const infPrecision &_infP2)
{
    _infP1 = _infP1 * _infP2;
    return _infP1;
}

/**
 * @brief Overloaded unary operator - to get the negative of
 * infPrecision object
 *
 * @param _infP The infinite precision integer
 * @return infPrecision The negated infinite precision integer
 */
infPrecision operator-(const infPrecision &_infP)
{
    vector<int64_t> temp = _infP.get_vec();
    temp[0] = temp[0] * (-1); // First digit multiplies -1
    infPrecision invert = infPrecision(temp);
    return invert;
}

/**
 * @brief Overloaded binary operator == to determine if the first
 * infPrecision object is equal to the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if two infPrecision are equal, false otherwise
 */
bool operator==(const infPrecision &_infP1, const infPrecision &_infP2)
{
    return (_infP1.get_vec() == _infP2.get_vec());
}

/**
 * @brief Overloaded binary operator != to determine if the first
 * infPrecision object is unequal to the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if two infPrecision are unequal, false otherwise
 */
bool operator!=(const infPrecision &_infP1, const infPrecision &_infP2)
{
    return (_infP1.get_vec() != _infP2.get_vec());
}

/**
 * @brief Overloaded binary operator < to determine if the first
 * infPrecision object is smaller than the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if first infPrecision is smaller than
 * second infPrecision, false otherwise
 */
bool operator<(const infPrecision &_infP1, const infPrecision &_infP2)
{
    vector<int64_t> vec1 = _infP1.get_vec();
    vector<int64_t> vec2 = _infP2.get_vec();

    // if x is negative, y is non-negative, then x<y
    if (vec1[0] < 0 && vec2[0] >= 0)
    {
        return true;
    }
    // if x,y are both negative
    else if (vec1[0] < 0 && vec2[0] < 0)
    {
        // the larger the negative integer's size, the smaller its value
        if (vec1.size() > vec2.size())
        {
            return true;
        }
        // if both negative integers have the same size
        else if (vec1.size() == vec2.size())
        {
            for (uint64_t i = 0; i < vec1.size(); i++)
            {
                // for negative numbers, if one of the digit's absolute value
                // in x is bigger than the corresponding digit's absolute value
                // in y, then x < y
                if ((abs)(vec1[i]) > (abs)(vec2[i]))
                {
                    return true;
                }
            }
        }
        return false;
    }
    // if x,y are both non-negative
    else if (vec1[0] >= 0 && vec2[0] >= 0)
    {
        // the smaller the positive integer's size, the smaller its value
        if (vec1.size() < vec2.size())
        {
            return true;
        }
        // if both positive integers have the same size
        else if (vec1.size() == vec2.size())
        {
            for (uint64_t i = 0; i < vec1.size(); i++)
            {
                // if one of the digit in x is smaller
                // than the corresponding digit in y, then x < y
                if (vec1[i] < vec2[i])
                {
                    return true;
                }
            }
        }
        return false;
    }
    return false;
}

/**
 * @brief Overloaded binary operator <= to determine if the first
 * infPrecision object is smaller than or equal to the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if first infPrecision is smaller than or equal to the
 * second infPrecision, false otherwise
 */
bool operator<=(const infPrecision &_infP1, const infPrecision &_infP2)
{
    if (_infP1 < _infP2 || _infP1 == _infP2)
    {
        return true;
    }
    return false;
}

/**
 * @brief Overloaded binary operator > to determine if the first
 * infPrecision object is greater than the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if first infPrecision is greater than the
 * second infPrecision, false otherwise
 */
bool operator>(const infPrecision &_infP1, const infPrecision &_infP2)
{
    vector<int64_t> vec1 = _infP1.get_vec();
    vector<int64_t> vec2 = _infP2.get_vec();

    // if x is non-negative, y is negative, then x>y
    if (vec1[0] >= 0 && vec2[0] < 0)
    {
        return true;
    }
    // if x,y are both negative
    else if (vec1[0] < 0 && vec2[0] < 0)
    {
        // the larger the negative integer's size, the smaller its value
        if (vec1.size() < vec2.size())
        {
            return true;
        }
        // if both negative integers have the same size
        else if (vec1.size() == vec2.size())
        {
            for (uint64_t i = 0; i < vec1.size(); i++)
            {
                // for negative numbers, if one of the digit's absolute value
                // in x is smaller than the corresponding digit's absolute value
                // in y, then x > y
                if ((abs)(vec1[i]) < (abs)(vec2[i]))
                {
                    return true;
                }
            }
        }
        return false;
    }
    // if x,y are both non-negative
    else if (vec1[0] >= 0 && vec2[0] >= 0)
    {
        // the greater the positive integer's size, the greater its value
        if (vec1.size() > vec2.size())
        {
            return true;
        }
        // if both positive integers have the same size
        else if (vec1.size() == vec2.size())
        {
            for (uint64_t i = 0; i < vec1.size(); i++)
            {
                // if one of the digit in x is greater
                // than the corresponding digit in y, then x > y
                if (vec1[i] > vec2[i])
                {
                    return true;
                }
            }
        }
        return false;
    }
    return false;
}

/**
 * @brief Overloaded binary operator >= to determine if the first
 * infPrecision object is greater than or equal to the second infPrecision object
 *
 * @param _infP1 The first infinite precision integer
 * @param _infP2 The second infinite precision integer
 * @return bool Return true if first infPrecision is greater than or equal to the
 * second infPrecision, false otherwise
 */
bool operator>=(const infPrecision &_infP1, const infPrecision &_infP2)
{
    if (_infP1 > _infP2 || _infP1 == _infP2)
    {
        return true;
    }
    return false;
}

/**
 * @brief Overloaded binary operator << to print the
 * integer stored in infPrecision object
 *
 * @param out An ostream object
 * @param _infP A infPrecision object
 * @return ostream& Returns a reference to an ostream object.
 */
ostream &operator<<(ostream &out, const infPrecision &_infP)
{
    vector<int64_t> temp = _infP.get_vec();
    for (uint64_t i = 0; i < temp.size(); i++)
        out << temp[i];
    return out;
}