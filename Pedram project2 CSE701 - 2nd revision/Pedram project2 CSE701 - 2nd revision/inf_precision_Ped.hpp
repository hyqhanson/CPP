#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

vector<uint64_t> addingTwoPos(vector<uint64_t> pp1, vector<uint64_t> pp2);

vector<uint64_t> subtracting(vector<uint64_t> pp1, vector<uint64_t> pp2);

int64_t tasavi(const vector<uint64_t> pp1, const vector<uint64_t> pp2);

class inf_precision_Ped
{

public:
    inf_precision_Ped(const string &_str_Ped)
    {
        // making sure the first char is '-' or '+' or a digit
        if ((isdigit(_str_Ped[0]) == 1)) //(1 first digit when it is positive)
        {
            int_vec.push_back(1);
            if (_str_Ped[0] != '0')
            {
                int_vec.push_back((_str_Ped[0] - '0'));
            }
            if (_str_Ped[0] == '0' && _str_Ped.size() == 1) // if it is only one zero
            {
                int_vec.push_back((_str_Ped[0] - '0'));
            }
        }
        else if (_str_Ped[0] == '-' && _str_Ped.size() == 1) // it is only "-"
        {
            throw not_integer;
        }
        else if (_str_Ped[0] == '+' && _str_Ped.size() == 1) // it is only "+"
        {
            throw not_integer;
        }
        else if (_str_Ped[0] == '-') //(put 0 first digit if it is negative)
        {
            int_vec.push_back(0);
        }
        else if (_str_Ped[0] == '+') //(put 1 first digit if it is positive)
        {
            int_vec.push_back(1);
        }
        else
        {
            throw not_integer;
            // cout << "Error: the value " << _str_Ped << " is not an integer.\n";
            // exit(0);
        }

        uint64_t pppp = 1;
        for (uint64_t i = 1; i < _str_Ped.size(); i++)
        {
            if (_str_Ped[i] == '0' && pppp == 1 && i != _str_Ped.size() - 1)
            {
                continue;
            }

            if (isdigit(_str_Ped[i]) == 1)
            {
                pppp = 0;
                // cout << _str_Ped[i];
                int_vec.push_back(_str_Ped[i] - '0');
            }
            else
            {
                throw not_integer;
                // cout << "Error: the value " << _str_Ped << " is not an integer.\n";
                // exit(0);
            }
        }
    }

    inf_precision_Ped(const vector<uint64_t> &_pedram)
    {
        int_vec = _pedram;
    }

    inf_precision_Ped()
    {
        throw not_integer;
    }

    vector<uint64_t> get_vec() const
    {
        return int_vec;
    }

    vector<uint64_t> changeV(const vector<uint64_t> &_pp)
    {
        int_vec.clear();
        int_vec = _pp;
        /*for (uint64_t i = 0; i < _pp.size(); i++)
        {
            int_vec.push_back(_pp[i]);
        }*/
        return int_vec;
    }

    // Reference : https://en.cppreference.com/w/cpp/language/operators
    inf_precision_Ped &operator++()
    {
        vector<uint64_t> inc_ped = {1, 1}; // y, first 1 is for +
        if (int_vec[0] == 1)               // int_vec is a positive number is like x+y
        {
            if (int_vec.size() >= inc_ped.size())
            {
                int_vec = addingTwoPos(int_vec, inc_ped);
                int_vec.insert(int_vec.begin(), 1); // it always will be positive
            }
            else
            {
                int_vec = addingTwoPos(inc_ped, int_vec);
                int_vec.insert(int_vec.begin(), 1); // it always will be positive
            }
        }
        else // like -x+y => -(x-y)
        {
            int_vec.erase(int_vec.begin());
            int_vec.insert(int_vec.begin(), 1); // absolute

            int64_t ind = tasavi(int_vec, inc_ped);
            if (+0 < ind) // abs(pp1)> abs(pp2)
            {
                // return pp;
                int_vec = subtracting(int_vec, inc_ped);
                int_vec.insert(int_vec.begin(), 0); // it always will be negative
            }
            else if (+0 > ind) // abs(pp1)< abs(pp2)
            {
                int_vec = subtracting(inc_ped, int_vec);
                int_vec.insert(int_vec.begin(), 1); // it always will be positive
            }
            else // abs(pp1) = abs(pp2)
            {
                vector<uint64_t> pp = {0};
                int_vec = pp;
                int_vec.insert(int_vec.begin(), 1); // +0
            }
        }

        return *this; // return new value by reference
    }

    inf_precision_Ped operator++(int)
    {
        inf_precision_Ped old = *this; // copy old value
        operator++();                  // prefix increment
        return *this;                  // return old value
    }

    // prefix decrement
    inf_precision_Ped &operator--()
    {
        vector<uint64_t> inc_ped = {0, 1}; // y, first 0 is for -
        if (int_vec[0] == 1)               // int_vec is a positive number is like x-y
        {
            inc_ped.erase(inc_ped.begin());
            inc_ped.insert(inc_ped.begin(), 1); // absolute

            int64_t ind = tasavi(int_vec, inc_ped);
            if (+0 < ind) // abs(pp1)> abs(pp2)
            {
                // return pp;
                int_vec = subtracting(int_vec, inc_ped);
                int_vec.insert(int_vec.begin(), 1); // it always will be positive
            }
            else if (+0 > ind) // abs(pp1)< abs(pp2)
            {
                int_vec = subtracting(inc_ped, int_vec);
                int_vec.insert(int_vec.begin(), 0); // it always will be negative
            }
            else // abs(pp1) = abs(pp2)
            {
                vector<uint64_t> pp = {0};
                int_vec = pp;
                int_vec.insert(int_vec.begin(), 1); // +0
            }
        }
        else // like -(x+y)
        {
            inc_ped.erase(inc_ped.begin());
            inc_ped.insert(inc_ped.begin(), 1); // absolute
            int_vec.erase(int_vec.begin());
            int_vec.insert(int_vec.begin(), 1); // absolute

            if (int_vec.size() >= inc_ped.size())
            {
                int_vec = addingTwoPos(int_vec, inc_ped);
                int_vec.insert(int_vec.begin(), 0); // it always will be negative
            }
            else
            {
                int_vec = addingTwoPos(inc_ped, int_vec);
                int_vec.insert(int_vec.begin(), 0); // it always will be negative
            }
        }

        return *this; // return new value by reference
    }

    // postfix decrement
    inf_precision_Ped operator--(int)
    {
        inf_precision_Ped old = *this; // copy old value
        operator--();                  // prefix decrement
        return *this;                  // return old value
    }

    inline static invalid_argument not_integer = invalid_argument("The entered value is not an integer!");

private:
    vector<uint64_t> int_vec;
    /* data */
};

vector<uint64_t> addingTwoPos(vector<uint64_t> pp1, vector<uint64_t> pp2)
{
    vector<uint64_t> pp;
    pp1.erase(pp1.begin());
    pp2.erase(pp2.begin());
    uint64_t bb = 0;
    int64_t aa;
    for (uint64_t i = 0; i < pp2.size(); i++)
    {
        aa = pp1[(pp1.size() - 1 - i)] + pp2[(pp2.size() - 1 - i)];

        if (aa + bb > 10)
        {
            pp.push_back(aa % 10 + bb);
            bb = 1;
        }
        else if (aa + bb == 10)
        {
            pp.push_back(0);
            bb = 1;
        }
        else // aa+bb<10
        {
            pp.push_back(aa + bb);
            bb = 0;
        }
    }

    // taking care of the untouched digits in the upper number
    for (uint64_t i = (pp1.size() - pp2.size()); i > 0; i--)
    {
        if ((pp1[i - 1] + bb) > 9)
        {
            pp.push_back(0);
            bb = 1;
            continue;
        }
        pp.push_back(pp1[i - 1] + bb);
        bb = 0;
    }
    if (bb != 0)
    {
        pp.push_back(1); // pushing the last 1
    }

    std::reverse(pp.begin(), pp.end());
    return pp;
}

vector<uint64_t> subtracting(vector<uint64_t> pp1, vector<uint64_t> pp2)
{
    vector<uint64_t> pp;
    pp1.erase(pp1.begin());
    pp2.erase(pp2.begin());
    int64_t bb = 0;
    int64_t aa;
    for (uint64_t i = 0; i < pp2.size(); i++)
    {
        aa = pp1[(pp1.size() - 1 - i)] - pp2[(pp2.size() - 1 - i)];
        if (aa + bb < 0)
        {
            pp.push_back(aa + bb + 10);
            bb = -1;
        }
        else if (aa + bb == 0)
        {
            pp.push_back(0);
            bb = 0;
        }
        else // aa+bb>0
        {
            pp.push_back(aa + bb);
            bb = 0;
        }
    }

    for (uint64_t i = (pp1.size() - pp2.size()); i > 0; i--)
    {
        if (+0 > (pp1[i - 1] + bb))
        {
            pp.push_back(10 + bb); // only when bb=-1 and pp1[i - 1]=0
            bb = -1;
            continue;
        }
        pp.push_back(pp1[i - 1] + bb);
        bb = 0;
    }
    std::reverse(pp.begin(), pp.end());

    vector<uint64_t> ppfinal;
    uint64_t stopCr = 0;
    for (uint64_t i = 0; i < pp.size(); i++) // if the beginning has zeros
    {

        if (pp[i] == 0 && stopCr == 0)
        {
            continue;
        }
        else // (pp[i] != 0)
        {
            ppfinal.push_back(pp[i]);
            stopCr = 1;
        }
    }
    return ppfinal;
}

int64_t tasavi(const vector<uint64_t> pp1, const vector<uint64_t> pp2)
{
    // eqaullity = 0, smaller = -1, larger = 1,
    int64_t ind = 0;

    if (pp1[0] == 1 && pp2[0] == 1)
    {
        if (pp1.size() > pp2.size())
        {
            ind = 1; // pp1 larger
        }
        else if (pp1.size() < pp2.size())
        {
            ind = -1; // pp1 smaller
        }
        else
        {
            for (uint64_t i = 1; i < pp1.size(); i++)
            {
                ind = 0; // Equal
                if (pp1[i] > pp2[i])
                {
                    ind = 1; // pp1 larger
                    break;
                }
                else if (pp1[i] < pp2[i])
                {
                    ind = -1; // pp1 smaller
                    break;
                }
                else
                {
                    continue;
                }
            }
        }
    }

    else if (pp1[0] == 0 && pp2[0] == 0)
    {
        if (pp1.size() < pp2.size())
        {
            ind = 1; // p1 larger
        }
        else if (pp1.size() > pp2.size())
        {
            ind = -1;
        }
        else
        {
            for (uint64_t i = 1; i < pp1.size(); i++)
            {
                ind = 0; // Equal
                if (pp1[i] < pp2[i])
                {
                    ind = 1; // pp1 larger
                    break;
                }
                else if (pp1[i] > pp2[i])
                {
                    ind = -1; // pp1 smaller
                    break;
                }
                else
                {
                    continue;
                }
            }
        }
    }

    else if (pp1[0] == 1 && pp2[0] == 0)
    {
        ind = 1; // pp1 larger
    }

    else if (pp1[0] == 0 && pp2[0] == 1)
    {
        ind = -1; // pp1 smaller
    }

    return ind;
}

vector<uint64_t> multiplying(vector<uint64_t> pp1, vector<uint64_t> pp2)
{
    vector<uint64_t> pp;
    vector<uint64_t> sum_pp = {0};
    pp1.erase(pp1.begin());
    pp2.erase(pp2.begin());
    uint64_t bb = 0;
    int64_t aa;

    for (uint64_t i = 0; i < pp2.size(); i++)
    {
        for (uint64_t j = 0; j < pp1.size(); j++)
        {
            aa = pp1[(pp1.size() - 1 - j)] * pp2[(pp2.size() - 1 - i)];

            if (aa + bb > 10)
            {
                pp.push_back((aa + bb) % 10);
                bb = (aa + bb) / 10;
            }
            else if (aa + bb == 10)
            {
                pp.push_back(0);
                bb = 1;
            }
            else // aa+bb<10
            {
                pp.push_back(aa + bb);
                bb = 0;
            }
        }
        if (bb != 0)
        {
            pp.push_back(bb); // pushing the last 1
            bb = 0;           // default for next iteration for i
        }

        std::reverse(pp.begin(), pp.end());

        // put the sign b since adding revmoves the first digit(so any number is ok) then send it for +
        pp.insert(pp.begin(), 1);
        sum_pp.insert(sum_pp.begin(), 1);

        if (pp.size() >= sum_pp.size())
        {
            sum_pp = addingTwoPos(pp, sum_pp); // give us back without first digit which is sign but always pos
        }
        else
        {
            sum_pp = addingTwoPos(sum_pp, pp);
        }

        pp.clear();
        for (uint64_t k = 0; k < i + 1; k++)
        {
            pp.push_back(0); // the next rows of multiplication
        }
    }
    return sum_pp;
}

inf_precision_Ped operator+(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator+=(inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator-(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator-=(inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator*(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator*=(inf_precision_Ped &_v, const inf_precision_Ped &_w);
inf_precision_Ped operator-(const inf_precision_Ped &_v);
bool operator==(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
bool operator<=(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
bool operator>=(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
bool operator!=(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
bool operator<(const inf_precision_Ped &_v, const inf_precision_Ped &_w);
bool operator>(const inf_precision_Ped &_v, const inf_precision_Ped &_w);

inf_precision_Ped operator+(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    vector<uint64_t> gg;
    // if both ++ is like x+y
    if (pp1[0] == 1 && pp2[0] == 1)
    {
        if (pp1.size() >= pp2.size())
        {
            gg = addingTwoPos(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else
        {
            gg = addingTwoPos(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
    }
    // both -- is like -(x+y)
    else if (pp1[0] == 0 && pp2[0] == 0)
    {
        if (pp1.size() >= pp2.size())
        {
            gg = addingTwoPos(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else
        {
            gg = addingTwoPos(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
    }
    // +- is like x-y
    else if (pp1[0] == 1 && pp2[0] == 0)
    {
        pp2.erase(pp2.begin());
        pp2.insert(pp2.begin(), 1); // absolute

        int64_t ind = tasavi(pp1, pp2);
        if (+0 < ind) // abs(pp1)> abs(pp2)
        {
            // return pp;
            gg = subtracting(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else if (+0 > ind) // abs(pp1)< abs(pp2)
        {
            gg = subtracting(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else // abs(pp1) = abs(pp2)
        {
            vector<uint64_t> pp = {0};
            gg = pp;
            gg.insert(gg.begin(), 1); // +0
        }
    }
    // -+ is like -(x-y)
    else // if (pp1[0] == 0 && pp2[0] == 1)
    {
        pp1.erase(pp1.begin());
        pp1.insert(pp1.begin(), 1); // absolute

        int64_t ind = tasavi(pp1, pp2);
        if (+0 < ind) // abs(pp1)> abs(pp2)
        {
            // return pp;
            gg = subtracting(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else if (+0 > ind) // abs(pp1)< abs(pp2)
        {
            gg = subtracting(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else // abs(pp1) = abs(pp2)
        {
            vector<uint64_t> pp = {0};
            gg = pp;
            gg.insert(gg.begin(), 1); // +0
        }
    }
    inf_precision_Ped temp(gg);
    return temp;
}

inf_precision_Ped operator-(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    vector<uint64_t> gg;
    // if both ++ is like x-y
    if (pp1[0] == 1 && pp2[0] == 1)
    {
        int64_t ind = tasavi(pp1, pp2);
        if (+0 < ind) // abs(pp1)> abs(pp2)
        {
            gg = subtracting(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else if (+0 > ind) // abs(pp1)< abs(pp2)
        {
            gg = subtracting(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else // abs(pp1) = abs(pp2)
        {
            vector<uint64_t> pp = {0};
            gg = pp;
            gg.insert(gg.begin(), 1); // +0
        }
    }
    // both -- is like -(x-y)
    else if (pp1[0] == 0 && pp2[0] == 0)
    {
        pp1.erase(pp1.begin());
        pp1.insert(pp1.begin(), 1); // absolute
        pp2.erase(pp2.begin());
        pp2.insert(pp2.begin(), 1); // absolute

        int64_t ind = tasavi(pp1, pp2);
        if (+0 < ind) // abs(pp1)> abs(pp2)
        {
            // return pp;
            gg = subtracting(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else if (+0 > ind) // abs(pp1)< abs(pp2)
        {
            gg = subtracting(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else // abs(pp1) = abs(pp2)
        {
            vector<uint64_t> pp = {0};
            gg = pp;
            gg.insert(gg.begin(), 1); // +0
        }
    }
    // +- is like x+y
    else if (pp1[0] == 1 && pp2[0] == 0)
    {
        pp2.erase(pp2.begin());
        pp2.insert(pp2.begin(), 1); // changing it to positive
        if (pp1.size() >= pp2.size())
        {
            gg = addingTwoPos(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else
        {
            gg = addingTwoPos(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
    }
    // -+ is like -(x+y)
    else // if (pp1[0] == 0 && pp2[0] == 1)
    {
        pp1.erase(pp1.begin());
        pp1.insert(pp1.begin(), 1); // changing it to positive
        if (pp1.size() >= pp2.size())
        {
            gg = addingTwoPos(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
        else
        {
            gg = addingTwoPos(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be negative
        }
    }

    inf_precision_Ped temp(gg);
    return temp;
}

bool operator>(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    if (ind == 1)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

bool operator<(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    int64_t a = -1;
    if (ind == a)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

bool operator==(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    if (ind == 0)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

bool operator<=(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    if (ind == 0 or ind < 0)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

bool operator>=(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    if (ind == 0 or ind > 0)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

bool operator!=(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    bool AA = true;
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    int64_t ind = tasavi(pp1, pp2);
    if (ind < 0 or ind > 0)
    {
        // cout << "true\n";
        AA = true;
    }
    else
    {
        // cout << "false\n";
        AA = false;
    }
    return AA;
}

inf_precision_Ped operator+=(inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp = (_v + _w).get_vec();
    pp = _v.changeV(pp);
    inf_precision_Ped temp(pp);
    return temp;
}

inf_precision_Ped operator-=(inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp = (_v - _w).get_vec();
    pp = _v.changeV(pp);
    inf_precision_Ped temp(pp);
    return temp;
}

inf_precision_Ped operator*=(inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp = (_v * _w).get_vec();
    pp = _v.changeV(pp);
    inf_precision_Ped temp(pp);
    return temp;
}

inf_precision_Ped operator*(const inf_precision_Ped &_v, const inf_precision_Ped &_w)
{
    vector<uint64_t> pp1 = _v.get_vec();
    vector<uint64_t> pp2 = _w.get_vec();
    vector<uint64_t> gg;
    // if both ++ is like x*y
    if (pp1[0] == 1 && pp2[0] == 1)
    {
        if (pp1.size() >= pp2.size())
        {

            gg = multiplying(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else
        {
            gg = multiplying(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
    }
    // both -- is like (x*y)
    else if (pp1[0] == 0 && pp2[0] == 0)
    {
        if (pp1.size() >= pp2.size())
        {
            gg = multiplying(pp1, pp2);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
        else
        {
            gg = multiplying(pp2, pp1);
            gg.insert(gg.begin(), 1); // it always will be positive
        }
    }
    // +- is like -(x*y)
    else if (pp1[0] == 1 && pp2[0] == 0)
    {
        if (pp1.size() >= pp2.size())
        {
            gg = multiplying(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be positive
        }
        else
        {
            gg = multiplying(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be positive
        }
    }
    // -+ is like -(x*y)
    else // if (pp1[0] == 0 && pp2[0] == 1)
    {
        if (pp1.size() >= pp2.size())
        {
            gg = multiplying(pp1, pp2);
            gg.insert(gg.begin(), 0); // it always will be positive
        }
        else
        {
            gg = multiplying(pp2, pp1);
            gg.insert(gg.begin(), 0); // it always will be positive
        }
    }
    inf_precision_Ped temp(gg);
    return temp;
}

inf_precision_Ped operator-(const inf_precision_Ped &_v)
{
    vector<uint64_t> pp = _v.get_vec();
    if (pp[0] == 0) // it is negative number
    {
        pp[0] = 1;
    }
    else
    {
        pp[0] = 0;
    }
    inf_precision_Ped temp(pp);
    return temp;
}

ostream &operator<<(ostream &out, const inf_precision_Ped &v)
{
    vector<uint64_t> pp = v.get_vec();
    if (pp.empty())
    {
        out << "There is no number.\n";
        return out;
    }
    else
    {
        if (pp[0] == 0 && pp[1] != 0) // if the number is negative the first digit is 0 so print "-", otherwise skip
        {
            out << '-';
            for (uint64_t i = 1; i < pp.size(); i++)
                out << pp[i];
        }
        else if (pp[1] == 0) // if the second digit is zero it is just a zero
        {
            out << pp[1];
        }
        else
        {
            for (uint64_t i = 1; i < pp.size(); i++)
                out << pp[i];
        }
        return out;
    }
}