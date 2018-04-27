#include "stringSearch.hpp"

#include <iostream>     // std::cout
#include <algorithm>    // std::min

UTILS_API size_t ComputeLevenshteinDistance(std::string str1, std::string str2)
{
    size_t n = str1.length();
    size_t m = str2.length();
    size_t **d = new size_t*[n + 1];
    for(size_t I = 0 ; I < n + 1 ; I ++){
        d[I] = new size_t[m + 1];
    }

    // Step 1
    if (n == 0)
    {
        return m;
    }

    if (m == 0)
    {
        return n;
    }

    // Step 2
    for (size_t i = 0; i <= n ; i++)
    {
        d[i][0] = i;
    }

    for (size_t j = 0; j <= m ; j++)
    {
        d[0][j] = j;
    }

    // Step 3
    for (size_t i = 1; i <= n; i++)
    {
        //Step 4
        for (size_t j = 1; j <= m; j++)
        {
            // Step 5
            size_t cost = (str2[j - 1] == str1[i - 1]) ? 0 : 1;

            // Step 6
            d[i][j] = std::min(
                std::min(d[i - 1][j] + 1, d[i][j - 1] + 1),
                d[i - 1][j - 1] + cost);
        }
    }
    // Returned value:
    size_t returned_value = d[n][m];
    // Free:
    for(size_t I = 0 ; I < n + 1 ; I ++){
        delete[] d[I];
    }
    delete[] d;
    // Step 7
    return returned_value;
}