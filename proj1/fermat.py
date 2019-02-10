import random


def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N,k), miller_rabin(N,k)


# This function is n^3
# assume n bits for input
# we repeat this log y times, or n times.
# thus, n^2 * n is O(n^3)
# space complexity is n calls of n bits of input so O(n^2)
def mod_exp(x, y, N):
    # check n bits = n
    if y == 0:
        # constant time
        return 1
    # dividing by two is O(n) - bit shift
    z = mod_exp(x, y // 2, N)
    # mod 2 is O(1) to check the last bit
    if y % 2 == 0:
        # Exponentiation is n^2 (z * z is n^2) and mod is n^2 (because it's division) makes it O(n^2)
        return z**2 % N
    else:
        # N^2 + 2N^2 or O(N^2)
        return (x * z**2) % N


def fprobability(k):
    # subtraction is n, division is n^2, and the exponentiation is n bit shifts so the line runs in O(n^2)
    # space complexity is constant since it can never be greater than 1 or less than 0.
    return 1 - (1 / 2**k)


def mprobability(k):
    # subtraction is n and division is n^2 and the exponentiation is 2n bit shifts so the line runs in O(n^2)
    # space complexity is constant since it can never be greater than 1 or less than 0.
    return 1 - (1 / 4**k)


def fermat(N,k):
    # k has to be less than N so we know at max k is n
    # function runs k times so it is O(k*n^3)
    for item in range(0, k):
        # call is constant
        a = random.randint(2, N - 1)
        # O(n^3) for mod_exp
        if mod_exp(a, N-1, N) == 1:
            continue
        else:
            return 'composite'
    # returns are constant
    return 'prime'


def miller_rabin(N,k):
    # repeat this loop k times
    # thus the algorithm is O(n^4*k)
    for rep in range(0, k):
        p = (N - 1) * 2  # multiplication is O(n), worst case multiply so we can divide by two on the first run
        # constant time
        a = random.randint(2, N - 1)
        # checks and mod 2 check are O(1)
        # repeats this look at most log(p) times which is equivalent to O(n)
        while p > 1 and p % 2 == 0:
            # bit shift it O(1)
            p = p / 2  # reduce the exponent every
            # mod_exp is O(n^3)
            mod = mod_exp(a, p, N)
            if mod != 1:
                if mod == N - 1:
                    # still could be prime
                    break
                else:
                    return 'composite'
    return 'prime'

#179426549
#Charmichal numbers:
# 1105, 29341, 10585