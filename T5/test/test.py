from random import randint

N1 = randint(3, 10)
N2 = randint(3, 10)
K = randint(3, 5)

def summ(a, b):
    res = 0
    for i in range(K):
        res += ((a // (N1 ** i)) % N1 + (b // (N2 ** i)) % N2) * (N1 + N2)**i

    return res

for j in range(int(input())):
    a = [randint(0, N1-1) for i in range(K)]
    b = [randint(0, N2-1) for i in range(K)]
    a_dec = sum(a[i]*(N1**i) for i in range(K))
    b_dec = sum(b[i]*(N2**i) for i in range(K))
    res = [a[i] + b[i] for i in range(K)]
    res_dec = sum([res[i] * (N1 + N2)**i for i in range(K)])
    print(f"Test number {j}:")
    print(f"A = {a} |-> {a_dec}")
    print(f"B = {b} |-> {b_dec}")
    print(f"res = {res} |-> {res_dec}")
    print(f"f(a, b) = {summ(a_dec, b_dec)}")
    print(res_dec == summ(a_dec, b_dec))
    print()
    