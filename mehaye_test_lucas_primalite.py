# Matrices 2x2 modulo n

def mat_mul(A, B, n):
    return [
        [(A[0][0]*B[0][0] + A[0][1]*B[1][0]) % n,
         (A[0][0]*B[0][1] + A[0][1]*B[1][1]) % n],
        [(A[1][0]*B[0][0] + A[1][1]*B[1][0]) % n,
         (A[1][0]*B[0][1] + A[1][1]*B[1][1]) % n]
    ]


def mat_pow(M, k, n):
    R = [[1, 0], [0, 1]]
    while k > 0:
        if k & 1:
            R = mat_mul(R, M, n)
        M = mat_mul(M, M, n)
        k >>= 1
    return R

# Suite de Lucas U_k modulo n


def lucas_U(P, Q, n, k):
    if k == 0:
        return 0
    if k == 1:
        return 1

    M = [[P % n, (-Q) % n],
         [1,     0       ]]

    Mk = mat_pow(M, k - 1, n)
    return Mk[0][0]


def test_lucas(n, P):
    u_nm1 = lucas_U(P, -P, n, n - 1) % n
    u_np1 = lucas_U(P, -P, n, n + 1) % n
    ok = (u_nm1 == 1 % n) and (u_np1 == 0)
    return ok, u_nm1, u_np1


# Test de Miller–Rabin

def miller_rabin(n):
    if n < 2:
        return False

    petits_premiers = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    for p in petits_premiers:
        if n == p:
            return True
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    bases = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False

    return True


# Programme principal

n = int(input("Entrez un entier n à tester : ").strip())

print("\nAnalyse de n =", n)
print("Classe modulo 10 :", n % 10)

if n % 10 in (3, 7):
    suites = [5]
elif n % 10 in (1, 9):
    suites = [10, 15, 20, 25, 30]
else:
    suites = []

if not suites:
    print("\nPas de suite de Lucas prévue pour cette classe modulo 10.")
else:
    for P in suites:
        print(f"\nTest de la suite de Lucas (P,Q) = ({P}, {-P})")

        ok, u_nm1, u_np1 = test_lucas(n, P)

        print("  Valeur de U_(n-1) mod n :", u_nm1)
        print("  Valeur de U_(n+1) mod n :", u_np1)
        print("  Congruences de Lucas satisfaites :", "OUI" if ok else "NON")

        if ok:
            break

mr = miller_rabin(n)
print(" Miller–Rabin :", "OUI" if mr else "NON")
