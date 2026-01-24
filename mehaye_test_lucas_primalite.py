# Calcul rapide de U_k modulo n
# Suite de Lucas avec (P, Q) = (P, -P)
# U0 = 0, U1 = 1, U_{k+1} = P*U_k + P*U_{k-1}

def lucasU_Qeq_minusP(P, n, k):
    if k == 0:
        return 0
    if k == 1:
        return 1 % n

    # Matrice compagnon M = [[P, P], [1, 0]]
    a, b, c, d = P % n, P % n, 1, 0

    # Matrice résultat R = identité
    e, f, g, h = 1, 0, 0, 1

    t = k - 1
    while t > 0:
        if t & 1:
            # R = R * M
            e, f, g, h = (
                (e*a + f*c) % n,
                (e*b + f*d) % n,
                (g*a + h*c) % n,
                (g*b + h*d) % n
            )

        # M = M^2
        a2 = (a*a + b*c) % n
        b2 = (a*b + b*d) % n
        c2 = (c*a + d*c) % n
        d2 = (c*b + d*d) % n
        a, b, c, d = a2, b2, c2, d2

        t >>= 1

    # U_k = coefficient (0,0) de M^(k-1)
    return e


# Test Lucas_Mehaye : U(n-1) = 1 et U(n+1) = 0 (mod n)
def test_lucas_Mehaye(n, P):
    u_nm1 = lucasU_Qeq_minusP(P, n, n - 1) % n
    u_np1 = lucasU_Qeq_minusP(P, n, n + 1) % n
    ok = (u_nm1 == 1) and (u_np1 == 0)
    return ok, u_nm1, u_np1


# Miller–Rabin déterministe (64 bits)
# Utilisé uniquement pour validation expérimentale
def miller_rabin(n):
    if n < 2:
        return False

    petits = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    for p in petits:
        if n == p:
            return True
        if n % p == 0:
            return False

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
n = int(input("Entrez un entier n a tester : ").strip())

print("\nAnalyse de n =", n)
print("Classe modulo 10 :", n % 10)

if n % 10 in (3, 7):
    suites = [5]
elif n % 10 in (1, 9):
    suites = [10, 15, 20, 25, 30]
else:
    suites = []

if not suites:
    print("\nAucune suite de Lucas pour cette classe.")
else:
    for P in suites:
        print(f"\nSuite de Lucas (P,Q)=({P},{-P})")
        ok, u_nm1, u_np1 = test_lucas_Mehaye(n, P)
        print("  U(n-1) mod n :", u_nm1)
        print("  U(n+1) mod n :", u_np1)
        print("  Lucas_Mehaye :", "OUI" if ok else "NON")
        if ok:
            break

print("Miller-Rabin (validation) :", "OUI" if miller_rabin(n) else "NON")
