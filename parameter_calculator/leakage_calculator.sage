from sage.all import sum


def pr_rank(N, s, rank):
    H_N = sum(1 / (k**s) for k in range(1, N + 1))
    return (1 / (H_N)) * (1 / (rank**s))


def expected_rank_count(N, s, rank, n):
    return n * pr_rank(N, s, rank)


def run(n, N, t, ell):
    leakage_count = 0

    num_reports_leaked = 0
    for i in range(1, 300):
        erc = expected_rank_count(N, 1.03, i, n)

        if erc < t and erc > ell:
            leakage_count += 1
            num_reports_leaked += erc

    print("Total Values Leaked (expected): ", leakage_count)
    print(f"Fraction of reports leaked (expected): {num_reports_leaked / n:0.2f}")
