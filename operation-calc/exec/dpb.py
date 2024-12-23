from condition import Condition
from mytype import Op
import fp, fp2
import montgomery as mont
import sympy
import mycgl
from parameter import to_str

"""
1. GenerateBasis
"""


# Calcualte x(P),x(Q) from montgomery coefficient A
# where ([f]P,[f]Q) is a basis of E[2^n] when p = 2^n f - 1
def montgomery_basis_x(c: Condition, bit_length_of_p: Op) -> Op:
    step0 = (
        2 * fp.square()
        + fp.add()
        + fp.exp(
            c,
            bit_length_of_n=bit_length_of_p - 2,
            stand_bit_of_n=(bit_length_of_p - 2) / 2,
        )
    )

    prob = 1 / 2
    step1 = (
        fp2.mul(c)
        + 2 * fp2.mul(c)
        + fp2.square(c)
        + fp.square()
        + fp.add()
        + fp.exp(
            c,
            bit_length_of_n=bit_length_of_p - 2,
            stand_bit_of_n=(bit_length_of_p - 2) / 2,
        )
        + fp.square()
    )
    return step0 + 1 / prob * step1


# Calculate x(P),x(Q) from montgomery coefficient A
# where (P,Q) is a basis of E[2^n] when p = 2^n f - 1
def montgomery_basis(c: Condition, bit_length_of_p: Op, n: Op) -> Op:
    return montgomery_basis_x(
        c, bit_length_of_p=bit_length_of_p
    ) + 2 * mont.point_ladder(c, bit_length_of_m=bit_length_of_p - n)


# Calcualte phi(Q) and generate P in E[2^n] such that (P, phi(Q)) is a basis.
def paired_basis(c: Condition, bit_length_of_p: Op, n: Op) -> Op:
    return (
        # Calcualte phi(Q)
        (n / 2) * mont.four_iso_eval(c)
        # Calcuate paired basis
        + mycgl.generate_basis(c, bit_length_of_p, n)
    )


# Calculate x(P-Q), x(P+Q) from x(P),x(Q)
def sub(c: Condition, bit_length_of_p: Op) -> Op:
    # Calculate on projective coordinate
    return (
        # t0 - t5
        4 * fp2.mul(c)
        + 3 * fp2.square(c)
        + 3 * fp2.add()
        # y0, y1
        + 6 * fp2.mul(c)
        + 2 * fp2.square(c)
        + 4 * fp2.add()
        # S, U, T
        + 6 * fp2.mul(c)
        + 3 * fp2.add()
        + fp2.sqrt(c, bit_length_of_p)
        # X, Z
        + 2 * fp2.add()
    )
    # Calculate on normal coordinate
    # return (
    #     2 * fp2.inv(c)
    #     # (Xp - Xq)^2
    #     + fp2.add()
    #     + fp2.square(c)
    #     # (Xp + Xq + A)(Xp - Xq)^2
    #     + 2 * fp2.add()
    #     + fp2.mul(c)
    #     # Yp = Xp^3 + AXp^2 + Xp, Yq = Xp^3 + AXp^2 + Xp
    #     + 2 * (fp2.square(c) + 2 * fp2.mul(c) + 2 * fp2.add())
    #     # Yp - √(YpYq) - √(YpYq) + Yq - (Xp + Xq + A)(Xp - Xq)^2
    #     + fp2.sqrt(c, bit_length_of_p)
    #     + 4 * fp2.add()
    # )


"""
2. DecideBacktrack
"""


# Use 2-isogeny in last step of 2^n-isogeny calculation, then
# use the j-invariant to check the backtrack point.
# Input must be x(P), x(Q), x(P-Q)
def backtrack_by_two_isog(c: Condition, bit_length_of_p: Op, n: Op) -> Op:
    prob_backtrack = 1 / 3

    # Cost to calculate 2^{n-1}P and 2-isogeny and j-invariant
    check_backtrack = (
        mont.xDBLe_proj(c, n - 1)
        + mont.two_iso_curve(c, bit_length_of_p)
        + mont.j_inv(c)
    )

    return (
        # Remove one 4-isogeny calculation
        -mont.two_e_iso(c, n, 0)
        + mont.two_e_iso(c, n - 2, 0)
        # Add 2times 2-isogeny calculation
        + 2 * mont.two_iso_curve(c, bit_length_of_p)
        + mont.two_iso_eval(c)
        # Backtrack j-invariant calculation
        + mont.j_inv(c)
        # If p is backtrack
        + prob_backtrack * check_backtrack
        # If q is backtrack
        + prob_backtrack * (2 * check_backtrack)
        # IF p - q is backtrack
        + prob_backtrack * (2 * check_backtrack)
    )


# Use the last j-invariant to check the backtrack point.
# Input must be x(P), x(Q), x(P-Q)
# P' = 2^{n-2}P, Q' = 2^{n-2}Q, R' = 2^{n-2}(P-Q)
# P' <- P', P' + 2Q' = (P' + Q') + Q'
# Q' <- Q', Q' + 2R' = P' + (P' - Q')
# R' <- R', R' + 2Q' = P' + Q'
def backtrack_by_four_isog(c: Condition, n: Op) -> Op:
    prob_backtrack = 1 / 3

    # Cost to calculate Q',Q'+2R' and 4-isogeny and j-invariant
    check_backtrack = 1 / 2 * (mont.four_iso_curve(c) + mont.j_inv(c)) + 1 / 2 * (
        mont.xDBLADD_proj(c) + mont.four_iso_curve(c) + mont.j_inv(c)
    )

    return (
        # Calculate last j-invarinant in the previous step
        mont.j_inv(c)
        # Calculate P',Q',R'
        + 3 * mont.xDBLe_proj(c, n - 2)
        # If q is backtrack
        + prob_backtrack * check_backtrack
        # If r is backtrack
        + prob_backtrack * (2 * check_backtrack)
        # If p is backtrack
        + prob_backtrack * (2 * check_backtrack)
    )


# Use phi(Q) to calcualte backtrack point.
# To use this method, you have to calcualte phi(Q)
def backtrack_by_phi_q(c: Condition, n: Op, phi_q_ready: bool = False) -> Op:
    prob_backtrack = 1 / 3

    # Cost to calculate 2^{n-1}P
    check_backtrack = mont.xDBLe_proj(c, n - 1)

    return (
        # Calculate phi(Q)
        ((n / 2) * mont.four_iso_eval(c) if phi_q_ready else 0)
        # Calculate 2^{n-1}phi(Q)
        + mont.xDBLe_proj(c, n - 1)
        # If p is backtrack
        + prob_backtrack * check_backtrack
        # If q is backtrack
        + prob_backtrack * (2 * check_backtrack)
        # If p-q is backtrack
        + prob_backtrack * (2 * check_backtrack)
    )


def main():
    c = Condition(is_p_equal_3_mod_4=True)
    p = sympy.var("p")
    bit_length_of_p = 256  # sympy.log(p)
    n = bit_length_of_p - 10

    # GenerateBasis
    gb_method0 = montgomery_basis(c, bit_length_of_p, n=n) + sub(c, bit_length_of_p)
    gb_method1 = paired_basis(c, bit_length_of_p, n) + sub(c, bit_length_of_p)
    print(f"1. GenerateBasis")
    print(
        f"    MontgomeryBasis: {to_str(gb_method0)} = {to_str(fp.reduce(gb_method0))}"
    )
    print(f"    PairedBasis: {to_str(gb_method1)} = {to_str(fp.reduce(gb_method1))}")

    # DecideBacktrack
    db_method0 = backtrack_by_two_isog(c, bit_length_of_p, n)
    db_method1 = backtrack_by_four_isog(c, n)
    db_method2 = backtrack_by_phi_q(c, n, phi_q_ready=False)
    print(f"2. DecideBacktrack")
    print(f"    ByTwoIsog: {to_str(db_method0)} = {to_str(fp.reduce(db_method0))}")
    print(f"    ByFourIsog: {to_str(db_method1)} = {to_str(fp.reduce(db_method1))}")
    print(f"    ByPhiQ: {to_str(db_method2)} = {to_str(fp.reduce(db_method2))}")

    # Ladder
    ld_method0 = mont.three_point_ladder(c, bit_length_of_m=n)
    print(f"3. Ladder")
    print(
        f"    MontgomeryLadder: {to_str(ld_method0)} = {to_str(fp.reduce(ld_method0))}"
    )

    # Isogeny Computation
    ic_method0 = mont.two_e_iso(c, n, 0)
    ic_method1 = mont.optimal_strategy(c, n, 0)
    print(f"4. IsogenyComputation")
    print(f"    NormalStrategy: {to_str(ic_method0)} = {to_str(fp.reduce(ic_method0))}")
    print(
        f"    OptimalStrategy: {to_str(ic_method1)} = {to_str(fp.reduce(ic_method1))}"
    )

    # Total cost
    tot_method0 = (gb_method0 + ld_method0 + ic_method1) / n
    tot_method1 = (gb_method0 + db_method0 + ld_method0 + ic_method1) / n
    print(f"5. TotalCost")
    print(f"    Previous: {to_str(tot_method0)} = {to_str(fp.reduce(tot_method0))}")
    print(
        f"    PreventBacktrack: {to_str(tot_method1)} = {to_str(fp.reduce(tot_method1))}"
    )


if __name__ == "__main__":
    main()
