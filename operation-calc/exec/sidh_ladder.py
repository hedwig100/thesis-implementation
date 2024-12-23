from condition import Condition as C
from mytype import Op
import fp2
import montgomery as mont
import sympy
import parameter as params


def upside_isogenies(c: C, bit_length_of_p: Op, w: Op, a: Op) -> Op:
    prob = 1 / 2
    return (
        1 / prob * mont.xDBLe_proj(c, bit_length_of_p - a) + mont.two_e_iso(c, a)
    ) * w


def leftside_isogenies(c: C, bit_length_of_p: Op, h: Op, b: Op) -> Op:
    prob = 2 / 3
    return (
        1 / prob * mont.xTPLe_proj(c, bit_length_of_p - b) + mont.three_e_iso(c, b)
    ) * h


def sidh_ladder(c: C, bit_length_of_p: Op, a: Op, b: Op, w: Op, h: Op) -> Op:
    # random isogenies
    step1 = upside_isogenies(c, bit_length_of_p, w, a) + leftside_isogenies(
        c, bit_length_of_p, h, b
    )

    # SIDH ladder
    step2 = mont.three_e_iso(c, b) * w * h + mont.two_e_iso(c, a) * h * w

    return step1 + step2


def main():
    c = C(is_p_equal_3_mod_4=True)
    p = sympy.symbols("p")
    bit_length_of_p = 512  # sympy.log(p)
    a, b, h, w = sympy.symbols("a b h w")
    cost = sidh_ladder(c, bit_length_of_p, a, b, w, h)
    cost = sympy.simplify(cost)
    print("Cost")
    params.print_MSAI(cost)
    print()

    cost = cost.subs([(b, 890 / h), (a, 705 / w)])
    print("Cost(a*w = 705, b*h = 890)")
    params.print_MSAI(cost)
    print()

    # Currenct cost
    print("Currenct cost")
    current_cost = cost.subs([(w, 4), (h, 7)])
    params.print_MSAI(current_cost)
    print()

    # Future cost
    print("Future cost")
    future_cost = cost.subs([(w, 445), (h, 353)])
    params.print_MSAI(future_cost)


if __name__ == "__main__":
    main()
