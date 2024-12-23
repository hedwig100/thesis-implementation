import sympy
from condition import Condition
import mycgl
import parameter as param


def main():
    c = Condition(is_p_equal_3_mod_4=True)
    p = sympy.symbols("p")
    bit_length_of_p = 256  # sympy.log(p)
    k = sympy.symbols("k")
    cost = mycgl.sidh_cgl(c, bit_length_of_p, k)
    cost = sympy.simplify(cost)
    print(f"k: {k}, cost: {cost}")
    for k in range(2, 100, 2):
        cost = mycgl.sidh_cgl(c, bit_length_of_p, k)
        cost = sympy.simplify(cost)
        print(f"k: {k}, cost: {cost}")


if __name__ == "__main__":
    main()
