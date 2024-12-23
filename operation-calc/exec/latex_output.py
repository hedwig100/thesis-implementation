import os

NUMBER_OF_HEADER_LINE = 2
NUMBER_OF_EXPERIMENT = 5
NUMBER_OF_MSAI_LINE = 5

LENGTH_UNTIL_NAME = 7
LENGTH_UNTIL_M = 21
LENGTH_UNTIL_SAI = 12


class Cost:
    def __init__(self, m: float, s: float, a: float, i: float) -> None:
        self.m = m
        self.s = s
        self.a = a
        self.i = i

    def __add__(self, other: "Cost") -> "Cost":
        return Cost(
            self.m + other.m, self.s + other.s, self.a + other.a, self.i + other.i
        )

    def to_m(
        self,
        reduce_map: dict[str, float] = {"S": 0.8, "I": 10, "A": 0.004},
    ) -> float:
        return (
            self.m
            + reduce_map["S"] * self.s
            + reduce_map["A"] * self.a
            + reduce_map["I"] * self.i
        )


def retrieve_costs(
    filepath: str,
    number_of_algorithm: int,
) -> list[tuple[str, Cost]]:
    with open(filepath, "r") as f:
        lines = f.readlines()
    assert (
        len(lines)
        == NUMBER_OF_HEADER_LINE
        + (NUMBER_OF_HEADER_LINE + NUMBER_OF_MSAI_LINE * number_of_algorithm)
        * (NUMBER_OF_EXPERIMENT + 1)
        - 1
    )

    total_result_line = (
        NUMBER_OF_HEADER_LINE
        + (NUMBER_OF_HEADER_LINE + NUMBER_OF_MSAI_LINE * number_of_algorithm)
        * NUMBER_OF_EXPERIMENT
        + 1
    )
    results: list[tuple[str, float]] = []
    for i in range(total_result_line, len(lines), NUMBER_OF_MSAI_LINE):
        name = lines[i][LENGTH_UNTIL_NAME:].strip()
        m = float(lines[i + 1][LENGTH_UNTIL_M:].strip())
        s = float(lines[i + 2][LENGTH_UNTIL_SAI:].strip())
        a = float(lines[i + 3][LENGTH_UNTIL_SAI:].strip())
        i = float(lines[i + 4][LENGTH_UNTIL_SAI:].strip())
        results.append((name, Cost(m, s, a, i)))

    return results


def size_of_p(filepath: str) -> int:
    if "256" in filepath:
        return 256
    elif "512" in filepath:
        return 512
    elif "1024" in filepath:
        return 1024
    elif "1536" in filepath:
        return 1536
    return 0


def output_msai(results: dict[str, list[tuple[str, Cost]]]):
    for filename, costs in results.items():
        print(filename)
        for name, cost in costs:
            print(name)
            print(
                f"& ${cost.m:.0f}\mathrm{{M}} + {cost.s:.0f}\mathrm{{S}} + {cost.a:.0f}\mathrm{{A}} + {cost.i:.4f}\mathrm{{I}}$"
            )
        print()


def output_m(
    results: dict[str, list[tuple[str, Cost]]], reduce_maps: dict[int, dict[str, float]]
):
    for filename, costs in results.items():
        print(filename)
        for name, cost in costs:
            print(name)
            print(f"& ${cost.to_m(reduce_maps[size_of_p(filename)]):.2f}\mathrm{{M}}$")
        print()


def main():
    print("Dirname after 'experiment/' ?", end=" ")
    dirname = "experiment/" + input()
    print("Number of algorithm ?", end=" ")
    number_of_algorithm = int(input())

    os.chdir(os.pardir)
    filenames = sorted(os.listdir(dirname))

    reduce_maps: dict[int, dict[str, float]] = {
        256: {"A": 0.24, "S": 1.02, "I": 105.39},
        512: {"A": 0.13, "S": 1.00, "I": 118.69},
        1024: {"A": 0.050, "S": 1.00, "I": 101.73},
        1536: {"A": 0.029, "S": 1.00, "I": 86.93},
    }

    results: dict[str, list[tuple[str, Cost]]] = {}
    for filename in filenames:
        # Unused files
        if "cost" in filename:
            continue
        elif filename == "remove_tonelli":
            continue

        filepath = os.path.join(dirname, filename)
        results[filename] = retrieve_costs(filepath, number_of_algorithm)

    output_msai(results)
    output_m(results, reduce_maps=reduce_maps)


if __name__ == "__main__":
    main()
