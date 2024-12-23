import os

NUMBER_OF_HEADER_LINE = 2
NUMBER_OF_EXPERIMENT = 5
NUMBER_OF_MSAI_LINE = 5

LENGTH_UNTIL_NAME = 7
LENGTH_UNTIL_M = 21
LENGTH_UNTIL_SAI = 12


def convert_to_m(
    filepath: str,
    number_of_algorithm: int,
    reduce_map: dict[str, float] = {"S": 0.8, "I": 10, "A": 0.004},
) -> list[tuple[str, float]]:
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
        cost = m + reduce_map["S"] * s + reduce_map["A"] * a + reduce_map["I"] * i
        results.append((name, cost))

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

    results: dict[str, list[tuple[str, float]]] = {}
    for filename in filenames:
        # Unused files
        if "cost" in filename:
            continue
        elif filename == "remove_tonelli":
            continue

        filepath = os.path.join(dirname, filename)
        result = convert_to_m(
            filepath, number_of_algorithm, reduce_map=reduce_maps[size_of_p(filename)]
        )
        results[filename] = result

    # Print costs
    print("Costs")
    for filename, result in results.items():
        print(filename)
        for name, cost in result:
            print(name, ": ", cost)
        print()

    # Compute accelration ratio in each security level
    print("Accelration ratio")
    for bit_length_of_p in [256, 512, 1024, 1536]:
        print(f"bit_length_of_p: {bit_length_of_p}")
        # Compute minimum cost of each function
        min_costs: dict[str, float] = {}
        for filename, result in results.items():
            if str(bit_length_of_p) not in filename:
                continue
            for name, cost in result:
                if name not in min_costs:
                    min_costs[name] = cost
                elif cost < min_costs[name]:
                    min_costs[name] = cost

        # Compute accelration ratio
        max_cost = min_costs["ModularPolynomial-CGL"]
        for name, cost in min_costs.items():
            if max_cost == min_costs[name]:
                print(f"{name}: 1.0")
            else:
                print(f"{name}: {max_cost / cost:.4f}")
        print()

    # Compute cost reduction from DPB-CGL-prevent-backtrack
    print("Cost Reduction from 'DPB-CGL-prevent-backtrack'")
    for bit_length_of_p in [256, 512, 1024, 1536]:
        print(f"bit_length_of_p: {bit_length_of_p}")
        # Compute minimum cost of each function
        min_costs: dict[str, float] = {}
        for filename, result in results.items():
            if str(bit_length_of_p) not in filename:
                continue
            for name, cost in result:
                if name not in min_costs:
                    min_costs[name] = cost
                elif cost < min_costs[name]:
                    min_costs[name] = cost

        # Compute accelration ratio
        dpb_eprint_cost = min_costs["DPB-CGL-prevent-backtrack"]
        for name, cost in min_costs.items():
            if max_cost == min_costs[name]:
                print(f"{name}: 1.0")
            else:
                print(f"{name}: {cost / dpb_eprint_cost:.4f}")
        print()


if __name__ == "__main__":
    main()
