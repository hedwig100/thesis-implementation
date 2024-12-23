import os

NUMBER_OF_HEADER_LINE = 2
NUMBER_OF_EXPERIMENT = 5
NUMBER_OF_MSAI_LINE = 5

LENGTH_UNTIL_NAME = 7
LENGTH_UNTIL_M = 21
LENGTH_UNTIL_SAI = 12

N_EXPERIMENT = 5


def compute_mean(
    filepath: str,
    number_of_algorithm: int,
) -> dict[str, list[float]]:
    with open(filepath, "r") as f:
        lines = f.readlines()
    if (
        len(lines)
        != NUMBER_OF_HEADER_LINE
        + (NUMBER_OF_HEADER_LINE + NUMBER_OF_MSAI_LINE * number_of_algorithm)
        * NUMBER_OF_EXPERIMENT
    ):
        # Dont have to compute mean
        return {}

    results: dict[str, list[float]] = {}
    for exp_id in range(N_EXPERIMENT):
        start_line = (
            NUMBER_OF_HEADER_LINE
            + (NUMBER_OF_HEADER_LINE + NUMBER_OF_MSAI_LINE * number_of_algorithm)
            * exp_id
            + NUMBER_OF_HEADER_LINE
        )
        for j in range(number_of_algorithm):
            name = lines[start_line][LENGTH_UNTIL_NAME:].strip()
            m = float(lines[start_line + 1][LENGTH_UNTIL_M:].strip())
            s = float(lines[start_line + 2][LENGTH_UNTIL_SAI:].strip())
            a = float(lines[start_line + 3][LENGTH_UNTIL_SAI:].strip())
            i = float(lines[start_line + 4][LENGTH_UNTIL_SAI:].strip())
            start_line += NUMBER_OF_MSAI_LINE

            if name in results:
                results[name][0] += m
                results[name][1] += s
                results[name][2] += a
                results[name][3] += i
            else:
                results[name] = [m, s, a, i]

    return {name: [c / N_EXPERIMENT for c in costs] for name, costs in results.items()}


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

    results: dict[str, list[tuple[str, float]]] = {}
    for filename in filenames:
        # Unused files
        if "cost" in filename:
            continue
        elif filename == "remove_tonelli":
            continue

        filepath = os.path.join(dirname, filename)
        result = compute_mean(filepath, number_of_algorithm)

        print(filepath)
        print("Total Results:")
        for name, cost in result.items():
            print(f"  name: {name}")
            print(f"  Count:   M_COUNTER: {cost[0]}")
            print(f"  S_COUNTER: {cost[1]}")
            print(f"  A_COUNTER: {cost[2]}")
            print(f"  I_COUNTER: {cost[3]}")
        print()


if __name__ == "__main__":
    main()
