import random
import argparse

ALPHABET = "ACGT"


def random_dna(length):
    return "".join(random.choice(ALPHABET) for _ in range(length))


def make_genome(motif, cn, flank):
    left_flank = random_dna(flank)
    right_flank = random_dna(flank)
    repeat = motif * cn
    return left_flank + repeat + right_flank


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("motif", type=str)
    parser.add_argument("cn", type=int)
    parser.add_argument("output", type=argparse.FileType("w"))

    parser.add_argument("-f", "--flank", type=int, default=10000)
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="test_genome",
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=42,
    )

    args = parser.parse_args()

    random.seed(args.seed)
    genome = make_genome(args.motif, args.cn, args.flank)

    print(f">{args.name}", file=args.output)
    for i in range(0, len(genome), 80):
        print(genome[i : i + 80], file=args.output)


if __name__ == "__main__":
    main()
