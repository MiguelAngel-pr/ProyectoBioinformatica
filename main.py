import argparse
from validate_sequences import validate_sequences
from p_distances import calculate_p_distances

def parse_args():
    parser = argparse.ArgumentParser(description="Proyecto de alineacion y distancias p")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_align = subparsers.add_parser("validate", help="Validar secuencias FASTA")
    parser_align.add_argument("--input", default="dataset/dataset.fasta", help="Archivo FASTA de entrada")

    parser_dist = subparsers.add_parser("distance", help="Calcular distancias p desde un alineamiento")
    parser_dist.add_argument("--input", default="dataset/dataset.fasta", help="Archivo FASTA alineado de entrada")
    parser_dist.add_argument("--output", default="dataset/p_distances.csv", help="Archivo CSV de distancias de salida")

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == "validate":
        validate_sequences(args.input)
    elif args.command == "distance":
        calculate_p_distances(args.input, args.output)

if __name__ == "__main__":
    main()
