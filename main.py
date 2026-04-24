import argparse
from validate_sequences import validate_sequences
from p_distances import calculate_p_distances
from build_matrix import distance_matrix, write_transition_matrix, distance_from_csv
from upgma import read_distance_matrix, upgma, print_tree

def parse_args():
    parser = argparse.ArgumentParser(description="Proyecto de alineacion y distancias p")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_align = subparsers.add_parser("validate", help="Validar secuencias FASTA")
    parser_align.add_argument("--input", default="dataset/dataset.fasta", help="Archivo FASTA de entrada")

    parser_dist = subparsers.add_parser("distance", help="Calcular distancias p desde un alineamiento")
    parser_dist.add_argument("--input", default="dataset/dataset.fasta", help="Archivo FASTA alineado de entrada")
    parser_dist.add_argument("--output", default="dataset/p_distances.csv", help="Archivo CSV de distancias de salida")

    parser_matrix = subparsers.add_parser("matrix", help="Construir matriz de distancias evolutivas y de transicion")
    parser_matrix.add_argument("--input", default="dataset/p_distances.csv")
    parser_matrix.add_argument("--output", default="dataset/evo_distances.csv")

    parser_tree = subparsers.add_parser("tree", help="Construir arbol UPGMA")
    parser_tree.add_argument("--input", default="dataset/evo_distances.csv")

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == "validate":
        validate_sequences(args.input)

    elif args.command == "distance":
        calculate_p_distances(args.input, args.output)

    elif args.command == "matrix":
        distance_matrix(args.input, args.output)
        d = distance_from_csv(args.output)
        write_transition_matrix(d, "dataset/transition_matrix.csv")


    elif args.command == "tree":
        matrix = read_distance_matrix(args.input)
        tree = upgma(matrix)
        print("Árbol UPGMA:")
        print_tree(tree)

if __name__ == "__main__":
    main()