import argparse
from validate_sequences import validate_sequences

def parse_args():
    parser = argparse.ArgumentParser(description="Proyecto de alineación y distancias p")
    subparsers = parser.add_subparsers(dest="command", required=True)

    parser_align = subparsers.add_parser("validate", help="Validar secuencias FASTA")
    parser_align.add_argument("--input", default="dataset/dataset.fasta", help="Archivo FASTA de entrada")

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == "validate":
        validate_sequences(args.input)

if __name__ == "__main__":
    main()
