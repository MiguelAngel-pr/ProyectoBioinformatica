import csv
from itertools import combinations
from Bio import SeqIO


def p_distance(seq_a: str, seq_b: str) -> float:
    if len(seq_a) != len(seq_b):
        raise ValueError("Las secuencias deben tener la misma longitud alineada")
    differences = 0
    comparable = 0
    for a, b in zip(seq_a, seq_b):
        if a == "-" or b == "-":
            continue
        comparable += 1
        if a != b:
            differences += 1
    if comparable == 0:
        return 0.0
    return differences / comparable


def calculate_p_distances(input_fasta: str, output_csv: str) -> None:
    records = list(SeqIO.parse(input_fasta, "fasta"))
    # Comprobamos la validez de las secuencias
    if len(records) < 2:
        raise ValueError("Se necesitan al menos dos secuencias alineadas para calcular distancias")

    lengths = {len(rec.seq) for rec in records}
    if len(lengths) != 1:
        raise ValueError("Todas las secuencias deben tener la misma longitud en el archivo alineado")

    # Formamos la matriz de distancias p
    names = [rec.id for rec in records]
    distances = {name: {} for name in names}

    for a, b in combinations(records, 2):
        dist = p_distance(str(a.seq), str(b.seq))
        distances[a.id][b.id] = dist
        distances[b.id][a.id] = dist

    for name in names:
        distances[name][name] = 0.0

    with open(output_csv, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequence"] + names)
        for name in names:
            row = [name] + [f"{distances[name][other]:.6f}" for other in names]
            writer.writerow(row)

    print(f"Matriz de distancias p escrita en: {output_csv}")

if __name__ == "__main__":
    calculate_p_distances("dataset/dataset.fasta", "dataset/p_distances.csv")
