import math
import csv


def jukes_cantor(p: float) -> float:
    if p >= 0.75: #Si p >= 0.75 la distancia es incalculable
        return float("inf")
    return -3/4 * math.log(1 - (4/3)*p)

def distance_from_csv(input_csv: str) -> float: # Calcula la media de las distancias evolutivas del csv para usarla en la matriz de transicion
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = list(csv.reader(f))

    values = []

    for i in range(1, len(reader)):
        for j in range(1, len(reader[i])):
            d = float(reader[i][j])
            if d != float("inf") and d != 0.0: # Excluimos los valores de la diagonal y valores no estimables (inf)
                values.append(d)

    return sum(values) / len(values)

def distance_matrix(input_csv: str, output_csv: str) -> None: 
    with open(input_csv, newline="", encoding="utf-8") as f:
        reader = list(csv.reader(f))

    header = reader[0]
    names = header[1:]

    evo_matrix = [header]

    for i in range(1, len(reader)):
        row_name = reader[i][0]
        new_row = [row_name]

        for j in range(1, len(reader[i])):
            p = float(reader[i][j])
            
            # Aplicamos la correccion de Jukes-Cantor a cada p-distancia
            d = jukes_cantor(p)
            new_row.append(f"{d:.6f}")

        evo_matrix.append(new_row)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows(evo_matrix)

    print(f"Matriz de distancias evolutivas escrita en: {output_csv}")

def transition_matrix(d: float):
    import math

    exp_term = math.exp((-4/3) * d)

    p_same = 0.25 + 0.75 * exp_term # probabilidad de permanecer en la misma base
    p_change = 0.25 - 0.25 * exp_term # probabilidad de cambiar a una base concreta

    bases = ["A", "C", "G", "T"]
    matrix = {}

    for b1 in bases: #dependiendo de si esta fuera o dentro de la diagonal, la base cambia o no.
        matrix[b1] = {}
        for b2 in bases:
            if b1 == b2:
                matrix[b1][b2] = p_same
            else:
                matrix[b1][b2] = p_change

    return matrix

def write_transition_matrix(d: float, output_csv: str):
    matrix = transition_matrix(d)
    bases = ["A", "C", "G", "T"]

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        # cabecera
        writer.writerow([""] + bases)

        # filas
        for b1 in bases:
            row = [b1] + [f"{matrix[b1][b2]:.6f}" for b2 in bases]
            writer.writerow(row)

    print(f"Matriz de transicion escrita en: {output_csv}")

if __name__ == "__main__":
    distance_matrix("dataset/p_distances.csv", "dataset/evo_distances.csv")
    d = distance_from_csv("dataset/evo_distances.csv")
    write_transition_matrix(d, "dataset/transition_matrix.csv")