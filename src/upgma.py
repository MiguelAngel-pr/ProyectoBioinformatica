import csv
import random


class Node:
    def __init__(self, name, left=None, right=None, height=0.0):
        self.name = name
        self.left = left
        self.right = right
        self.height = height

    def is_leaf(self):
        return self.left is None and self.right is None


def read_distance_matrix(file_path: str):
    with open(file_path, newline="", encoding="utf-8") as f:
        reader = list(csv.reader(f))

    names = reader[0][1:]
    matrix = {}

    for i in range(1, len(reader)):
        row_name = reader[i][0]
        matrix[row_name] = {}
        for j in range(1, len(reader[i])):
            matrix[row_name][names[j - 1]] = float(reader[i][j])

    return matrix


def upgma(distance_matrix: dict):
    clusters = {name: Node(name) for name in distance_matrix.keys()}

    # Copia profunda de la matriz para no modificar el original
    distances = {a: dict(row) for a, row in distance_matrix.items()}

    while len(clusters) > 1:
        # --- Encontrar el par con distancia mínima ---
        min_dist = float("inf")
        pair = None

        cluster_keys = list(clusters.keys())
        for i in range(len(cluster_keys)):
            for j in range(i + 1, len(cluster_keys)):  # solo triángulo superior → evita duplicados
                a, b = cluster_keys[i], cluster_keys[j]
                d = distances[a][b]
                if d < min_dist:
                    min_dist = d
                    pair = (a, b)

        a, b = pair

        # --- Guardar claves ANTES de modificar clusters ---
        # Excluimos a y b; el nuevo nodo todavía no existe
        remaining = [k for k in clusters if k != a and k != b]

        # --- Crear nuevo nodo interno ---
        new_name = f"({a},{b})"
        new_node = Node(new_name, clusters[a], clusters[b], height=min_dist / 2)

        # --- Actualizar distancias al nuevo clúster ---
        distances[new_name] = {}
        for k in remaining:
            dist = (distances[a][k] + distances[b][k]) / 2
            distances[new_name][k] = dist
            distances[k][new_name] = dist
        distances[new_name][new_name] = 0.0

        # --- Limpiar entradas antiguas ---
        distances.pop(a, None)
        distances.pop(b, None)
        for k in distances:
            distances[k].pop(a, None)
            distances[k].pop(b, None)

        # --- Actualizar clusters ---
        del clusters[a]
        del clusters[b]
        clusters[new_name] = new_node

    return list(clusters.values())[0]

def worst_tree_extremo(distance_matrix: dict):
    species = list(distance_matrix.keys())
    # Aleatoriedad real
    random.shuffle(species)

    # Empezamos con la primera especie aleatoria
    current = Node(species[0], height=0)

    for i in range(1, len(species)):
        new_leaf = Node(species[i], height=0)
        
        # Distancia aleatoria o extrema para peor visualización
        dist_extrema = random.uniform(10, 100) 

        # Creamos un ancestro común que siempre se inclina a un lado
        current = Node(
            name=f"ancestro_{i}",
            left=current,
            right=new_leaf,
            height=dist_extrema
        )

    return current


def print_tree(node, level=0):
    indent = "  " * level
    if node.is_leaf():
        print(f"{indent}{node.name}  (h={node.height:.4f})")
    else:
        print(f"{indent}[nodo interno h={node.height:.4f}]")
    if node.left:
        print_tree(node.left, level + 1)
    if node.right:
        print_tree(node.right, level + 1)


if __name__ == "__main__":
    matrix = read_distance_matrix("dataset/jc_distances.csv")
    tree = upgma(matrix)
    print("Árbol UPGMA:")
    print_tree(tree)