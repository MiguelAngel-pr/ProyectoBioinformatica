import csv


class Node:
    def __init__(self, name, left=None, right=None, height=0.0):
        self.name = name
        self.left = left
        self.right = right
        self.height = height

    def is_leaf(self):
        return self.left is None and self.right is None


def read_distance_matrix(file_path: str): #a partir del csv forma la matriz de distancias evolutivas para usarla en el algoritmo UPGMA
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


def upgma(distance_matrix: dict): #Construye un arbol filogenetico jerarquico uniendo en cada paso los dos clusteres mas cercanos.
    
    # Empezamos con un nodo hoja por cada secuencia
    clusters = {name: Node(name) for name in distance_matrix.keys()}

    # Copia profunda de la matriz para no modificar el original
    distances = {a: dict(row) for a, row in distance_matrix.items()}

    while len(clusters) > 1:
        # Encontramos el par con distancia mínima
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

        # Guardamos claves antes de modificar clusters
        remaining = [k for k in clusters if k != a and k != b]

        # Creamos nuevo nodo interno que une a y b
        new_name = f"({a},{b})"
        new_node = Node(new_name, clusters[a], clusters[b], height=min_dist / 2)

        # Actualizamos distancias al nuevo clúster
        distances[new_name] = {}
        for k in remaining:
            dist = (distances[a][k] + distances[b][k]) / 2
            distances[new_name][k] = dist
            distances[k][new_name] = dist

        distances[new_name][new_name] = 0.0

        # Limpiamos entradas antiguas, quitando a y b
        distances.pop(a, None)
        distances.pop(b, None)
        for k in distances:
            distances[k].pop(a, None)
            distances[k].pop(b, None)

        # Actualizamos los clusters
        del clusters[a]
        del clusters[b]
        clusters[new_name] = new_node

    return list(clusters.values())[0]


def print_tree(node, level=0):
    indent = "  " * level
    if node.is_leaf(): #Nombre Secuencia
        print(f"{indent}{node.name}  (h={node.height:.4f})")
    else: #Altura del nodo interno
        print(f"{indent}[nodo interno h={node.height:.4f}]")

    # Recorremos recursivamente los hijos
    if node.left:
        print_tree(node.left, level + 1)
    if node.right:
        print_tree(node.right, level + 1)


if __name__ == "__main__":
    matrix = read_distance_matrix("dataset/evo_distances.csv")
    tree = upgma(matrix)
    print("Árbol UPGMA:")
    print_tree(tree)