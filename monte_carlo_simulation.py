import csv
import math
import random
import copy
import matplotlib
# Configuración para que no pete en entornos sin GUI
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from upgma import Node

# LECTURA DE DATOS
def read_all_data(file_path):
    distancias_dict = {}
    matrices_dict = {}
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            data = f.read()
        
        # Separamos por bloques usando el separador 'Par ' que vimos en el CSV
        bloques = data.split("Par ")
        for b in bloques:
            if "vs" in b and "d =" in b:
                lineas = [l.strip() for l in b.strip().split("\n") if l.strip()]
                cabecera = lineas[0]
                
                # Parseo manual de la línea
                nombres = cabecera.split(":")[1].split("(")[0].strip()
                sp1, sp2 = [s.strip() for s in nombres.split(" vs ")]
                d_val = float(cabecera.split("d =")[1].replace(")", "").strip())
                
                par_ordenado = tuple(sorted([sp1, sp2]))
                distancias_dict[par_ordenado] = d_val

                # Buscamos la matriz de frecuencias que viene debajo
                try:
                    for idx, linea in enumerate(lineas):
                        if "A,C,G,T" in linea:
                            m = {}
                            # Leemos las 4 filas (A, C, G, T)
                            for offset in range(1, 5):
                                row = lineas[idx + offset].split(",")
                                m[row[0]] = {
                                    "A": float(row[1]), "C": float(row[2]),
                                    "G": float(row[3]), "T": float(row[4])
                                }
                            matrices_dict[par_ordenado] = m
                            break
                except Exception:
                    continue 
        return distancias_dict, matrices_dict
    except Exception as e:
        print(f"Error al leer el CSV: {e}")
        return {}, {}

# LÓGICA DE DISTANCIAS EN EL ÁRBOL
def obtener_hojas(nodo):
    if nodo.is_leaf():
        return [nodo.name]
    return obtener_hojas(nodo.left) + obtener_hojas(nodo.right)

def encontrar_lca(nodo, a, b):
    # Buscamos el ancestro común más cercano (LCA) para medir la distancia
    def buscar(n, t):
        if n is None: return False
        if n.name == t: return True
        return buscar(n.left, t) or buscar(n.right, t)

    if nodo is None: return None
    if nodo.name == a or nodo.name == b: return nodo
    
    izq_tiene = buscar(nodo.left, a) or buscar(nodo.left, b)
    der_tiene = buscar(nodo.right, a) or buscar(nodo.right, b)
    
    if izq_tiene and der_tiene: return nodo
    return encontrar_lca(nodo.left, a, b) if izq_tiene else encontrar_lca(nodo.right, a, b)

def distancia_al_ancestro(nodo, target, acumulado=0):
    if nodo is None: return None
    if nodo.name == target: return acumulado
    
    # Probamos por la izquierda
    res = distancia_al_ancestro(nodo.left, target, acumulado + (nodo.left.height if nodo.left else 0))
    if res is not None: return res
    
    # Si no, por la derecha
    return distancia_al_ancestro(nodo.right, target, acumulado + (nodo.right.height if nodo.right else 0))

def distancia_entre_hojas(arbol, a, b):
    ancestro = encontrar_lca(arbol, a, b)
    if not ancestro: return 0
    return distancia_al_ancestro(ancestro, a) + distancia_al_ancestro(ancestro, b)

# FUNCIÓN DE COSTO (SCORE)
def score(arbol, d_reales, m_probs):
    hojas = obtener_hojas(arbol)
    error_acumulado = 0.0
    num_pares = 0
    
    for i in range(len(hojas)):
        for j in range(i + 1, len(hojas)):
            p = tuple(sorted([hojas[i], hojas[j]]))
            if p in d_reales:
                # Distancia del árbol actual
                d_calc = distancia_entre_hojas(arbol, hojas[i], hojas[j])
                
                # Probabilidad Jukes-Cantor
                p_tree = 0.25 + 0.75 * math.exp(-(4/3) * d_calc)
                
                # Probabilidad real del CSV
                m = m_probs[p]
                p_real = (m["A"]["A"] + m["C"]["C"] + m["G"]["G"] + m["T"]["T"]) / 4
                
                # Sumamos el error cuadrático de distancia y de probabilidad
                error_acumulado += (d_calc - d_reales[p])**2 + (p_tree - p_real)**2
                num_pares += 1
                
    return error_acumulado / num_pares if num_pares > 0 else 999.0

# MUTACIONES (SPR e INSERT/REMOVE)
def remove_node(nodo, objetivo):
    if nodo is None or nodo.name == objetivo:
        return None
    
    if nodo.left and nodo.left.name == objetivo:
        return nodo.right # El hermano sube de nivel
    if nodo.right and nodo.right.name == objetivo:
        return nodo.left # El hermano sube de nivel
        
    nodo.left = remove_node(nodo.left, objetivo)
    nodo.right = remove_node(nodo.right, objetivo)
    return nodo

def insert_node(actual, target_name):
    # Altura aleatoria para la nueva rama (ajustada para JC)
    nueva_h = random.uniform(0.01, 0.35)
    
    # Decidimos si insertamos aquí o seguimos bajando
    if actual.is_leaf() or random.random() < 0.3:
        return Node("nodo_int", actual, Node(target_name), height=nueva_h)
    
    if random.random() < 0.5:
        actual.left = insert_node(actual.left, target_name)
    else:
        actual.right = insert_node(actual.right, target_name)
    return actual

def generate_neighbor(tree):
    # Hacemos copia profunda para no modificar el original antes de tiempo
    nuevo_arbol = copy.deepcopy(tree)
    
    # Tenemos dos tipos de mutación: de alturas (70%) o de estructura (30%)
    if random.random() < 0.7:
        def mutar_alturas(n):
            if n:
                # Pequeño ajuste a la rama
                n.height += random.uniform(-0.05, 0.05)
                if n.height < 0.001: n.height = 0.001 # No puede ser negativa
                if not n.is_leaf():
                    mutar_alturas(n.left)
                    mutar_alturas(n.right)
        mutar_alturas(nuevo_arbol)
    else:
        hojas = obtener_hojas(nuevo_arbol)
        if len(hojas) > 2:
            sp = random.choice(hojas)
            arbol_podado = remove_node(nuevo_arbol, sp)
            if arbol_podado:
                nuevo_arbol = insert_node(arbol_podado, sp)
                
    return nuevo_arbol

# ALGORITMO PRINCIPAL (METROPOLIS)
def metropolis(tree, d_reales, m_probs, iteraciones=100000, beta=500.0, paciencia=5000):
    actual = tree
    score_actual = score(actual, d_reales, m_probs)
    
    mejor_arbol = copy.deepcopy(actual)
    mejor_score = score_actual
    
    historial_scores = []
    espera = 0
    ultimos_mejores = [] # Para el historial de los 5 mejores

    for i in range(iteraciones):
        vecino = generate_neighbor(actual)
        score_vecino = score(vecino, d_reales, m_probs)
        
        diff = score_vecino - score_actual
        
        # Criterio de aceptación de Metrópolis
        if diff < 0 or random.random() < math.exp(-beta * diff):
            actual = vecino
            score_actual = score_vecino
            
        historial_scores.append(score_actual)
        
        # Si encontramos una mejora, actualizamos y reseteamos paciencia
        if score_actual < mejor_score:
            mejor_score = score_actual
            mejor_arbol = copy.deepcopy(actual)
            espera = 0
            
            # Guardamos para el historial de variantes
            ultimos_mejores.append((copy.deepcopy(mejor_arbol), mejor_score))
            if len(ultimos_mejores) > 5:
                ultimos_mejores.pop(0)
        else:
            espera += 1
            
        if espera >= paciencia:
            print(f"\n[Terminado por paciencia en iteración {i}]")
            break
            
        if i % 10000 == 0:
            print(f"Iteración {i} | Score actual: {score_actual:.6f}")
            
    return mejor_arbol, mejor_score, historial_scores, ultimos_mejores

def graficar_convergencia(historial):
    """
    Genera una gráfica de convergencia optimizada para ser visible en 
    simulaciones largas (muestreo n=100).
    """
    plt.figure(figsize=(10, 6))
    
    # Si graficamos 1.000.000 de puntos no se verá nada.
    # Aplicamos un salto (step) de 100 para limpiar el ruido visual.
    paso = 100
    historial_reducido = historial[::paso]
    eje_x = list(range(0, len(historial), paso))
    
    # 1. Dibujamos el historial con transparencia para ver la exploración
    plt.plot(eje_x, historial_reducido, color='#3498db', alpha=0.4, 
             label=f'Score Muestreado (n={paso})', linewidth=0.8)
    
    # 2. Calculamos media móvil sobre la muestra para resaltar la convergencia
    window = max(1, len(historial_reducido) // 50)
    suavizado = [sum(historial_reducido[max(0, i-window):i+1]) / (min(i, window)+1) 
                 for i in range(len(historial_reducido))]
    
    plt.plot(eje_x, suavizado, color='#e74c3c', linewidth=2.5, label='Tendencia Global')
    
    # Configuración de ejes (Escala logarítmica es clave para ver mejoras pequeñas)
    plt.yscale('log') 
    plt.title("Análisis de Convergencia: Optimización del Árbol Filogenético", fontsize=14)
    plt.xlabel("Número de Iteraciones (Simulación Monte Carlo)", fontsize=12)
    plt.ylabel("Logaritmo del Error (Score)", fontsize=12)
    
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(loc='upper right')
    
    plt.tight_layout()
    
    plt.savefig("dataset/convergencia_mc.png")
    print(f"\n Gráfica de convergencia guardada como convergencia_mc.png.")
    plt.close()

def graficar_distribucion_scores(historial):
    plt.figure(figsize=(10, 6))
    plt.hist(historial, bins=50, color='#9b59b6', edgecolor='white', alpha=0.7)
    plt.title("Distribución de Scores Aceptados")
    plt.xlabel("Valor del Score (Error)")
    plt.ylabel("Frecuencia (Iteraciones)")
    plt.grid(axis='y', alpha=0.3)
    plt.savefig("dataset/distribucion_scores.png")
    print(f"\n Gráfica de distribución de scores guardada como distribucion_scores.png.")
    plt.close()

def get_tree_distances(node):
    def find_leaves(n):
        if n.is_leaf(): return [n]
        return find_leaves(n.left) + find_leaves(n.right)

    hojas = find_leaves(node)
    especies = [h.name for h in hojas]
    n_sp = len(especies)
    distancias_arbol = {}
    
    for i in range(n_sp):
        for j in range(i + 1, n_sp):
            sp1, sp2 = hojas[i].name, hojas[j].name
            # Usa la función que ya tienes para medir distancias reales del árbol
            d = distancia_entre_hojas(node, sp1, sp2) 
            llave = tuple(sorted([sp1, sp2]))
            distancias_arbol[llave] = d
            
    return distancias_arbol, especies

def graficar_comparativa_calor(matriz_real_dict, arbol_final):
    # Crea dos mapas de calor: uno con los datos originales y otro con lo que el árbol 'logró' representar.

    # Obtener distancias del árbol y nombres de especies
    dist_arbol_dict, especies = get_tree_distances(arbol_final)
    especies = sorted(especies)
    n = len(especies)
    
    # Inicializar matrices numéricas
    m_real = np.zeros((n, n))
    m_arbol = np.zeros((n, n))
    
    for i, sp1 in enumerate(especies):
        for j, sp2 in enumerate(especies):
            if i == j: continue
            llave = tuple(sorted([sp1, sp2]))
            m_real[i, j] = matriz_real_dict.get(llave, 0)
            m_arbol[i, j] = dist_arbol_dict.get(llave, 0)

    # Graficas
    fig,(ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    sns.heatmap(m_real, annot=True, fmt=".3f", cmap="YlGnBu", 
                xticklabels=especies, yticklabels=especies, ax=ax1)
    ax1.set_title("Distancias Reales (Entrada)")
    
    sns.heatmap(m_arbol, annot=True, fmt=".3f", cmap="YlOrRd", 
                xticklabels=especies, yticklabels=especies, ax=ax2)
    ax2.set_title("Distancias del Árbol (Optimizado)")
    
    plt.tight_layout()
    plt.savefig("dataset/comparativa_calor.png")
    print("\n Mapa de calor guardado como 'comparativa_calor.png'")
    plt.close()

def worst_tree(especies):
    #Crea un árbol inicial 'malo' (tipo peine/degenerado), para que el Monte Carlo tenga de donde empezar a optimizar.
    
    # Empezamos con la primera especie como raíz/hoja
    arbol = Node(especies[0])
    
    # Vamos colgando el resto de especies una a una para crear un árbol muy desequilibrado
    for i in range(1, len(especies)):
        # Creamos un nodo intermedio que une el árbol que ya tenemos con la nueva especie
        # Ponemos una altura inicial genérica (0.1)
        arbol = Node("nodo_int", arbol, Node(especies[i]), height=0.1)
        
    return arbol

def print_tree(node, level=0):
    indent = "  " * level
    print(f"{indent}{node.name} (h={node.height:.4f})")
    if node.left: print_tree(node.left, level + 1)
    if node.right: print_tree(node.right, level + 1)

if __name__ == "__main__":
    ruta_csv = "dataset/transition_matrix.csv" 
    matrix_datos = read_all_data(ruta_csv)
    
    if matrix_datos:
        especies = list(set([sp for par in matrix_datos.keys() for sp in par]))
        print(f"Especies detectadas: {especies}")
        tree_ini = worst_tree(especies)
        
        print("\n--- ÁRBOL INICIAL ---")
        print_tree(tree_ini)
        print(f"Score inicial: {score(tree_ini, matrix_datos):.6f}")
        
        # Subimos Beta para que sea más estricto con el error cuadrático
        best_tree, best_score = metropolis(tree_ini, matrix_datos, beta=500.0)
        
        print("\n--- ÁRBOL OPTIMIZADO ---")
        print_tree(best_tree)
        print(f"Score final: {best_score:.6f}")
    else:
        print("No se encontraron datos. Verifica el formato del CSV.")