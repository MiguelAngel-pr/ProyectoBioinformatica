import argparse
from validate_sequences import validate_sequences
from p_distances import calculate_p_distances
from build_matrix import distance_matrix, write_transition_matrix
from upgma import read_distance_matrix, upgma, print_tree, worst_tree_extremo
from monte_carlo_simulation import graficar_distribucion_scores, metropolis, worst_tree, read_all_data, graficar_convergencia, graficar_comparativa_calor

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

    parser_tree = subparsers.add_parser("tree", help="Construir árbol UPGMA")
    parser_tree.add_argument("--input", default="dataset/evo_distances.csv")

    parser_mc = subparsers.add_parser("montecarlo", help="Optimizar árbol con Metropolis")
    parser_mc.add_argument("--input", default="dataset/evo_distances.csv")
    parser_mc.add_argument("--matrices", default="dataset/transition_matrix.csv")
    parser_mc.add_argument("--iter", type=int, default=1000000)
    parser_mc.add_argument("--beta", type=float, default=0.1)
    parser_mc.add_argument("--paciencia", type=int, default=50000)

    return parser.parse_args()

def main():
    args = parse_args()
    
    if args.command == "validate":
        validate_sequences(args.input)

    elif args.command == "distance":
        calculate_p_distances(args.input, args.output)

    elif args.command == "matrix":
        distance_matrix(args.input, args.output)
        write_transition_matrix("dataset/evo_distances.csv", "dataset/transition_matrix.csv")

    elif args.command == "tree":
        matrix = read_distance_matrix(args.input)
        tree = upgma(matrix)
        print("Árbol UPGMA:")
        print_tree(tree)

    elif args.command == "montecarlo":
        matrix_dist_bloques, matrix_prob = read_all_data(args.matrices)
        matrix_tabla = read_distance_matrix(args.input)

        #nombres_especies = list(matrix_tabla.keys())
        
        #tree_ini = worst_tree(nombres_especies)

        tree_ini = worst_tree_extremo(matrix_tabla) 

        #tree_ini = upgma(matrix_tabla)

        # Ejecutamos la simulación recuperando los nuevos valores
        best_tree, best_score, historial, ultimos_5 = metropolis(
            tree_ini,
            matrix_dist_bloques,
            matrix_prob,
            iteraciones=args.iter,
            beta=args.beta,
            paciencia=args.paciencia
        )

        # --- MOSTRAR LOS ÚLTIMOS 5 ÁRBOLES ---
        print("\n" + "="*40)
        print("HISTORIAL: ÚLTIMOS 5 MEJORES ÁRBOLES")
        print("="*40)
        for idx, (arbol, s) in enumerate(ultimos_5):
            print(f"\n[Variante {idx+1}] Score: {s:.6f}")
            print_tree(arbol)

        print("\n" + "="*40)
        print("RESULTADO FINAL")
        print("="*40)
        print_tree(best_tree)
        print(f"\nMejor Score: {best_score:.6f}")

        # --- GENERAR GRÁFICA ---
        graficar_convergencia(historial)
        graficar_distribucion_scores(historial)
        graficar_comparativa_calor(matrix_dist_bloques, best_tree)

if __name__ == "__main__":
    main()