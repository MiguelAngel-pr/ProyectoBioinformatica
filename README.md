Ejecucion del proyecto;
 1º; Instalar librerias necesarias -> install requirements.txt
 2º; Validar fasta -> python main.py validate
 3º; Calcular matriz de p-distancias apartir del alineamiento FASTA -> python main.py distance
 4º; Construir matriz de distancias evolutivas y de transición -> python main.py matrix
 5º; Construir árbol inicial -> python main.py tree
 6º; Optimizar árbol con Metropolis y generar gráficas -> python main.py montecarlo
