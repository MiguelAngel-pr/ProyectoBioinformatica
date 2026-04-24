from Bio import SeqIO

def validate_sequences(input_fasta: str) -> None:
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if len(records) == 0: 
        raise ValueError(f"No se encontraron secuencias en {input_fasta}")

    lengths = [len(rec.seq) for rec in records]
    if len(set(lengths)) != 1:
        raise ValueError("Las secuencias no tienen la misma longitud")

    valid_bases = set("ACGT")
    for rec in records: # Bases diferentes a A, C, G, T
        if not set(str(rec.seq)).issubset(valid_bases):
            raise ValueError(f"Secuencia {rec.id} contiene caracteres invalidos")

    if len(set(lengths)) == 1: # Todas las secuencias tienen la misma longitud
        print(f"Archivo validado")
        return

    if __name__ == "__main__":
        validate_sequences("dataset/dataset.fasta")