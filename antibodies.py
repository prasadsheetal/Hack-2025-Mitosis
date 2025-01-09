from Bio import SeqIO
from rich.console import Console
from rich.text import Text

EXAMPLE_ANTIBODIES = [
    # placeholder
    "ATGCTAGCTAGCTAGCTA",
    "TCGATCGTAGCTAGATCT",
    "GCTAGCTCGATCGATCGA",
]




def search_and_highlight(query, fasta_file):
    # Initialize a console for printing
    console = Console()

    # Read sequences from the FASTA file
    sequences = SeqIO.parse(fasta_file, "fasta")

    # Loop through each sequence
    for record in sequences:
        sequence_str = str(record.seq)
        if query in sequence_str:
            # Highlight the matching part in the sequence
            highlighted_text = Text(sequence_str)
            start = 0
            while True:
                start = sequence_str.find(query, start)
                if start == -1:
                    break
                # Add highlight for each match
                highlighted_text.stylize("bold red", start, start + len(query))
                start += len(query)  # Move past this match

            # Print the result
            console.print(f"> {record.id} {record.description}")
            console.print(highlighted_text)


# Example usage
search_query = input("Enter the sequence query: ").upper()
fasta_file_path = "sequences.fasta"
search_and_highlight(search_query, fasta_file_path)


''' reference code from the demo

def smiles_to_svg(smiles: str, width: int = 400, height: int = 400) -> bytes:
    """
    makes an SVG image of a molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES")

    Chem.rdCoordGen.AddCoords(mol)
    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    # set drawing options on drawer.getOptions()
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()


def get_substructure_fingerprint(mol):
    """
    TODO: substructure fingerprints
    returns substructure fingerprint for mol
    """
    fp = ExplicitBitVect(SUBSTRUCTURE_FP_SIZE, True)
    return fp  # this is currently empty and useless


def get_highlighted_image(
    target_smiles: str, query_smiles: str, width: int = 400, height: int = 400
):
    """
    TODO: Creates image of highlighted molecule
    """
    target_mol = Chem.MolFromSmiles(target_smiles)
    query_mol = Chem.MolFromSmiles(query_smiles)
    match = target_mol.GetSubstructMatch(query_mol)

    Chem.rdCoordGen.AddCoords(target_mol)
    Chem.rdCoordGen.AddCoords(query_mol)

    drawer = Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(target_mol, highlightAtoms=match)
    drawer.FinishDrawing()
    return drawer.GetDrawingText().encode()


def search_compounds(
    query_smiles: str, compound_list: List[str] = EXAMPLE_COMPOUNDS
) -> List[List[int]]:
    """
    search the list of smiles and return substructure match indices

    is empty list if not match and that index

    it would be nice to fingerprint and store fingerprints in memory
    """
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        raise RuntimeError("Invalid query SMILES")

    compounds = [Chem.MolFromSmiles(s) for s in compound_list]
    matches = []
    for m in compounds:
        if m is None:
            matches.append([])
        else:
            matches.append(m.GetSubstructMatch(query_mol))
    return matches

'''
# print(search_compounds("C"))
# print(search_compounds("C1CCCCC1"))

# import rdkit
# print(rdkit.__version__)

# print(smiles_to_svg("CC"))


#function to assign locations to different parts of the antibodies..