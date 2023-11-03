def read_file(name):
    """Read the contents of a file from the 'data' folder"""
    with open(f"data/{name}") as f:
        return f.read().strip()

def read_FASTA(name):
    """Read the lines of a FASTA file from the 'data' folder"""
    entries = read_file(name).split(">")[1:]
    data = {}

    for entry in entries:
        split_entry = entry.split("\n")
        data[split_entry[0]] = ''.join([line for line in split_entry[1:]])

    return data


