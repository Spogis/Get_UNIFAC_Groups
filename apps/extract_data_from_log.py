import pandas as pd


def create_dataset_from_log(log_path):
    with open(log_path, 'r') as file:
        log_content = file.read()

    # Identificar as seções de início e fim para a extração
    start_marker = "Fragmentation from the algorithm:"
    end_marker = "Fragmentation from the reference database:"

    # Extrair a seção relevante do log
    start_index = log_content.find(start_marker) + len(start_marker)
    end_index = log_content.find(end_marker)
    extracted_text = log_content[start_index:end_index].strip()

    # Processar o texto extraído para criar o dataset
    data_lines = extracted_text.strip().split("\n")
    data_rows = [line.split(maxsplit=2) for line in data_lines]

    # Criar o DataFrame
    df = pd.DataFrame(data_rows, columns=['Subgroup Name', 'Subgroup Number', 'Count'])
    df = df[['Subgroup Number', 'Subgroup Name', 'Count']]


    return df
