import dash
from dash import dcc, html, dash_table, Input, Output, State
from dash.dash_table.Format import Format, Scheme
import dash_bootstrap_components as dbc
import pandas as pd
import pubchempy as pcp

from apps.get_unifac_groups import *

# Inicializa o app Dash
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.title = "UNIFAC Groups"

# Definir o layout do aplicativo
app.layout = html.Div([
    html.Br(),
    html.Br(),
    html.Div([
        html.Img(src='assets/logo.png', style={'height': '100px', 'margin-left': 'auto', 'margin-right': 'auto'}),
    ], style={'text-align': 'center', 'margin-bottom': '10px'}),
    html.Br(),
    html.Div([
        html.A("Find data on PubChem", href="https://pubchem.ncbi.nlm.nih.gov/", target="_blank",
               style={'fontSize': 20}),
    ], style={'textAlign': 'center', 'marginBottom': '20px'}),

    html.Div([
        html.Div([
            html.Label("PubChem CID:", style={'fontSize': 20}),
            dcc.Input(id='input-cid', type='text', placeholder='PubChem CID',
                      style={'width': '100%', 'padding': '10px', 'fontSize': 18},),
        ], style={'marginBottom': '20px'}),

        html.Div([
            html.Label("IUPAC Name:", style={'fontSize': 20}),
            dcc.Textarea(id='IUPAC-Name',
                         style={'width': '100%', 'padding': '10px', 'fontSize': 18,
                                'resize': 'none' },
                         disabled=True),
        ], style={'marginBottom': '20px'}),

        html.Div([
            html.Label("Molecular Formula:", style={'fontSize': 20}),
            dcc.Textarea(id='Molecular-Formula',
                         style={'width': '100%', 'padding': '10px', 'fontSize': 18,
                                'resize': 'none' },
                         disabled=True),
        ], style={'marginBottom': '20px'}),

        html.Button('Determine UNIFAC Groups', id='submit-button', n_clicks=0, style={'fontSize': 20, 'padding': '10px 20px'}),
    ], style={'padding': '20px', 'margin': 'auto', 'width': '50%', 'boxShadow': '0px 0px 10px #ccc',
              'borderRadius': '15px'}),

    html.Br(),
    html.Div([
        html.Div(id='table-container', style={'margin': '0 auto'})
    ], style={'padding': '20px', 'margin': 'auto', 'width': '50%', 'boxShadow': '0px 0px 10px #ccc',
              'borderRadius': '15px'}),

], style={'fontFamily': 'Arial, sans-serif'})


# Callback para alternar a visibilidade da tabela
@app.callback(Output('table-container', 'children'),
              Output('IUPAC-Name', 'value'),
              Output('Molecular-Formula', 'value'),
              Input('submit-button', 'n_clicks'),
              State('input-cid', 'value'))
def display_table(n_clicks, cid):
    if n_clicks > 0 and all([cid]):
        c = pcp.Compound.from_cid(cid)
        inchikey = c.inchikey
        smiles = c.canonical_smiles
        molecular_formula = c.molecular_formula
        iupac_name = c.iupac_name

        # Definir o texto que será gravado no arquivo CSV
        csv_content = inchikey + "," + smiles + "," + cid + ",1:0|2:0|3:0"

        # Definir o caminho do arquivo CSV onde o conteúdo será gravado
        csv_file_path = 'apps/reference_DB.csv'  # Ajuste o caminho conforme necessário

        # Gravar o texto no arquivo CSV
        with open(csv_file_path, 'w') as file:
            file.write(csv_content)

        print(f"Conteúdo gravado com sucesso em {csv_file_path}")

        df = Get_UNIFAC_Groups()
        Table = dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            data=df.to_dict('records'),
            style_cell={'textAlign': 'left', 'padding': '10px', 'fontSize': '18px'},
            style_header={
                'backgroundColor': 'lightgrey',
                'fontWeight': 'bold',
                'textAlign': 'center',
                'fontSize': '20px'
            },
            style_data_conditional=[
                {
                    'if': {'row_index': 'odd'},
                    'backgroundColor': 'rgb(248, 248, 248)'
                }
            ],
            style_table={'margin': 'auto'},  # Assegura que a tabela será centralizada
            fill_width=False  # Impede a tabela de automaticamente preencher a largura
        )
        return Table, iupac_name, molecular_formula

    else:  # Caso contrário, não mostrar nada
        return "", "", ""

# Rodar o aplicativo
if __name__ == '__main__':
    app.run_server(debug=False)
