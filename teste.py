import pandas as pd

with open("./output/arquivo1.txt") as arq:
    f= pd.DataFrame(data= arq ,columns=["TIPO_INTERAÇÃO, TIPO_ATOMO1, RESIDUO_AMINOACIDO1, CADEIA_AMINOACIDO1, TIPO_ATOMO2, RESIDUO_AMINOACIDO2, CADEIA_AMINOACIDO2, DISTANCIA"])

print(f.head())


