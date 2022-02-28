import pandas as pd
import BeautifulSoup
import requests


df = pd.read_excel("VBPO.xls")

print(df.code)

codeSeries = df.code
codeList = codeSeries.tolist()

for i in range(len(codeList)):
    url = "https://www.ncbi.nlm.nih.gov/protein/" + codeList[i] + "?report=fasta"
    print(url)
    r = requests.get(url, stream=True)
    soup = BeautifulSoup.BeautifulSoup(r.raw)
    sequence = soup.find("div", {"id": "articlebody"})
    print(sequence)
