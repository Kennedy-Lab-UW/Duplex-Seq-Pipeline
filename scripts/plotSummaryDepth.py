import sys
from argparse import ArgumentParser
import pandas as pd
from PIL import Image
import img2pdf

parser = ArgumentParser()
parser.add_argument("config")
o = parser.parse_args()
samples = pd.read_csv(o.config).set_index("sample", drop=False)
mySamples = [f'{samples.loc[x, "baseDir"]}/Stats/plots/{x}.dcs.targetCoverage.png' for x in samples.index]
with open(f'{o.config}.summaryDepth.pdf', 'wb') as f:
    f.write(img2pdf.convert(mySamples))
