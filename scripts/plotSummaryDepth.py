from PIL import Image
from argparse import ArgumentParser
import pandas as pd

parser = ArgumentParser()
parser.add_argument("config")
o = parser.parse_args()
samples = pd.read_csv(o.config).set_index("sample", drop=False)


myImages = []
for sampIter in samples.index:
    myImages.append(Image.open(f'{samples.loc[sampIter,"baseDir"]}/Stats/plots/{sampIter}.dcs.targetCoverage.png').convert('RGB'))

myImages[0].save(f'{o.config}.summaryDepth.pdf', save_all=True, append_images=myImages[1:])