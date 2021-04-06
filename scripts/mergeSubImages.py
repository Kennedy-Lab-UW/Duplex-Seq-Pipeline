from argparse import ArgumentParser
import math
from PIL import Image


def main():
    parser = ArgumentParser()
    parser.add_argument("inBed"), 
    parser.add_argument("inPrefix")
    o = parser.parse_args()

    bedLinesPerImage = 4

    inBed = open(o.inBed, 'r')
    numLines = 0.
    for line in inBed:
        if line.strip() != "":
            numLines += 1
    numImages = math.ceil(numLines / bedLinesPerImage)
    firstImage = Image.open(f"{o.inPrefix}.1.targetCoverage.png")
    imageStart = 0
    my_images = []
    totalSize = 0
    for imgIter in range(numImages):
        my_images.append(Image.open(f"{o.inPrefix}.{imgIter + 1}.targetCoverage.png"))
    x_size = my_images[0].size[0]
    totalSize = sum([x.size[1] for x in my_images])
    outImage = Image.new('RGB', (x_size, totalSize), (250,250,250))
    for imgIter in range(numImages):
        nextImage = Image.open(f"{o.inPrefix}.{imgIter + 1}.targetCoverage.png")
        outImage.paste(nextImage, (0,imageStart))
        imageStart += nextImage.size[1]
    outImage.save(f"{o.inPrefix}.targetCoverage.png")

if __name__ == "__main__":
    main()
