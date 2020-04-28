from PIL import Image
from argparse import ArgumentParser, ArgumentTypeError

def isPNG(fileName):
    if fileName[-4:].upper() != '.PNG':
        msg = f"{fileName} is not a PNG file"
        raise ArgumentTypeError(msg)
    return fileName

def isTIFF(fileName):
    if (
            fileName[-5:].upper() != '.TIFF' 
            and fileName[-4:].upper() != '.TIF'
            ):
        msg = f"{fileName} is not a TIFF file"
        raise ArgumentTypeError(msg)
    return fileName

def main():
    parser = ArgumentParser()
    parser.add_argument('inTiff', type=isTIFF)
    parser.add_argument('outPng', type=isPNG)
    o = parser.parse_args()
    
    Image.open(o.inTiff).save(o.outPng)

if __name__ == "__main__":
    main()