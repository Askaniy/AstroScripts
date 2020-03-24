
# Average brightness calculator of a planetocentric spheroid texture
# The algorithm was proposed by Artyom Volgin (Zemlyanin) and was implemented in June 2019, when I began to study Python


from math import pi, sqrt, cos
from PIL import Image, ImageDraw

name = "_._" # file name
path = "_"   # file path

obl = 0.06487 # Oblateness (of Jupiter for example)
cum = 0 # Cumulative color conversion for the entire texture
avg = 0 # The average color value for the entire texture
n = 0 # The total number of texture pixels to be counted

img = Image.open(path + "/" + name)
w = img.size[0]
h = img.size[1]
pix = img.load()

for y in range(h - 1):
    lat = pi * (0.5 - (y + 0.5) / h)
    k = (1 - obl) * cos(lat) / sqrt(1 - (cos(lat)) ** 2 * (2 * obl - obl ** 2))
    c_y = 0 
    Ny = w
    for x in range(w - 1):
        b = pix[x, y]
        if b != 0 and b != 255:
            c_y += b
        else:
            Ny -= 1
    cum += c_y * k
    n += Ny * k
avg = cum / n

print(avg)