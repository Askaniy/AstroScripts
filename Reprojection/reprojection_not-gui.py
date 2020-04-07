from PIL import Image, ImageDraw
from math import ceil, floor, pi, cos, atan, tan
from scipy.interpolate import interp1d

obl = 0.06487  # of Jupiter, to change
open_path = "_/Reprojection/Jupiter test/jupiter-for-test.png" # to change
save_path = "_/Reprojection/Jupiter test/script-result.png"    # to change
mode = "pl-centric_pl-graphic" # "pl-graphic_pl-centric" "lambert_pl-centric"

img_old = Image.open(open_path)
pix = img_old.load()

w = img_old.size[0]
h_old = img_old.size[1]
if mode == "lambert_pl-centric":
    h_new = floor(h_old * pi / 2)
else:
    h_new = h_old
    
img_new = Image.new("RGB", (w, h_new), (0, 0, 0))
draw = ImageDraw.Draw(img_new)

if mode != "lambert_pl-centric":
    n = (1 - obl)**2
    hpi = h_old / pi
    pih = pi / h_old

for x in range(w):
    line_r = []
    line_g = []
    line_b = []
    for y in range(h_old):
        line_r.append(pix[x, y][0])
        line_g.append(pix[x, y][1])
        line_b.append(pix[x, y][2])
    interp_r = interp1d(range(h_old), line_r, kind = "cubic")
    interp_g = interp1d(range(h_old), line_g, kind = "cubic")
    interp_b = interp1d(range(h_old), line_b, kind = "cubic")
    for y in range(h_new - 1):
        if mode == "pl-centric_pl-graphic":
            y_old = hpi * (atan(tan(pih * (y + 0.5)) / n)) # planetoCENTRIC to planetoGRAPHIC
        elif mode == "pl-graphic_pl-centric":
            y_old = hpi * (atan(tan(pih * (y + 0.5)) * n)) # planetoGRAPHIC to planetoCENTRIC
        else:
            y_old = h_old / 2 * (1 - cos(2 * y / h_old))   # lambert to planetoCENTRIC
        if y_old < 0:
            y_old += h_old - 1
        if y_old > h_old - 1:
            y_old = h_old - 1
        r = int(interp_r(y_old))
        g = int(interp_g(y_old))
        b = int(interp_b(y_old))
        draw.point((x, y), (r, g, b))

img_new.save(save_path)