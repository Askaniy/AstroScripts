from PIL import Image, ImageDraw
from math import ceil, floor, pi, cos, atan, tan
from scipy.interpolate import interp1d
import PySimpleGUI as sg

layout = [
    [sg.Text("Import file:")], [sg.Input(), sg.FileBrowse()],
    [sg.Text("Export folder:")], [sg.Input(), sg.FolderBrowse()],
    [sg.Text("Export file name:")], [sg.InputText("script-result")],
    [sg.Text("Choose projection mode:")], 
    [sg.InputCombo(("Planetocentric to planetographic", "Planetographic to planetocentric", "Lambert to planetocentric"), size=(30, 1))],
    [sg.Text("Oblateness (if required, by default - Jupiter's):")], [sg.InputText("0.06487")],
    [sg.OK(), sg.Cancel()],
    [sg.Text("Progress:")], [sg.ProgressBar(1000, orientation = "h", size = (30, 20), key = "progress")]
]

window = sg.Window("Map reprojection", layout)
event, values = window.Read()
mode = values[3]
obl = float(values[4])
img_old = Image.open(values[0])
pix = img_old.load()

w = img_old.size[0]
h_old = img_old.size[1]
if mode == "Lambert to planetocentric":
    h_new = floor(h_old * pi / 2)
else:
    h_new = h_old
    
img_new = Image.new("RGB", (w, h_new), (0, 0, 0))
draw = ImageDraw.Draw(img_new)

if mode != "Lambert to planetocentric":
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
        if mode == "Planetocentric to planetographic":
            y_old = hpi * (atan(tan(pih * (y + 0.5)) / n)) # planetoCENTRIC to planetoGRAPHIC
        elif mode == "Planetographic to planetocentric":
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
    window["progress"].update_bar(x / w * 1000)

window.Close()
img_new.save(values[1] + "/" + values[2] + ".png")

# for EXE print: pyinstaller -wF "x:/Documents/GitHub/AstroScripts/Map calculations/Reprojection/reprojection_gui.py"