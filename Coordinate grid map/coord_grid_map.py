from PIL import Image, ImageDraw

def check(d, h, x, y):
    k = h / d
    if x % k == 0 or y % k == 0 or x % k == k - 1 or y % k == k - 1:
        return True

name = "coord_grid_map"
path = "_"
form = "_"

q = int(input("The quality of the final map (from 1 to 100): "))
if q < 1 or q > 100:
    q = 50
    print ("Error. The value is set to 50.")
h = 18 * q
w = 36 * q

img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
draw = ImageDraw.Draw(img)

for x in range(w):
    for y in range(h):
        a = 255
        g = 255
        b = 255
        if check(2, h, x, y):
            True
        elif check(6, h, x, y):
            b = 0
        elif check(18, h, x, y):
            g = 0
            b = 0
        else:
            a = 0
        draw.point((x, y), (255, g, b, a))

print("Done!")
img.save(path + "/" + name + "-script_result-q_" + str(q) + "." + form)