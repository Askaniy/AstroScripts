from PIL import Image, ImageDraw

save_path = "_/"                  # to change
save_name = "coordinate-grid-map" # to change
save_type = "png"


def check(d, h, x, y):
    k = h / d
    if x % k == 0 or y % k == 0 or x % k == k - 1 or y % k == k - 1:
        return True

q = int(input("Quality of the result (from 1 to 100): "))
if q < 1 or q > 100:
    q = 50
    print ("The limit is exceeded. The value is set to 50.")
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
img.save(save_path + save_name + "_quality-" + str(q) + "." + save_type)