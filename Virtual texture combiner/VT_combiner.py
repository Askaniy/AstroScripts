from PIL import Image, ImageDraw
Image.MAX_IMAGE_PIXELS = None

open_path = "_/"      # to change
save_path = "_/_.png" # to change

level = 6             # to change

if level == 0:
    l = 1024
else:
    l = 512
n = 2**(level + 1)

result = Image.new("RGB", (int(n * l), int(n/2 * l)))
for x in range(n):
    for y in range(int(n/2)):
        img = Image.open(open_path + "tx_{}_{}.png".format(x, y))
        result.paste(img, (x * l, y * l))

result.save(save_path)