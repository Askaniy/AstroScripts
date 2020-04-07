from math import sqrt
from PIL import Image, ImageDraw, ImageFilter
Image.MAX_IMAGE_PIXELS = None

# If you want to use this script not on Pluto or Charon height maps, you have to adapt the code a bit

open_path = "_/_.png" # to change
save_path = "_/_.png" # to change

body = "Charon" #"Pluto"
result_bit = 8 #16


if body == "Pluto":
    b_min = -4101
    b_max = 6491
elif body == "Charon":
    b_min = -14133
    b_max = 6911

if result_bit == 8:
    r_type = "L"
    r_mid = 127
    r_end = 255
elif result_bit == 16:
    r_type = "I"
    r_mid = 32767
    r_end = 65535

img0 = Image.open(open_path)
img1 = Image.new(r_type, img0.size, 0)
mask = Image.new(r_type, img0.size, 0)
back = Image.new(r_type, img0.size, r_mid)

if img0.getbands()[0] == "L":
    sourse_bit = 8
    s_type = "L"
    s_mid = 127
    s_end = 255
elif img0.getbands()[0] == "I":
    sourse_bit = 16
    s_type = "I"
    s_mid = 32767
    s_end = 65535

base = Image.new(s_type, (img0.size[0] * 2, img0.size[1]), 0)
base.paste(img0, (int(-img0.size[0] / 2), 0))
base.paste(img0, (int(img0.size[0] / 2), 0))
base.paste(img0, (int(3 * img0.size[0] / 2), 0))

pixb = base.load()
draw_surf = ImageDraw.Draw(img1)
draw_mask = ImageDraw.Draw(mask)

k = 1/42 # radius of effect relative to card length
n_max = int(img0.size[0] * k)

for x in range(base.size[0]):
    for y in range(base.size[1]):
        b0 = pixb[x, y]
        if b0 == 0:
            a = 1
            n = 0
            flag = True
            while flag:
                n += 1
                if x + n >= base.size[0] or y + n >= base.size[1] or x - n <= 0 or y - n <= 0:
                    break
                if n >= n_max:
                    break
                for i in range(1, n+1):
                    if pixb[x+i, y+n] != 0:
                        flag = False
                        b0 = pixb[x+i, y+n]
                    elif pixb[x+i, y-n] != 0:
                        flag = False
                        b0 = pixb[x+i, y-n]
                    elif pixb[x-i, y+n] != 0:
                        flag = False
                        b0 = pixb[x-i, y+n]
                    elif pixb[x-i, y-n] != 0:
                        flag = False
                        b0 = pixb[x-i, y-n]
                    elif pixb[x+n, y+i] != 0:
                        flag = False
                        b0 = pixb[x+n, y+i]
                    elif pixb[x+n, y-i] != 0:
                        flag = False
                        b0 = pixb[x+n, y-i]
                    elif pixb[x-n, y+i] != 0:
                        flag = False
                        b0 = pixb[x-n, y+i]
                    elif pixb[x-n, y-i] != 0:
                        flag = False
                        b0 = pixb[x-n, y-i]
                    if flag == False:
                        l = sqrt(n**2 + i**2)
                        a = (1 + l / n_max) / 2
                        break
        else:
            a = 0
            n = 0
            flag = True
            while flag:
                n += 1
                if x + n >= base.size[0] or y + n >= base.size[1] or x - n <= 0 or y - n <= 0:
                    break
                if n >= n_max:
                    break
                for i in range(1, n+1):
                    if 0 in [pixb[x+i, y+n], pixb[x+i, y-n], pixb[x-i, y+n], pixb[x-i, y-n], pixb[x+n, y+i], pixb[x+n, y-i], pixb[x-n, y+i], pixb[x-n, y-i]]:
                        flag = False
                        l = sqrt(n**2 + i**2)
                        a = (1 - l / n_max) / 2
        if result_bit == sourse_bit:
            b1 = b0
        else:
            b1 = int(r_end * (0.5 + (b0 - s_mid) / (2 * max(abs(b_min), b_max))))
        x1 = x - img0.size[0] / 2
        if 0 <= x1 < img0.size[0]:
            draw_surf.point((x1, y), b1)
            draw_mask.point((x1, y), int(a * r_end))
    p = x * 200 / img0.size[0]
    if p % 5 < 0.01:
        print(str(int(p / 4)) + " %")

img1.paste(back, (0, 0), mask)
img1_blur = img1.filter(ImageFilter.GaussianBlur(n_max))
img1.paste(img1_blur, (0, 0), mask)
img1.save(save_path)
print("100 %")