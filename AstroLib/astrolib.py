
# AstroLib, version of June 28, 2020


from scipy import pi, sin, cos, arcsin, arccos, sqrt
import astropy.constants as c
import astropy.units as u
import codecs
import time
import re
import os
import sys
import discord
import json


bot = False
database_path = ".../Databases"          # to change
#discord_path = ".../Discord data"       # to change
token_path = ".../discord_bot_token.txt" # to change


supported_requests = {
    "help": "You really need help to request help?", 
    "create": """This module allows you to create your own object that can be processed as object from database. Input example:
    `create Russell's Teapot, teapot
    radius = 10 cm
    mass = 1 kg`""", 
    "return": """This module searches for the object in the available databases and return data without any processing or formalizing. Input example:
    return Jupiter X""", 
    "search": """This module searches for the object in the available databases and return formalized data. Input example:
    `search Jupiter X`""",
    "find": """This module searches for parameter of specified celestial body in the available databases. 
Specify the parameter before “of” and the celestial body after, or first a planet, and then a parameter, if each consists of one word. 
If you need an epoch, write after “on”. Input examples:
    `find surfase area of Jupiter X in m2`
    `find Jovian density in cgs`""",
    "generate": """This module generates an encoding in the specified format. Specify the format before “for” and the celestial body after.
If you need an epoch, write after “on”. Input example:
    `generate SSC for Jupiter X`"""}

formats = ["ssc", "json", "decor"]

links = {
    "WGCCRE2015": "https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf",
    "IAU2015 Resolution B3": "http://www.astro.osu.edu/~pogge/Ast2292/Docs/NominalConstants.pdf",
    "JPL Planetary Satellite Physical Parameters": "https://ssd.jpl.nasa.gov/?sat_phys_par",
    "johnstonsarchive.net Satellite data": "http://www.johnstonsarchive.net/astro/sssatellitedata.html"
}

def surface(a, b, c):
    if a == b == c:
        return 4 * pi * a**2
    elif a == b:
        e = sqrt(1 - c / a)
        return 2 * pi * a**2 * (1 + a * arcsin(e) / (c * e))
    else:
        p = 1.6075
        return 4 * pi * ((a**p * b**p + a**p * c**p + b**p * c**p) / 3)**(1 / p)

prmtr = {
    "mean radius": {
        "alt_names": set(["radius", "mean radius"]), 
        "formulas": [
            (["mean diameter"], lambda d: d / 2), 
            (["mean equatorial radius", "mean polar radius"], lambda a, c: (2 * a + c) / 3)
            ]},
    "mean diameter": {
        "alt_names": set(["diameter", "mean diameter"]),
        "formulas": [(["mean radius"], lambda r: r * 2)]},
    "mean equatorial radius": {
        "alt_names": set(["equatorial radius", "mean equatorial radius"]),
        "formulas": [
            (["mean equatorial diameter"], lambda d: d / 2), 
            (["long equatorial radius", "short equatorial radius"], lambda a, b: (a + b) / 2), 
            (["mean polar radius", "mean oblateness"], lambda c, obl: c / (1 - obl))
            ],
        "included": "mean radius"},
    "mean equatorial diameter": {
        "alt_names": set(["equatorial diameter", "mean equatorial diameter"]),
        "formulas": [(["mean equatorial radius"], lambda r: r * 2)],
        "included": "mean diameter"},
    "long equatorial radius": {
        "alt_names": set([
            "x radius", "x axis", "long radius", "long equatorial radius", 
            "subplanetary radius", "subplanetary equatorial radius"]),
        "formulas": [
            (["long equatorial diameter"], lambda d: d / 2),
            (["short equatorial radius", "equatorial oblateness"], lambda b, obl: b / (1 - obl))
            ],
        "included": "mean equatorial radius"},
    "long equatorial diameter": {
        "alt_names": set([
            "x diameter", "long diameter", "long equatorial diameter", 
            "subplanetary diameter", "subplanetary equatorial diameter"]),
        "formulas": [(["long equatorial radius"], lambda r: r * 2)],
        "included": "mean equatorial diameter"},
    "short equatorial radius": {
        "alt_names": set([
            "y radius", "y axis", "short radius", "short equatorial radius", 
            "along orbit equatorial radius", "along orbit radius"]),
        "formulas": [
            (["short equatorial diameter"], lambda d: d / 2),
            (["long equatorial radius", "equatorial oblateness"], lambda a, obl: a * (1 - obl))
            ],
        "included": "mean equatorial radius"},
    "short equatorial diameter": {
        "alt_names": set([
            "y diameter", "short diameter", "short equatorial diameter", 
            "along orbit equatorial diameter", "along orbit diameter"]),
        "formulas": [(["short equatorial radius"], lambda r: r * 2)],
        "included": "mean equatorial diameter"},
    "mean polar radius": {
        "alt_names": set(["z radius", "z axis", "polar radius", "mean polar radius"]),
        "formulas": [
            (["mean polar diameter"], lambda d: d / 2), 
            (["north polar radius", "south polar radius"], lambda n, s: (n + s) / 2),
            (["mean equatorial radius", "mean oblateness"], lambda a, obl: a * (1 - obl))
            ],
        "included": "mean radius"},
    "mean polar diameter": {
        "alt_names": set(["z diameter", "polar diameter", "mean polar diameter"]),
        "formulas": [(["mean polar radius"], lambda r: r * 2)],
        "included": "mean diameter"},
    "north polar radius": {
        "alt_names": set(["north radius", "north polar radius"]),
        "formulas": [(["mean polar diameter", "south polar radius"], lambda d, s: d - s)],
        "included": "mean polar radius"},
    "south polar radius": {
        "alt_names": set(["south radius", "south polar radius"]),
        "formulas": [(["mean polar diameter", "north polar radius"], lambda d, n: d - n)],
        "included": "mean polar radius"},
    "mean oblateness": {
        "alt_names": set(["oblateness", "mean oblateness"]),
        "formulas": [(["mean equatorial radius", "mean polar radius"], lambda ab, c: 1 - c / ab)]},
    "equatorial oblateness": {
        "alt_names": set(["equatorial oblateness"]),
        "formulas": [(["long equatorial radius", "short equatorial radius"], lambda a, b: 1 - b / a)]},
    "surface area": {
        "alt_names": set(["surface", "area", "surface area"]),
        "formulas": [(["long equatorial radius", "short equatorial radius", "mean polar radius"], surface)]},
    "volume": {
        "alt_names": set(["volume"]),
        "formulas": [
            (["long equatorial radius", "short equatorial radius", "mean polar radius"], lambda a, b, c: 4*pi * a*b*c /3),
            (["mass", "mean density"], lambda m, d: m / d)
            ]},
    "mean density": {
        "alt_names": set(["density", "mean density"]),
        "formulas": [(["mass", "volume"], lambda m, v: m / v)]},
    "mass": {
        "alt_names": set(["mass"]),
        "formulas": [
            (["standard gravitational parameter"], lambda gm: gm / c.G),
            (["volume", "mean density"], lambda v, d: v * d)
            ]},
    "standard gravitational parameter": {
        "alt_names": set(["mass parameter", "gravitational parameter", "standard gravitational parameter"]),
        "formulas": [(["mass"], lambda m: c.G * m)]},
    "gravitational acceleration": {
        "alt_names": set(["gravitational acceleration"]),
        "formulas": [(["standard gravitational parameter", "mean radius"], lambda gm, r: gm / r**2)]},
    "minimum gravitational acceleration": {
        "alt_names": set(["minimum gravitational acceleration"]),
        "formulas": [(["standard gravitational parameter", "long equatorial radius", "rotation period"], lambda gm, r, p: gm / r**2 - 4 * pi**2 * r / p**2)]},
    "maximum gravitational acceleration": {
        "alt_names": set(["maximum gravitational acceleration"]),
        "formulas": [(["standard gravitational parameter", "mean polar radius"], lambda gm, r: gm / r**2)]},
    "north pole longitude": {
        "alt_names": set(["north pole longitude"])},
    #    "formulas": [([ , ], lambda  , )]},
    "north pole latitude": {
        "alt_names": set(["north pole latitude"])},
    #    "formulas": [([ , ], lambda  , )]},
    "prime meridian angle": {
        "alt_names": set(["meridian angle", "prime meridian angle"])},
    #    "formulas": [([ , ], lambda  , )]},
    "rotation period": {
        "alt_names": set(["rotation period"])},
    #    "formulas": [([ , ], lambda  , )]},
    "spin–orbit resonance": {
        "alt_names": set(["spin–orbit resonance", "tidal locking", "gravitational locking", "spin-orbit locking"])},
    #    "formulas": [([ , ], lambda  , )]},
    "orbital period": {
        "alt_names": set(["orbital period"]),
        "horizons": "P"},
    #    "formulas": [([ , ], lambda  , )]},
    "semimajor axis": {
        "alt_names": set(["semimajor axis"]),
        "horizons": "a"},
    #    "formulas": [([ , ], lambda  , )]},
    "eccentricity": {
        "alt_names": set(["eccentricity"]),
        "horizons": "e"},
    #    "formulas": [([ , ], lambda  , )]},
    "inclination": {
        "alt_names": set(["inclination"]),
        "horizons": "incl"},
    #    "formulas": [([ , ], lambda  , )]},
    "longitude of ascending node": {
        "alt_names": set(["longitude of ascending node"])},
    #    "formulas": [([ , ], lambda  , )]},
    "argument of periapsis": {
        "alt_names": set(["argument of periapsis", "argument of pericenter"])},
    #    "formulas": [([ , ], lambda  , )]},
    "mean anomaly": {
        "alt_names": set(["mean anomaly"])},
    #    "formulas": [([ , ], lambda  , )]},
}


# Replacement block:

abbrevs = {
    "a": "semimajor axis",
    "e": "eccentricity",
    "i": "inclination",
    "r": "radius",
    "d": "diameter",
    "n": "north",
    "s": "south",
    "p": "pole",
    "w": "prime meridian angle",
    "m": "mass",
    "g": "gravitational acceleration",
    "gm": "standard gravitational parameter",
    "ly": "light years",
    "ec": "ecliptic",
    "eq": "equatorial",
    "pl": "polar",
    "ra": "right ascension",
    "dec": "declination",
    "obl": "oblateness",
    "max": "maximum",
    "min": "minimum",
    "lat": "latitude",
    "lon": "longitude",
    "grav": "gravitational",
    "tons": "tonnes",
    "radii": "radius",
    "masses": "mass"
}

def full(name):
    if name in abbrevs:
        return abbrevs[name]
    return name

adjectives = {
    "solar": "sun",
    "mercurian": "mercury",
    "venerian": "venus",
    "earthly": "earth",
    "lunar": "moon",
    "martian": "mars",
    "phobian": "phobos",
    "deimosian": "deimos",
    "cererian": "ceres",
    "jovian": "jupiter",
    "ionian": "io",
    "europan": "europa",
    "ganymedean": "ganymede",
    "callistoan": "callisto",
    "saturnian": "saturn",
    "titanian": "titan",
    "uranian": "uranus",
    "neptunian": "neptune",
    "tritonian": "triton",
    "plutonian": "pluto"
}

def noun(name):
    if name in adjectives:
        return adjectives[name]
    return name

units = {
    "seconds": u.s,
    "minutes": u.min,
    "hours": u.h,
    "days": u.d,
    "weeks": u.wk,
    "years": u.yr,
    "nanometers": u.nm,
    "millimeters": u.mm,
    "centimeters": u.cm,
    "meters": u.m,
    "kilometers": u.km,
    "grams": u.g,
    "kilograms": u.kg,
    "tonnes": u.t,
    "astronomical units": u.AU,
    "light years": u.lyr,
    "sun radius": u.solRad,
    "jupiter radius": u.jupiterRad,
    "earth radius": u.earthRad,
    "sun mass": u.solMass,
    "jupiter mass": u.jupiterMass,
    "earth mass": u.earthMass
}

def unit(name):
    if name in units:
        return units[name]
    return u.Unit(name)


# Deep layers of the script:

def check_request(request):
    request = [x.lower() for x in request.split(" ")]
    if request[0] in supported_requests:
        return read_request(request)

def read_request(request):
    meaning = {}
    if "the" in request:
        request.remove("the")
    if request[0] == "help":
        if len(request) == 1:
            meaning.update({"to_do": "help", "meta": "info"})
        else:
            if "with" in request:
                request.remove("with")
            if "in" in request:
                request.remove("in")
            meaning.update({"to_do": "help", "meta": request[1:]})
        return meaning
    elif request[0] == "create":
        meaning.update({"to_do": "create", "body": " ".join(request[1:])})
        return meaning
    elif request[0] == "return":
        meaning.update({"to_do": "return", "body": " ".join(list(map(noun, request[1:])))})
        return meaning
    elif request[0] == "search":
        meaning.update({"to_do": "search", "body": " ".join(list(map(noun, request[1:])))})
        return meaning
    elif request[0] == "find":
        request.remove("find")
        meaning.update({"to_do": "find"})
        if "on" in request:
            x = request.index("on")
            meaning.update({"epoch": request[x+1]})
            request = request[:x]
        if request[-1] in ["value", "error", "unit", "source", "comment", "info", "data"]:
            meaning.update({"meta": request[-1]})
            request = request[:-1]
        else:
            meaning.update({"meta": "info"})
        if "in" in request:
            x = request.index("in")
            meaning.update({"unit": " ".join(list(map(lambda i: noun(full(i)), request[x+1:])))})
            request = request[:x]
        if "of" in request:
            x = request.index("of")
            meaning.update({"body": " ".join(request[x+1:])})
            meaning.update({"parameter": " ".join(list(map(full, request[:x])))})
            return meaning
        elif len(request) == 2:
            meaning.update()
            meaning.update({"body": noun(request[0]), "parameter": full(request[1])})
            return meaning
        else:
            print("Error. Try using “of” construct.")
    elif request[0] == "generate":
        request.remove("generate")
        meaning.update({"to_do": "generate"})
        if "on" in request:
            x = request.index("on")
            meaning.update({"epoch": request[x+1]})
            request = request[:x]
        if "for" in request:
            meaning.update({"body": " ".join(request[2:]), "format": request[0]})
        else:
            meaning.update({"body": " ".join(request[:-1]), "format": request[-1]})
        if meaning["format"] in formats:
            return meaning
    #print("- result of read_request: " + str(meaning))

def quantity(request):
    spl = request.split(" ")
    try:
        value = float(spl[0])
    except ValueError:
        return request
    else:
        return value * unit(" ".join(spl[1:]))

def prmt_reader(line):
    if ":" in line:
        parameter = line.split(":")
        name = parameter[0].strip()
        value = parameter[1].strip()
    elif "=" in line:
        parameter = line.split("=")
        name = parameter[0].strip()
        value = quantity(parameter[1].strip())
    else:
        name = line.strip()
        value = "use `:` for text and `=` for quantity"
    return {"name": name, "value": value}

def find_body(request, database_path):
    global cannot_calc
    cannot_calc = []
    body = {}
    files = [f.path for f in os.scandir(database_path) if f.is_file()]
    files.extend([f.path for f in os.scandir(discord_path) if f.is_file()])
    files.reverse()
    print(files)
    for fl in files:
        if os.path.splitext(fl)[1] == ".askaniy":
            with codecs.open(fl) as f:
                f_list = list(f)
                f_list.extend(["\n", "Endgame"])
                obj_level = 0
                part_of = []
                last_parameter = ""
                check_system_flag = False
                system = False
                read_obj_name = True
                read_obj_prmt = False
                read_group = False
                read_head = False
                head = {}
                for line in f_list:
                    if not system:
                        if not read_obj_prmt:
                            if read_obj_name and not line.isspace():
                                obj_level = line.count("\t")
                                names = line.strip().split(", ")
                                noun_version = list(map(lambda i: " ".join(list(map(lambda j: noun(j.lower()), i.split(" ")))), names))
                                if request in noun_version:
                                    read_obj_prmt = True
                                    body["name"] = names
                                else:
                                    part_of.insert(obj_level, names[0])
                            read_obj_name = False
                            if line.isspace():
                                read_obj_name = True
                        else:
                            if check_system_flag:
                                if not line.isspace():
                                    if line.count("\t") > obj_level:
                                        system = True
                                        includes = [line.strip()]
                                    else:
                                        return body
                            else: 
                                if line.isspace():
                                    body.update({"part of": part_of[:obj_level]})
                                    check_system_flag = True
                                else:
                                    # Parameters reading
                                    if not read_group:
                                        if "<" in line:
                                            read_group = True
                                            read_head = True
                                            head.update(prmt_reader(line.replace("<", "")))
                                        else:
                                            parameter = prmt_reader(line) # < parameter reading block
                                            if line.count("\t") == obj_level:
                                                body[parameter["name"]] = {"value": parameter["value"]}
                                                last_parameter = parameter["name"]
                                            else:
                                                body[last_parameter].update({parameter["name"]: parameter["value"]}) # >
                                    else:
                                        if read_head and line.count("\t") == obj_level + 1:
                                            head.update(prmt_reader(line))
                                        else:
                                            read_head = False
                                            if ">" in line:
                                                read_group = False
                                                line = line.replace(">", "")
                                            parameter = prmt_reader(line) # < parameter reading block
                                            if line.count("\t") == obj_level:
                                                body[parameter["name"]] = {"value": parameter["value"]}
                                                body[parameter["name"]].update({head["name"]: head["value"]})
                                                last_parameter = parameter["name"]
                                            else:
                                                body[last_parameter].update({parameter["name"]: parameter["value"]}) # >
                    else:
                        if read_obj_name:
                            if line.count("\t") == obj_level + 1:
                                includes.append(line.strip())
                            elif line.count("\t") <= obj_level:
                                body.update({"includes": includes})
                                return body
                        read_obj_name = False
                        if line.isspace():
                            read_obj_name = True


calculating = []
cannot_calc = []

def calculate(name, body):
    global calculating
    global cannot_calc
    print("- request to calculate: " + name)
    calculating.append(name)
    if "formulas" in prmtr[name]:
        for formula in prmtr[name]["formulas"]:
            need_parameter_values = []
            for need_parameter in formula[0]:
                if need_parameter not in calculating and need_parameter not in cannot_calc:
                    need_parameter_values.append(process(need_parameter, body))
                else:
                    need_parameter_values = []
                    break
            print("- need_parameter_values: " + str(need_parameter_values))
            if need_parameter_values != [] and None not in need_parameter_values:
                to_formula_value = []
                result = {"database name": name, "comment": "*calculated*"}
                for need_parameter_value in need_parameter_values:
                    for component in body[need_parameter_value["database name"]]:
                        if component not in ["value", "error", "database name"]:
                            if component in result:
                                if type(result[component]) != list:
                                    result[component] = [result[component]]
                                if type(body[need_parameter_value["database name"]][component]) != list:
                                    result[component].append(body[need_parameter_value["database name"]][component])
                                else:
                                    result[component].extend(body[need_parameter_value["database name"]][component])
                                result.update({component: list(set(result[component]))})
                            else:
                                result.update({component: body[need_parameter_value["database name"]][component]})
                    value = body[need_parameter_value["database name"]]["value"]
                    to_formula_value.append(value)
                try:
                    print("- to_formula_value: " + str(to_formula_value))
                    result.update({"value": formula[1](*to_formula_value).decompose()})
                    body.update({name: result})
                    print("- result of calculate (if was calculated): " + str(result))
                    return result
                except IndexError:
                    pass
                finally:
                    calculating.remove(name)
    calculating.remove(name)
    cannot_calc.append(name)
    print(name.capitalize() + " is unknown and cannot be calculated.")

def process(name, body):
    try:
        print("- request to process: " + name)
        cross = list(set(body) & prmtr[name]["alt_names"])[0]
    except IndexError:
        result = calculate(name, body)
        if result:
            print("- result of process (if was calculated): " + str(result))
            return result
        elif "included" in prmtr[name]:
            #result.update({"comment": ["*teoretical*"]})
            #body.update({name: result})
            return process(prmtr[name]["included"], body)
        #elif "standart" in prmtr[name]:
    else:
        body[cross].update({"database name": cross})
        print("- result of process (if was found): " + str(body[cross]))
        return body[cross]

def find_class(body):
    if "includes" in body:
        return "barycenter"
    else:
        return "work in progress"

def find_parameter(request, body):
    if request in ["part of", "parent", "parents"]:
        return {"database name": "parents", "value": ", ".join(body["part of"])}
    elif request in ["system", "include", "includes"]:
        return {"database name": "system", "value": ", ".join(body["includes"])}
    elif request == "class":
        return find_class(body)
    else:
        for name in prmtr:
            if request in prmtr[name]["alt_names"]:
                return process(name, body)


# Second layer of the script:

def hlp(request):
    global em
    if request == "info":
        anwsr = "Here is what I support:"
        if bot:
            em = discord.Embed(color = discord.Colour.from_rgb(127, 127, 255))
            em.add_field(name = "Requests:", value = "help " + "\nhelp ".join(list(supported_requests)[1:]), inline = True)
            em.add_field(name = "Lists:", value = "help " + "\nhelp ".join(["parameters", "adjectives", "abbreviations", "formats"]), inline = True)
        else:
            anwsr += "\nRequests:\n* help " + "\n* help ".join(list(supported_requests)[1:])
            anwsr += "\nLists:\n* help " + "\n* help ".join(["parameters", "adjectives", "abbreviations", "formats"])
        return anwsr
    elif request[0] in supported_requests:
        return supported_requests[request[0]]
    elif request[0] == "parameters":
        if bot:
            em = discord.Embed(title = "Supported parameters:", description = "\n".join(prmtr), color = discord.Colour.from_rgb(127, 127, 255))
        else:
            return "Supported parameters:\n* " + "\n* ".join(prmtr)
    elif request[0] == "adjectives":
        if bot:
            em = discord.Embed(title = "Supported adjectives:", color = discord.Colour.from_rgb(127, 127, 255))
            for key, value in adjectives.items():
                em.add_field(name = value, value = key, inline = True)
        else:
            r = ""
            for key, value in adjectives.items():
                r += "\n{} - {}".format(value, key)
            return "Supported adjectives:" + r
    elif request[0] == "abbreviations":
        if bot:
            em = discord.Embed(title = "Supported abbreviations and replacements:", color = discord.Colour.from_rgb(127, 127, 255))
            for key, value in abbrevs.items():
                em.add_field(name = value, value = key, inline = True)
        else:
            r = ""
            for key, value in abbrevs.items():
                r += "\n{} - {}".format(key, value)
            return "Supported abbreviations and replacements:" + r
    elif request[0] == "formats":
        if bot:
            em = discord.Embed(title = "Supported formats:", description = "\n".join(formats), color = discord.Colour.from_rgb(127, 127, 255))
        else:
            return "Supported formats:\n* " + "\n* ".join(formats)
    else:
        return "I can't help with " + " ".join(request)

def create(body):
    file_name = time.strftime("%Y-%m-%d_%H-%M-%S")
    with codecs.open("{}/{}.askaniy".format(discord_path, file_name), "w") as f:
        f.write(body)
    body_name = body.split("\n")[0]
    print(body_name.split(", ")[0])
    try:
        print(find_body(body_name.split(", ")[0], database_path))
        return "I saved your object successfully, you can work with it"
    except Exception:
        return "Some error in your syntax, I can't read it"

def output(parameter):
    try:
        anwser = str(parameter["value"].value)
        #anwser += str(round(parameter["value"].value, 5))
        if "error" in parameter:
            anwser += " ± " + str(parameter["error"])
        anwser += " " + str(parameter["value"].unit)
    except AttributeError:
        anwser = parameter["value"]
    if "source" in parameter:
        if type(parameter["source"]) == list:
            if len(parameter["source"]) == 1:
                anwser += "; source: " + str(parameter["source"][0])
            else:
                anwser += "; sources: " + ", ".join(parameter["source"])
        else:
            anwser += "; source: " + str(parameter["source"])
    if "comment" in parameter:
        if type(parameter["comment"]) == list:
            anwser += " ({})".format(", ".join(parameter["comment"]))
        else:
            anwser += " ({})".format(parameter["comment"])
    return anwser

def short_list(system):
    res = []
    for i in system:
        res.append(i.split(", ")[0])
    if len(res) >= 20:
        res = res[:19]
        res.append("...")
    return res

def embed(body):
    global em
    em = discord.Embed(
        title = ", ".join(body["name"]),
        description = "class: " + find_class(body),
        color = discord.Colour.from_rgb(127, 255, 127)
        )
    if body["part of"] != []:
        em.add_field(name = "part of: ", value = ", ".join(body["part of"]), inline = False)
    if "includes" in body:
        em.add_field(name = "includes: ", value = ", ".join(short_list(body["includes"])), inline = False)
    for i in body:
        if type(body[i]) != list:
            em.add_field(name = i, value = output(body[i]), inline = False)
    return "I found this:"

def find(request, body):
    parameter = find_parameter(request["parameter"], body)
    if parameter == None:
        return "I can't understand this parameter"
    else:
        if "unit" in request:
            if request["unit"] == "si":
                parameter.update({"value": parameter["value"].si})
            elif request["unit"] == "cgs":
                parameter.update({"value": parameter["value"].cgs})
            else:
                parameter.update({"value": parameter["value"].to(unit(request["unit"]))})
            if "error" in parameter:
                parameter.pop("error")
        if request["meta"] == "info":
            anwser = "{} of {} is ".format(parameter["database name"].capitalize(), body["name"][0].capitalize())
            return anwser + output(parameter)
        elif request["meta"] == "value":
            return parameter["value"]
        elif request["meta"] == "error" and "error" in parameter:
            return parameter["error"]
        elif request["meta"] == "unit" and type(parameter["value"]) == u.quantity.Quantity:
            return parameter["value"].unit
        elif request["meta"] == "source" and "source" in parameter:
            return parameter["source"]
        elif request["meta"] == "comment" and "comment" in parameter:
            return parameter["comment"]
        else:
            return request["meta"].title() + " is not in the requested parameter"

def generate(request, body):
    if request == "ssc":
        code = '"{}" "{}"\n'.format(":".join(body["name"]), "/".join(body["part of"])) + "{"
        code += "\n\t"
        return code
    elif request == "json":
        for i in body:
            try: 
                body[i].update({"value": str(body[i]["value"])})
            except Exception:
                pass
        return json.dumps(body)
    else:
        decor = "\n{}\n{}\n".format(", ".join(body["name"]), "class: " + find_class(body))
        if body["part of"] != []:
            decor += "\npart of:\n{}\n".format(", ".join(body["part of"]))
        if "includes" in body:
            decor += "\nincludes:\n{}\n".format(", ".join(short_list(body["includes"])))
        for key, value in body.items():
            if type(value) != list:
                decor += ("\n{}\n{}\n".format(key, output(value)))
        return decor


# First layer of the script:

def anwser(request):
    print(request)
    global em
    em = None
    if request:
        if request["to_do"] == "help":
            return hlp(request["meta"])
        elif request["to_do"] == "create":
            return create(request["body"])
        body = find_body(request["body"], database_path)
        if body == None:
            return "I can't find {} in the available databases".format(request["body"])
        elif request["to_do"] == "return":
            return "`" + str(body) + "`"
        elif request["to_do"] == "search":
            if bot:
                return embed(body)
            else:
                return generate("decor", body)
        elif request["to_do"] == "find":
            return find(request, body)
        elif request["to_do"] == "generate":
            return generate(request["format"], body)
    else:
        return "I can't do it for now"


# Surface of the script:

if bot:
    class BotClient(discord.Client):
        async def on_ready(self):
            print("Logged in as " + str(self.user))
        async def on_message(self, message):
            request = [x.lower() for x in message.content.split(" ")]
            if request[0] in supported_requests and message.author.discriminator != "1451":
                await message.channel.send(anwser(read_request(request)), embed = em)
    BotClient().run(open(token_path).read())
else:
    request = check_request(input("Write your request here: "))
    print(anwser(request))
