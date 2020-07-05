
# AstroLib, version of July 4, 2020


from scipy import pi, sin, cos, arcsin, arccos, sqrt, log10
import astropy.constants as c
import astropy.units as u
import time
import re
import os
import sys
import discord
import json
from astropy.units import cds
#cds.enable()


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

formats = ["askaniy", "ttarrants", "ssc", "json", "decor"]

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
    "minimum mass": {
        "alt_names": set(["minimum mass", "m*sin(i)", "msini"])},
    #    "formulas": [([ , ], lambda  , )]},
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
        "alt_names": set(["rotation period", "p_rot", "prot"])},
    #    "formulas": [([ , ], lambda  , )]},
    "spin–orbit resonance": {
        "alt_names": set(["spin–orbit resonance", "tidal locking", "gravitational locking", "spin-orbit locking"])},
    #    "formulas": [([ , ], lambda  , )]},
    "orbital period": {
        "alt_names": set(["orbital period"]),
        "horizons": "P"},
    #    "formulas": [([ , ], lambda  , )]},
    "semimajor axis": {
        "alt_names": set(["semimajor axis", "semiaxis"]),
        "horizons": "a"},
    #    "formulas": [([ , ], lambda  , )]},
    "eccentricity": {
        "alt_names": set(["eccentricity"]),
        "horizons": "e"},
    #    "formulas": [([ , ], lambda  , )]},
    "inclination": {
        "alt_names": set(["inclination", "incl", "inc"]),
        "horizons": "incl"},
    #    "formulas": [([ , ], lambda  , )]},
    "longitude of ascending node": {
        "alt_names": set(["longitude of ascending node", "ascendingnode", "node"])},
    #    "formulas": [([ , ], lambda  , )]},
    "longitude of periapsis": {
        "alt_names": set(["longitude of periapsis", "longitude of pericenter", "longperi"])},
    #    "formulas": [([ , ], lambda  , )]},
    "argument of periapsis": {
        "alt_names": set(["argument of periapsis", "argument of pericenter"])},
    #    "formulas": [([ , ], lambda  , )]},
    "mean anomaly": {
        "alt_names": set(["mean anomaly", "meananomaly"])},
    #    "formulas": [([ , ], lambda  , )]},
    "temperature": {
        "alt_names": set(["temperature", "temp"])},
    #    "formulas": [([ , ], lambda  , )]},
    "luminosity": {
        "alt_names": set(["luminosity", "lum"])},
    #    "formulas": [([ , ], lambda  , )]},
    "metallicity": {
        "alt_names": set(["metallicity", "metal" "[fe/h]", "fe/h"])},
    #    "formulas": [([ , ], lambda  , )]},
    "minimum velocity": {
        "alt_names": set(["minimum velocity", "v*sin(i)", "vsini"])},
    #    "formulas": [([ , ], lambda  , )]},
    "distance": {
        "alt_names": set(["distance"])},
    #    "formulas": [([ , ], lambda  , )]},
    "parallax": {
        "alt_names": set(["parallax"])}
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
    "sun mass": u.Msun,
    "msun": u.Msun,
    "ms": u.Msun,
    "jupiter mass": u.Mjup,
    "mj": u.Mjup,
    "earth mass": u.Mearth,
    "me": u.Mearth,
    "sun radius": u.Rsun,
    "rs": u.Rsun,
    "jupiter radius": u.Rjup,
    "rj": u.Rjup,
    "earth radius": u.Rearth,
    "re": u.Rearth,
    "sec": u.s,
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
    "ppm": cds.ppm
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

def tt_class(tt_type):
    if tt_type == "system":
        return ["barycenter"]
    elif tt_type == "star":
        return ["celestial object", "stellar object", "star"]
    elif tt_type == "planet":
        return ["celestial object", "planemo", "exoplanet"]
    elif tt_type == "moon":
        return ["celestial object", "planemo", "exoplanet"]
    elif tt_type == "asteroid":
        return ["celestial object", "planemo", "exoplanet"]
    elif tt_type == "disk":
        return ["celestial object", "disk"]

def tt_unit(prmt, clss):
    if prmt in ["temp", "t_eql"]:
        return u.K
    elif prmt in ["vsini", "k_rv"]:
        return u.km / u.s
    elif prmt in ["semiaxis", "separation"]:
        return u.AU
    elif prmt in ["rotinc", "inc", "longperi", "meananomaly", "ascendingnode", "node"]:
        return u.deg
    elif prmt in ["period", "transit_dur"]:
        return u.d
    elif prmt in ["parallax", "k_astrometry"]:
        return u.mas
    elif prmt.endswith("mag"):
        return u.mag
    elif prmt == "distance":
        return u.pc
    elif prmt == "age":
        return u.yr * 10E9
    elif prmt == "lum":
        return u.Lsun
    elif prmt == "b":
        return u.Rsun
    elif prmt == "logg":
        return u.cm / u.s**2 #log10(u.cm / u.s**2)
    elif prmt == "msini":
        return u.Mjup
    if clss == "star":
        if prmt == "mass":
            return u.Msun
        elif prmt == "radius":
            return u.Rsun
        elif prmt == "p_rot":
            return u.d
    elif clss == "exoplanet":
        if prmt == "mass":
            return u.Mjup
        elif prmt == "radius":
            return u.Rjup
        elif prmt == "p_rot":
            return u.h
    return u.dimensionless_unscaled

def prmt_reader(line):
    if ":" in line:
        parameter = line.split(":")
        return {"name": parameter[0].strip(), "value": parameter[1].strip()}
    elif "=" in line:
        parameter = line.split("=")
        name = parameter[0].strip()
        value = parameter[1].strip()
        try:
            value_list = value.split(" ")
            q = unit(value_list[-1])
            value = "".join(value_list[:-1])
            if "±" in value:
                value_list = value.split("±")
                if value_list[0][-1].isalpha:
                    value_list[0] = value_list[0][:-1]
                value = float(value_list[0].strip())
                error = float(value_list[1].strip())
                return {"name": name, "value": value * q, "error": error * q}
            elif "+" in value and value[0] != "+":
                value_list = value.split("+")
                value = float(value_list[0])
                error = value_list[1].split("-")
                if len(error) == 2:
                    plus = float(error[0])
                    minus = float(error[1])
                elif len(error) == 4:
                    plus = float(error[0] + "-"+error[1])
                    minus = float(error[2] + "-"+error[3])
                return {"name": name, "value": value * q, "error": (plus * q, minus * q)}
            elif "-" in value and value[0] != "-":
                value_list = value.split("-")
                min_v = float(value_list[0])
                max_v = float(value_list[1])
                value = (min_v + max_v) / 2
                error = (max_v - min_v) / 2
                return {"name": name, "value": value * q, "error": error * q}
            elif value[-1] == "%":
                value = float(value.replace("%", "").strip())
                return {"name": name, "value": value / 100}
            else:
                return {"name": name, "value": float(value) * q}
        except IndexError:
            return {"name": parameter[0].strip(), "value": float(parameter[1].strip())}
    else:
        name = line.strip()
        value = "use `:` for text and `=` for quantity"
    return {"name": name, "value": value}

def find_body(request, database_path):
    global cannot_calc
    cannot_calc = []
    body = {}
    files = [f.path for f in os.scandir(database_path) if f.is_file()]
    #files.extend([f.path for f in os.scandir(discord_path) if f.is_file()])
    files.reverse()
    for fl in files:
        if os.path.splitext(fl)[1] == ".askaniy":
            with open(fl, errors="ignore") as f:
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
                head = []
                for line in f_list:
                    if not system:
                        if not read_obj_prmt:
                            if read_obj_name and not line.isspace():
                                obj_level = line.count("\t")
                                names = line.strip().split(", ")
                                noun_version = list(map(lambda i: " ".join(list(map(lambda j: noun(j.lower()), i.split(" ")))), names))
                                if request in noun_version:
                                    read_obj_prmt = True
                                    body["names"] = names
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
                                        head = []
                                        if "<" in line:
                                            read_group = True
                                            read_head = True
                                            head.append(prmt_reader(line.replace("<", "")))
                                        else:
                                            parameter = prmt_reader(line) # < parameter reading block
                                            if line.count("\t") == obj_level:
                                                body[parameter["name"]] = {"value": parameter["value"]}
                                                if "error" in parameter:
                                                    body[parameter["name"]].update({"error": parameter["error"]})
                                                last_parameter = parameter["name"]
                                            else:
                                                body[last_parameter].update({parameter["name"]: parameter["value"]}) # >
                                    else:
                                        if read_head and line.count("\t") == obj_level + 1:
                                            head.append(prmt_reader(line))
                                        else:
                                            read_head = False
                                            if ">" in line:
                                                read_group = False
                                                line = line.replace(">", "")
                                            parameter = prmt_reader(line) # < parameter reading block
                                            if line.count("\t") == obj_level:
                                                body[parameter["name"]] = {"value": parameter["value"]}
                                                if "error" in parameter:
                                                    body[parameter["name"]].update({"error": parameter["error"]})
                                                for h in range(len(head)):
                                                    body[parameter["name"]].update({head[h]["name"]: head[h]["value"]})
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
        
        elif os.path.splitext(fl)[1] == ".ttarrants":
            with open(fl, errors="ignore") as f:
                f_list = list(f)
                f_list.extend(["\nSystem:Endgame"])
                part_of = []
                includes = []
                aliases = [] # system names
                names = [] # object names
                collecting_subsystem = False
                lastnotspace = False
                system = False
                for line in f_list:
                    if line.isspace():
                        if lastnotspace:
                            if system:
                                if request in [aliase.lower() for aliase in aliases]:
                                    body.update({"names": aliases})
                                    collecting_subsystem = True
                                system = False
                            else:
                                includes.append(names[0])
                                if request in [name.lower() for name in names]:
                                    return body
                            lastnotspace = False
                    else:
                        dedoted = line.strip().split(":")
                        if dedoted[0].lower() == "system":
                            if not collecting_subsystem:
                                body = {"class": tt_class("system")}
                                system = True
                                aliases = [dedoted[1]]
                                part_of = []
                                includes = []
                            else:
                                collecting_subsystem = False
                                body.update({"part of": [], "includes": includes})
                                return body
                        elif dedoted[0].lower() == "aliases":
                            aliases.extend(dedoted[1].split(","))
                        elif dedoted[0].lower() in ["star", "planet", "moon", "asteroid", "disk"]:
                            confirmed = True
                            if dedoted[-1][-1] == "?":
                                confirmed = False
                                dedoted[-1] = dedoted[-1].replace("?", "")
                            i = " " + " ".join(dedoted[2:]) # with star name variant
                            names = [aliase + i for aliase in aliases]
                            if len(dedoted) > 3:
                                j = " " + " ".join(dedoted[3:]) # w/o star name variant
                                names.extend([aliase + j for aliase in aliases])
                            if len(dedoted) > 4:
                                part_of = [dedoted[1] + " system", dedoted[1] + " " + dedoted[-2] + " system"]
                            else:
                                part_of = [dedoted[1] + " system"]
                            if not collecting_subsystem:
                                body = {"names": names, "part of": part_of, "class": tt_class(dedoted[0].lower())}
                                if not confirmed:
                                    body["class"].append("not confirmed")
                        elif dedoted[0].lower() == "altname":
                            i = " " + " ".join(dedoted[2:])
                            names = [aliase + i for aliase in aliases]
                            if not collecting_subsystem:
                                body["names"].extend(names)
                        elif dedoted[0].lower() == "type":
                            if not collecting_subsystem:
                                body["class"].extend(dedoted[1].lower())
                        elif dedoted[0].lower() in ["reference", "published"]:
                            if not collecting_subsystem:
                                body.update({dedoted[0].lower(): {"value": line.strip().split(":", 1)[1]}})
                        elif dedoted[0].lower() in [
                            "note", "association", "field", "paper", "announce", "announced", "announcement", 
                            "update", "discovery", "discoverer", "discoverers", "dmethod", "detectionmethods", 
                            "confirm", "resonance", "circumbinary", "epoch", "t_peri", "ra", "dec", "spectral", "atmosphere"]:
                            if not collecting_subsystem:
                                if len(dedoted) == 2:
                                    body.update({dedoted[0].lower(): {"value": dedoted[1]}})
                                else:
                                    body.update({dedoted[0].lower(): {"value": dedoted[1], "comment": dedoted[2:]}})
                        else:
                            if not collecting_subsystem:
                                prmt_name = dedoted[0].lower()
                                if len(dedoted) == 2:
                                    body.update({prmt_name: {}})
                                else:
                                    body.update({prmt_name: {"comment": dedoted[2:]}})
                                check_text = dedoted[1].replace(" ", "").replace("+", "").replace("-", "").replace("/", "")
                                if check_text.isalpha():
                                    body[prmt_name].update({"value": dedoted[1]})
                                else:
                                    if "_" in dedoted[1]:
                                        q = unit(dedoted[1].lower().split("_")[0])
                                        prmt_value = dedoted[1].lower().replace("=", "").split("_")[1]
                                    else:
                                        q = tt_unit(dedoted[0].lower(), body["class"][-1])
                                        prmt_value = dedoted[1].lower().replace("=", "")
                                    if "?" in prmt_value:
                                        body[prmt_name].update({"comment": "?"})
                                        prmt_value = prmt_value.replace("?", "")
                                    if prmt_value[0] == ">":
                                        body[prmt_name].update({"comment": "min value"})
                                        prmt_value = prmt_value.replace(">", "")
                                    elif prmt_value[0] == "<":
                                        body[prmt_name].update({"comment": "max value"})
                                        prmt_value = prmt_value.replace("<", "")
                                    if "±" in prmt_value:
                                        text = prmt_value.split("±")
                                        if text[0][-1].isalpha:
                                            text[0] = text[0][:-1]
                                        if text[1] != "":
                                            body[prmt_name].update({"value": float(text[0]) * q, "error": float(text[1]) * q})
                                        else:
                                            body[prmt_name].update({"value": float(text[0]) * q})
                                    elif "+" in prmt_value and prmt_value[0] != "+":
                                        text = prmt_value.split("+")
                                        error = text[1].split("-")
                                        if len(error) == 2:
                                            try:
                                                plus = float(error[0]) * q
                                            except ValueError:
                                                plus = 0 * q
                                            try:
                                                minus = float(error[1]) * q
                                            except ValueError:
                                                minus = 0 * q
                                        elif len(error) == 4:
                                            plus = float(error[0]+"-"+error[1]) * q
                                            minus = float(error[2]+"-"+error[3]) * q
                                        body[prmt_name].update({"value": float(text[0]) * q, "error": (plus, minus)})
                                    elif "-" in prmt_value and prmt_value[0] != "-":
                                        text = prmt_value.split("-")
                                        min_v = float(text[0])
                                        max_v = float(text[1])
                                        value = (min_v + max_v) / 2
                                        error = (max_v - min_v) / 2
                                        body[prmt_name].update({"value": value * q, "error": error * q})
                                    elif prmt_value[-1] == "%":
                                        prmt_value = float(prmt_value.replace("%", ""))
                                        body[prmt_name].update({"value": prmt_value / 100})
                                    else:
                                        body[prmt_name].update({"value": float(prmt_value) * q})
                        lastnotspace = True

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
            #result.update({"comment": ["*theoretical*"]})
            #body.update({name: result})
            return process(prmtr[name]["included"], body)
        #elif "standart" in prmtr[name]:
    else:
        body[cross].update({"database name": cross})
        print("- result of process (if was found): " + str(body[cross]))
        return body[cross]

def find_parameter(request, body):
    if request in ["part of", "parent", "parents"]:
        return {"database name": "parents", "value": ", ".join(body["part of"])}
    elif request in ["system", "include", "includes"]:
        return {"database name": "system", "value": ", ".join(body["includes"])}
    elif request == "class":
        return find_class("askaniy", body)
    else:
        for name in prmtr:
            if request in prmtr[name]["alt_names"]:
                return process(name, body)

def find_class(standard, body):
    if "class" in body:
        if standard == "askaniy":
            return ", ".join(body["class"])
        elif standard == "ttarrants":
            if "barycenter" in body["class"]:
                return "System"
            elif "star" in body["class"]:
                return "Star"
            elif "planemo" in body["class"]:
                return "Planet"
            elif "disk" in body["class"]:
                return "Disk"
            else:
                return "UFO"
        elif standard == "ssc":
            radius = find_parameter("mean radius", body)["value"]
            if radius < u.Quantity("450 km"):
                return "asteroid"
            elif radius < u.Quantity("2400 km"):
                return "dwarfplanet"
            else:
                return "planet"
    else:
        return "Work in progress"


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
                r += f'\n{value} - {key}'
            return f'Supported adjectives:{r}'
    elif request[0] == "abbreviations":
        if bot:
            em = discord.Embed(title = "Supported abbreviations and replacements:", color = discord.Colour.from_rgb(127, 127, 255))
            for key, value in abbrevs.items():
                em.add_field(name = value, value = key, inline = True)
        else:
            r = ""
            for key, value in abbrevs.items():
                r += f'\n{key} - {value}'
            return f'Supported abbreviations and replacements:{r}'
    elif request[0] == "formats":
        if bot:
            em = discord.Embed(title = "Supported formats:", description = "\n".join(formats), color = discord.Colour.from_rgb(127, 127, 255))
        else:
            return "Supported formats:\n* " + "\n* ".join(formats)
    else:
        return f'I can’t help with {" ".join(request)}'

def create(body):
    file_name = time.strftime("%Y-%m-%d_%H-%M-%S")
    with open("{}/{}.askaniy".format(discord_path, file_name), "w") as f:
        f.write(body)
    body_name = body.split("\n")[0]
    print(body_name.split(", ")[0])
    try:
        print(find_body(body_name.split(", ")[0], database_path))
        return "I saved your object successfully, you can work with it"
    except Exception:
        return "Some error in your syntax, I can’t read it"

def output(parameter):
    try:
        anwser = str(parameter["value"].value)
        if "error" in parameter:
            try:
                if type(parameter["error"]) == tuple:
                    anwser += f' +{parameter["error"][0].value} -{parameter["error"][1].value}'
                else:
                    anwser += f' ± {parameter["error"].value}'
            except AttributeError:
                if type(parameter["error"]) == tuple:
                    anwser += f' +{parameter["error"][0]} -{parameter["error"][1]}'
                else:
                    anwser += f' ± {parameter["error"]}'
        anwser += " " + str(parameter["value"].unit)
    except AttributeError:
        anwser = str(parameter["value"])
    if "source" in parameter:
        if type(parameter["source"]) == list:
            if len(parameter["source"]) == 1:
                anwser += f'; source: {parameter["source"][0]}'
            else:
                anwser += f'; sources: {", ".join(parameter["source"])}'
        else:
            anwser += f'; source: {parameter["source"]}'
    if "comment" in parameter:
        if type(parameter["comment"]) == list:
            anwser += f' ({", ".join(parameter["comment"])})'
        else:
            anwser += f' ({parameter["comment"]})'
    return anwser

def short_list(system):
    res = []
    for component in system:
        res.append(component.split(", ")[0])
    if len(res) >= 20:
        res = res[:19]
        res.append("…")
    return res

def embed(body):
    global em
    em = discord.Embed(
        title = ", ".join(body["names"]),
        description = "class: " + find_class("askaniy", body),
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
        return "I can’t understand this parameter"
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
            return f'{parameter["database name"].capitalize()} of {body["names"][0].capitalize()} is {output(parameter)}'
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
    if request == "askaniy":
        level = ""
        systems = ""
        for system in body["part of"]:
            systems += f'{level}{system}\n\n'
            level += "\t"
        code = f'{systems}{level}{", ".join(body["names"])}\n'
        for parameter in body:
            if parameter not in ["names", "part of", "includes"]:
                code += f'{level}{parameter} = {body[parameter]["value"]}\n'
        return code
    elif request == "ttarrants":
        code = f'{find_class("ttarrants", body)}:{",".join(body["names"])}\n'
        for parameter in body:
            if parameter not in ["names", "part of", "includes"]:
                code += f'{parameter}:{body[parameter]["value"]}\n'
        return code
    elif request == "ssc":
        try:
            code = '"{}" "{}"\n'.format(":".join(body["names"]), "/".join(body["part of"]))
            code += "{\n"
            code += f'\tClass\t"{find_class("ssc", body)}"\n'
            code += f'\tTexture\t"{body["names"][0].lower()}.*"\n'
            code += f'\tRadius\t{find_parameter("mean radius", body)["value"].to(u.km).value}\n'
            code += "}\n"
            return code
        except Exception:
            return f'{body["names"][0]} doesn’t have all the basic parameters for creating a SSC file.'
    elif request == "json":
        for i in body:
            for j in body[i]:
                try: 
                    body[i].update({j: str(body[i][j])})
                except Exception:
                    pass
        return json.dumps(body)
    else:
        decor = f'\n{", ".join(body["names"])}\nclass: {find_class("askaniy", body)}\n'
        if body["part of"] != []:
            decor += f'\npart of:\n{", ".join(body["part of"])}\n'
        if "includes" in body:
            decor += f'\nincludes:\n{", ".join(short_list(body["includes"]))}\n'
        for key, value in body.items():
            if type(value) != list:
                decor += (f"\n{key}\n{output(value)}\n")
        return decor


# First layer of the script:

def anwser(request):
    print(f'- result of read_request: {request}')
    global em
    em = None
    if request:
        if request["to_do"] == "help":
            return hlp(request["meta"])
        elif request["to_do"] == "create":
            return create(request["body"])
        body = find_body(request["body"], database_path)
        if body == None:
            return f'I can’t find {request["body"]} in the available databases'
        elif request["to_do"] == "return":
            if bot:
                return f'`{body}`'
            else:
                return str(body)
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
        return "I can’t do it for now"


# Surface of the script:

if bot:
    class BotClient(discord.Client):
        async def on_ready(self):
            print(f'Logged in as {self.user}')
        async def on_message(self, message):
            request = [x.lower() for x in message.content.split(" ")]
            if request[0] in supported_requests and message.author.discriminator != "1451":
                await message.channel.send(anwser(read_request(request)), embed = em)
    BotClient().run(open(token_path).read())
else:
    request = check_request(input("Write your request here: "))
    print(anwser(request))
