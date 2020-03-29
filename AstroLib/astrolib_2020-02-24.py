
# AstroLib, version of February 24, 2020


import astropy.units as u
import codecs
import re
import os
import sys
import discord

bot = False

path = "X:/Documents/GitHub/AstroScripts/AstroLib/Database/universe.txt"
token_path = "X:/Documents/GitHub/ExternalData/discord_bot_token.txt"

links = {
    "WGCCRE2015": "https://astropedia.astrogeology.usgs.gov/download/Docs/WGCCRE/WGCCRE2015reprint.pdf",
    "IAU2015 Resolution B3": "http://www.astro.osu.edu/~pogge/Ast2292/Docs/NominalConstants.pdf",
    "JPL Planetary Satellite Physical Parameters": "https://ssd.jpl.nasa.gov/?sat_phys_par"
}

alt_names = {
    "mean radius": set(["radius", "mean radius"]),
    "mean diameter": set(["diameter", "mean diameter"]),
    "mean equatorial radius": set(["equatorial radius", "mean equatorial radius"]),
    "mean equatorial diameter": set(["equatorial diameter", "mean equatorial diameter"]),
    "long equatorial radius": set([
        "x radius", "x axis", "long radius", "long equatorial radius", 
        "subplanetary radius", "subplanetary equatorial radius"]),
    "long equatorial diameter": set([
        "x diameter", "long diameter", "long equatorial diameter", 
        "subplanetary diameter", "subplanetary equatorial diameter"]),
    "short equatorial radius": set([
        "y radius", "y axis", "short radius", "short equatorial radius", 
        "along orbit equatorial radius", "along orbit radius"]),
    "short equatorial diameter": set([
        "y diameter", "short diameter", "short equatorial diameter", 
        "along orbit equatorial diameter", "along orbit diameter"]),
    "mean polar radius": set(["z radius", "z axis", "polar radius", "mean polar radius"]),
    "mean polar diameter": set(["z diameter", "polar diameter", "mean polar diameter"]),
    "north polar radius": set(["north radius", "north polar radius"]),
    "south polar radius": set(["south radius", "south polar radius"]),
    "mean oblateness": set(["oblateness", "mean oblateness"]),
    "equatorial oblateness": set(["equatorial oblateness"]),
    #"surface area": set(["surface area"])
    #"volume": set(["volume"])
    #"standard gravitational parameter": set(["mass parameter", "gravitational parameter", "standard gravitational parameter"])
}

abbrevs = {
    "r": "radius",
    "d": "diameter",
    "eq": "equatorial",
    "pl": "polar",
    "obl": "oblateness",
    "max": "maximum",
    "min": "minimum"
}

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

def read_request(request):
    request = [x.lower() for x in request.split(" ")]
    meaning = {}
    if "=" in request:
        request.remove("=")
    if "→" in request:
        request.remove("→")
    if "the" in request:
        request.remove("the")
    if request[0] == "return":
        meaning.update({"to_do": "return", "body": " ".join(request[1:])})
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
            meaning.update({"unit": request[x+1]})
            request = request[:x]
        if "of" in request:
            x = request.index("of")
            meaning.update({"body": " ".join(request[x+1:])})
            parameter = request[:x]
            for i in range(len(parameter)):
                if parameter[i] in abbrevs:
                    parameter[i] = abbrevs[parameter[i]]
            meaning.update({"parameter": " ".join(parameter)})
            return meaning
        elif len(request) == 2:
            meaning.update({"body": noun(request[0])})
            if request[1] in abbrevs:
                request[1] = abbrevs[request[1]]
            meaning.update({"parameter": request[1]})
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
        if meaning["format"] in ["ssc"]:
            return meaning
    #print("- result of read_request: " + str(meaning))

def find_body(request, path):
    body = {}
    with codecs.open(path) as f:
        obj_level = 0
        parents = []
        last_parameter = ""
        flag = False
        read_flag = False
        for line in f:
            if read_flag:
                if line.isspace():
                    body.update({"parent": parents[:obj_level]})
                    #print("- result of find_body: " + str(body))
                    return body
                else:
                    parameter = line.strip().split(" = ")
                    if line.count("\t") == obj_level:
                        body[parameter[0]] = {"value": parameter[1]}
                        last_parameter = parameter[0]
                    else:
                        body[last_parameter].update({parameter[0]: parameter[1]})
            else:
                if flag:
                    obj_level = line.count("\t")
                    names = line.strip().split(", ")
                    noun_version = list(map(lambda i: " ".join(list(map(lambda j: noun(j.lower()), i.split(" ")))), names))
                    if request in noun_version:
                        read_flag = True
                        body["name"] = names
                    else:
                        parents.insert(obj_level, names[0])
                flag = False
                if line.isspace():
                    flag = True

def find_unit(unit):
    if unit in ["cm", "centimeter", "centimeters"]:
        return u.cm
    elif unit in ["m", "meter", "meters"]:
        return u.m
    elif unit in ["km", "kilometer", "kilometers"]:
        return u.km
    elif unit in ["au", "astronomical unit", "astronomical units"]:
        return u.au
    elif unit in ["pc", "parsec", "parsecs"]:
        return u.pc

def diameter_to_radius(x):
    return x[0] / 2.0
def radius_to_diameter(x):
    return x[0] * 2.0
def mean_value_2d(x):
    return (x[0] + x[1]) / 2.0
def mean_value_3d(x):
    return (2.0 * x[0] + x[1]) / 3.0
def all_minus_part(x):
    return x[0] - x[1]
def oblateness(x):
    return 1.0 - x[1]/x[0]
def long_axis(x):
    return x[0] / (1.0 - x[1])
def short_axis(x):
    return x[0] * (1.0 - x[1])

def teor_model(parameter, body):
    print("- request to teor_model: " + parameter)
    if parameter == "oblateness":
        return {"value": 0.}
    elif parameter == "equatorial oblateness":
        return {"value": 0.}
    else:
        return {}

searched = []

def calculate(body, name, formulas):
    print("- request to calculate: " + name, formulas)
    global searched
    crossing = list(set(body) & alt_names[name])
    if crossing != []:
        value = body[crossing[0]]["value"]
        unit = body[crossing[0]]["unit"]
        body[crossing[0]].update({"value": float(value) * find_unit(unit), "database name": crossing[0]})
        print("- result of calculate (if was found): " + str(body[crossing[0]]))
        return body[crossing[0]]
    else:
        for formula in formulas:
            need_parameter_values = []
            for need_parameter in formula[0]:
                if need_parameter not in searched:
                    searched.append(need_parameter)
                    need_parameter_values.append(find_parameter(need_parameter, body))
            if need_parameter_values != []:
                to_formula_value = []
                result = {"database name": name}
                for i in range(len(need_parameter_values)):
                    for j in ["source", "comment"]:
                        result.update({j: ["*calculated*"]})
                        if j in body[need_parameter_values[i]["database name"]]:
                            if body[need_parameter_values[i]["database name"]][j] not in result[j]:
                                if type(body[need_parameter_values[i]["database name"]][j]) == list:
                                    result[j].extend(body[need_parameter_values[i]["database name"]][j][1:])
                                else:
                                    result[j].append(body[need_parameter_values[i]["database name"]][j])
                    value = body[need_parameter_values[i]["database name"]]["value"]
                    to_formula_value.append(value)
                result.update({"value": formula[1](to_formula_value)})
                body.update({name: result})
                print("- result of calculate (if was calculated): " + str(result))
                return result
        result = teor_model(name, body)
        return result.update({"comment": "*teoretical*"})
        #print(name.capitalize() + " is unknown and cannot be calculated.")

def find_parameter(request, body):
    print("- request to find_parameter: " + request)
    if request in ["parent", "parents"]:
        return {"database name": "parent", "value": "/".join(body["parent"])}
    
    elif request in alt_names["mean radius"]:
        return calculate(body, "mean radius", [
            (["mean equatorial radius", "mean polar radius"], mean_value_3d), 
            (["mean diameter"], diameter_to_radius)])
    elif request in alt_names["mean diameter"]:
        return calculate(body, "mean diameter", [(["mean radius"], radius_to_diameter)])
    
    elif request in alt_names["mean equatorial radius"]:
        return calculate(body, "mean equatorial radius", [
            (["mean equatorial diameter"], diameter_to_radius), 
            (["long equatorial radius", "short equatorial radius"], mean_value_2d), 
            (["mean polar radius", "oblateness"], long_axis)])
    elif request in alt_names["mean equatorial diameter"]:
        return calculate(body, "mean equatorial diameter", [(["mean equatorial radius"], radius_to_diameter)])
    
    elif request in alt_names["long equatorial radius"]:
        return calculate(body, "long equatorial radius", [
            (["long equatorial diameter"], diameter_to_radius),
            (["short equatorial radius", "equatorial oblateness"], long_axis)])
    elif request in alt_names["long equatorial diameter"]:
        return calculate(body, "long equatorial diameter", [(["long equatorial radius"], radius_to_diameter)])
    
    elif request in alt_names["short equatorial radius"]:
        return calculate(body, "short equatorial radius", [
            (["short equatorial diameter"], diameter_to_radius),
            (["long equatorial radius", "equatorial oblateness"], short_axis)])
    elif request in alt_names["short equatorial diameter"]:
        return calculate(body, "short equatorial diameter", [(["short equatorial radius"], radius_to_diameter)])
    
    elif request in alt_names["mean polar radius"]:
        return calculate(body, "mean polar radius", [
            (["mean polar diameter"], diameter_to_radius), 
            (["long polar radius", "short polar radius"], mean_value_2d), 
            (["mean equatorial radius", "oblateness"], short_axis)])
    elif request in alt_names["mean polar diameter"]:
        return calculate(body, "mean polar diameter", [(["mean polar radius"], radius_to_diameter)])
    
    elif request in alt_names["north polar radius"]:
        return calculate(body, "north polar radius", [(["mean polar diameter", "south polar radius"], all_minus_part)])
    elif request in alt_names["south polar radius"]:
        return calculate(body, "south polar radius", [(["mean polar diameter", "north polar radius"], all_minus_part)])

    elif request in alt_names["mean oblateness"]:
        return calculate(body, "mean oblateness", [(["mean equatorial radius", "mean polar radius"], oblateness)])
    elif request in alt_names["equatorial oblateness"]:
        return calculate(body, "equatorial oblateness", [(["long equatorial radius", "short equatorial radius"], oblateness)])

    #elif request in alt_names["surface area"]:
    #    return calculate(body, "surface area", [(["mean polar diameter", "north polar radius"], all_minus_part)])
    #elif request in alt_names["volume"]:
    #    return calculate(body, "volume", [(["mean polar diameter", "north polar radius"], all_minus_part)])
    
    else:
        print("Unknown parameter.")

def find_anwser(request):
    body = find_body(request["body"], path)
    if body == None:
        print("- was requested " + str(request["body"]))
        return "Invalid request of celestial body."
    elif request["to_do"] == "return":
        return "`" + str(body) + "`"
    elif request["to_do"] == "find":
        parameter = find_parameter(request["parameter"], body)
        if parameter == None:
            return "Invalid request of parameter."
        else:
            if "unit" in request:
                parameter.update({"value": parameter["value"].to(find_unit(request["unit"]))})
                if "error" in parameter:
                    parameter.pop("error")
            if request["meta"] == "info":
                anwser = "{} of {} is {}".format(parameter["database name"].capitalize(), body["name"][0], parameter["value"])
                if "error" in parameter:
                    anwser += " ± " + str(parameter["error"])
                if "source" in parameter:
                    if type(parameter["source"]) == list:
                        anwser += " ({})".format(", ".join(parameter["source"]))
                    else:
                        anwser += " ({})".format(parameter["source"])
                if "comment" in parameter:
                    if type(parameter["comment"]) == list:
                        anwser += "  // " + ", ".join(parameter["comment"])
                    else:
                        anwser += "  // " + parameter["comment"]
                return anwser
            elif request["meta"] == "value":
                return parameter["value"]
            elif request["meta"] == "error" and "error" in parameter:
                return parameter["error"]
            elif request["meta"] == "unit" and "unit" in parameter:
                return parameter["unit"]
            elif request["meta"] == "source" and "source" in parameter:
                return parameter["source"]
            elif request["meta"] == "comment" and "comment" in parameter:
                return parameter["comment"]
            else:
                return request["meta"].title() + " is not in the requested parameter."
    elif request["to_do"] == "generate":
        if request["format"] == "ssc":
            code = '"{}" "{}"\n'.format(":".join(body["name"]), "/".join(body["parent"])) + "{"
            code += "\n\t"
            return code

if bot:
    class BotClient(discord.Client):
        async def on_ready(self):
            print("Logged in as " + str(self.user))
        async def on_message(self, message):
            if message.content == "`BOT FINISH`":
                sys.exit()
            elif message.content == "`BOT RESTART`":
                os.execl(sys.executable, sys.executable, *sys.argv)
            elif message.author.discriminator != "1451" and "→" in message.content:
                await message.channel.send(find_anwser(read_request(message.content)))
    BotClient().run(open(token_path).read())
else:
    request = read_request(input("Write your request here: "))
    print(find_anwser(request))