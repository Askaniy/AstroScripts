
# AstroLib, version of February 14, 2020
# Discord bot release and many improvements


import astropy.units as u
import codecs
import re
import os
import sys
import discord

token = open(".../discord_bot_token.txt").read()
path = ".../database.txt"

adjectivals = {
    "solar": "sun",
    "mercurian": "mercury",
    "venerian": "venus",
    "earthly": "earth",
    "martian": "mars",
    "cererian": "ceres",
    "jovian": "jupiter",
    "saturnian": "saturn",
    "uranian": "uranus",
    "neptunian": "neptune",
    "plutonian": "pluto"
}

abbrevs = {
    "r": "radius",
    "d": "diameter",
    "eq": "equatorial",
    "pl": "polar"
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
    "south polar radius": set(["south radius", "south polar radius"])
}

def read_request(request):
    request = [x.lower() for x in request.split(" ")]
    meaning = {}
    if "=" in request:
        request.remove("=")
    if "the" in request:
        request.remove("the")
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
        meaning.update({"body": " ".join(request[x+1:]).title()})
        parameter = request[:x]
        for i in range(len(parameter)):
            if parameter[i] in abbrevs:
                parameter[i] = abbrevs[parameter[i]]
        meaning.update({"parameter": " ".join(parameter)})
        return meaning
    elif len(request) == 2:
        if request[0] in adjectivals:
            request[0] = adjectivals[request[0]]
        meaning.update({"body": request[0].title()})
        if request[1] in abbrevs:
            request[1] = abbrevs[request[1]]
        meaning.update({"parameter": request[1]})
        return meaning
    else:
        print("Error. Try using “of” construct.")
    #print("- result of read_request: " + str(meaning))

def find_body(request, path):
    name = request
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
                    if name in names:
                        read_flag = True
                        body["name"]=names
                    else:
                        parents.insert(obj_level, names[0])
                flag = False
                if line.isspace():
                    flag = True

def find_unit(unit):
    if unit in ["km", "kilometer"]:
        return u.km
    if unit in ["m", "meter"]:
        return u.m
    if unit in ["au", "astronomical unit"]:
        return u.au
    if unit in ["cm", "centimeter"]:
        return u.cm
    if unit in ["pc", "parsec"]:
        return u.pc
    if unit in ["ly", "light years"]:
        return u.ly

def diameter_to_radius(x):
    return x[0] / 2.0
def radius_to_diameter(x):
    return x[0] * 2.0
def mean_value(x):
    return (x[0] + x[1]) / 2.0

def calculate(body, name, formuls):
    print("- request to calculate: " + name, formuls)
    crossing = list(set(body) & alt_names[name])
    if crossing != []:
        value = body[crossing[0]]["value"]
        unit = body[crossing[0]]["unit"]
        body[crossing[0]].update({"value": float(value) * find_unit(unit), "database name": crossing[0]})
        print("- result of calculate (if was found): " + str(body[crossing[0]]))
        return body[crossing[0]]
    else:
        for formula in formuls:
            need_parameter_values = []
            for need_parameter in formula[0]:
                need_parameter_values.append(find_parameter(need_parameter, body)) 
            if need_parameter_values != []:
                to_formula_value = []
                result = {"database name": name}
                for i in range(len(need_parameter_values)):
                    for j in ["source", "comment"]:
                        result.update({j: ["*calculated*"]})
                        if j in body[need_parameter_values[i]["database name"]]:
                            if body[need_parameter_values[i]["database name"]][j] not in result[j]:
                                result[j].append(body[need_parameter_values[i]["database name"]][j])
                    value = body[need_parameter_values[i]["database name"]]["value"]
                    to_formula_value.append(value)
                result.update({"value": formula[1](to_formula_value)})
                print("- result of calculate (if was calculated): " + str(result))
                return result
        print(name.capitalize() + " is unknown and cannot be calculated.")

def find_parameter(request, body):
    print("- request to find_parameter: " + request)
    if request in ["data", "parameters"]:
        return body
    elif request in ["parent", "parents"]:
        return {"value": "/".join(body["parent"])}
    elif request in alt_names["mean radius"]:
        return calculate(body, "mean radius", [(["mean diameter"], diameter_to_radius)])
    elif request in alt_names["mean diameter"]:
        return calculate(body, "mean diameter", [(["mean radius"], radius_to_diameter)])
    elif request in alt_names["mean equatorial radius"]:
        return calculate(body, "mean equatorial radius", [(["mean equatorial diameter"], diameter_to_radius)])
    elif request in alt_names["mean equatorial diameter"]:
        return calculate(body, "mean equatorial diameter", [(["mean equatorial radius"], radius_to_diameter)])
    elif request in alt_names["long equatorial radius"]:
        return calculate(body, "long equatorial radius", [(["long equatorial diameter"], diameter_to_radius)])
    elif request in alt_names["long equatorial diameter"]:
        return calculate(body, "long equatorial diameter", [(["long equatorial radius"], radius_to_diameter)])
    elif request in alt_names["short equatorial radius"]:
        return calculate(body, "short equatorial radius", [(["short equatorial diameter"], diameter_to_radius)])
    elif request in alt_names["short equatorial diameter"]:
        return calculate(body, "short equatorial diameter", [(["short equatorial radius"], radius_to_diameter)])
    elif request in alt_names["mean polar radius"]:
        return calculate(body, "mean polar radius", [(["mean polar diameter"], diameter_to_radius)])
    elif request in alt_names["mean polar diameter"]:
        return calculate(body, "mean polar diameter", [(["mean polar radius"], radius_to_diameter)])
    #elif request in alt_names["north polar radius"]:
    #    return calculate(body, "north polar radius", ["south polar radius"], radius_to_diameter)
    #elif request in alt_names["south polar radius"]:
    #    return calculate(body, "south polar radius", ["north polar radius"], radius_to_diameter)
    else:
        print("Unknown parameter.")

def find_anwser(request):
    body = find_body(request["body"], path)
    if body == None:
        return "Invalid request of celestial body."
    else:
        parameter = find_parameter(request["parameter"], body)
        if parameter == None:
            return "Invalid request of parameter."
        elif parameter == body:
            return body
        else:
            if "unit" in request:
                parameter.update({"value": parameter["value"].to(find_unit(request["unit"]))})
                if "error" in parameter:
                    parameter.pop("error")
            if request["meta"] == "info":
                anwser = "{} of {} is {}".format(parameter["database name"].capitalize(), request["body"], parameter["value"])
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

#request = read_request(input("Write your request here: "))
#print(find_anwser(request))


class BotClient(discord.Client):
    async def on_ready(self):
        print("Logged in as " + str(self.user))

    async def on_message(self, message):
        if message.content == "`BOT FINISH`":
            sys.exit()
        elif message.content == "`BOT RESTART`":
            os.execl(sys.executable, sys.executable, *sys.argv)
        elif message.author.discriminator != "1451" and "=" in message.content:
            await message.channel.send(find_anwser(read_request(message.content)))

BotClient().run(token)