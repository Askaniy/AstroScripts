
# AstroLib, version of February 11, 2020
# The first completed version of the script, without support for the bot mode


import codecs
import re

path = ".../Database/database.txt"

def read_request(request):
    request = [x.lower() for x in request]
    if "in" in request:
        request.remove("in")
    if "the" in request:
        request.remove("the")
    if "of" in request:
        request.remove("of")
        request[0], request[1] = request[1], request[0]
    if len(request) == 2:
        request.insert(2, "info")
    request[0] = request[0].title()
    return request

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

def find_parameter(request, body):
    if request == "data" or request == "parameters":
        return body
    elif request == "parent" or request == "parents":
        return {"value": "/".join(body["parent"])}
    elif request == "radius":
        if "radius" in body:
            return body["radius"]
        elif "radius" not in body and "diameter" in body:
            radius = {"value": float(body["diameter"]["value"])/2}
            if "error" in body["radius"]:
                radius.update({"error": float(body["diameter"]["error"])/2})
            if "unit" in body["radius"]:
                radius.update({"unit": body["diameter"]["unit"]})
            if "source" in body["radius"]:
                radius.update({"source": body["diameter"]["source"]})
            return radius
        else:
            print("Radius cannot be calculated.")
    elif request == "diameter":
        if "diameter" in body:
            return body["diameter"]
        elif "diameter" not in body and "radius" in body:
            diameter = {"value": float(body["radius"]["value"])*2}
            if "error" in body["radius"]:
                diameter.update({"error": float(body["radius"]["error"])*2})
            if "unit" in body["radius"]:
                diameter.update({"unit": body["radius"]["unit"]})
            if "source" in body["radius"]:
                diameter.update({"source": body["radius"]["source"]})
            return diameter
        else:
            print("Diameter cannot be calculated.")
    else:
        print("Unknown parameter.")

def find_anwser(request):
    body = find_body(request[0], path)
    if body == None:
        return "Invalid request of celestial body."
    parameter = find_parameter(request[1], body)
    if parameter == None:
        return "Invalid request of parameter."
    elif parameter == body:
        return body
    else:
        if request[2] == "info":
            anwser = "{} of {} is {}".format(request[1].title(), request[0], parameter["value"])
            if "error" in parameter:
                anwser += " ± " + str(parameter["error"])
            if "unit" in parameter:
                anwser += " " + parameter["unit"]
            if "source" in parameter:
                anwser += " ({})".format(parameter["source"])
            return anwser
        elif request[2] == "value":
            return parameter["value"]
        elif request[2] == "error" and "error" in parameter:
            return parameter["error"]
        elif request[2] == "unit" and "unit" in parameter:
            return parameter["unit"]
        elif request[2] == "source" and "source" in parameter:
            return parameter["source"]
        elif request[2] in ["m", "km", "au", "pc"]:
            return "Can't work with units now."
        else:
            return request[2].title() + " is not in the requested parameter."

request = read_request(input("Write request here: ").split(" "))
print(find_anwser(request))