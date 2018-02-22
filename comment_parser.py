from hemo_phrases import YES_HEMOLYTIC, LITTLE_HEMOLYTIC, NO_HEMOLYTIC, YES_HEMOLYTIK, LITTLE_HEMOLYTIK, NO_HEMOLYTIK
import re

def parser(comment):
    for phrase in YES_HEMOLYTIC:
        if phrase in comment:
            return "YES"
    for phrase in LITTLE_HEMOLYTIC:
        if phrase in comment:
            return "LITTLE"
    for phrase in NO_HEMOLYTIC:
        if phrase in comment:
            return "NO"

    return hemolytik_parser(comment)

    
def hemolytik_parser(comment):
    comment = comment.replace("  "," ")
        
    if "% hemolysis" in comment and not re.match("^0% hemolysis$", comment):
        return "YES"
    if "% hemolytic" in comment and not re.match("^0% hemolytic$", comment):
        return "YES"
    if "% of hemolysis" in comment and not re.match("^0% of hemolysis$", comment):
        return "YES"

    if "fractions hemolysis" in comment and not re.match("^0.0 fractions hemolysis$", comment):
        return "YES"
    
    for phrase in YES_HEMOLYTIK:
        if phrase in comment:
            return "YES"
    
    for phrase in LITTLE_HEMOLYTIK:
        if phrase in comment:
            return "LITTLE"

    for phrase in NO_HEMOLYTIK:
        if phrase in comment:
            return "NO"
    return None


def dbaasp_parser_old(comments):
    hemolytic = False
    for comment in comments:
        print(comment)
        commentx = comment["lysis"]
        if parser(commentx) == "YES":
            hemolytic = True
        elif parser(commentx) == "NO":
            pass
        elif parser(commentx) == None:
            print('"'+commentx+'",')
    
    if hemolytic == True:
        return "YES"
    elif hemolytic == False:
        return "NO"

def ugml2um(ug, mw):
    #ug in units of ug / mL
    #mw in units of g/M == ug/uM
    uM = 1000 * float(ug) / float(mw) 
    return uM

def avg(a, b):
    a, b = float(a), float(b)
    return a+b/2

def dbaasp_parser(comments, mw):
    ##Consistent with the boundaries set in the HemoPi Paper
    
    hemolytic = False
    for comment in comments:
        commentx = comment["lysis"]
        #print(comment["concentration"] + "   " + comment["unit"])
        #print(commentx)
        #print("____________________")
        
        conc = comment["concentration"].replace(">=","").replace("<=","").replace(">","").replace("<","")
        conc = conc.replace("up to","")
        conc = conc.split("±")[0]
        
        if "-" in conc:
            conc = avg(conc.split("-")[0], conc.split("-")[1])
        
        if conc == "NA":
            conc = 0
        conc = float(conc)
        
        unit = comment["unit"]
        
        if unit == "µg/ml":
            conc = ugml2um(conc, mw)
            unit = "µM"
        
        assert unit == "µM"
        print(conc)
        if conc <= 100:
            if "Hemolysis" in commentx:
                hemopercent = float(commentx.lower().replace("hemolysis","").replace("%","").replace("<","").replace(">","").split("±")[0].replace("(","").replace(")",""))
                if hemopercent <= 10:
                    pass
                else:
                    hemolytic = True
        elif conc > 100:
            pass
    
    print(hemolytic)
    if hemolytic == True:
        return "YES"
    elif hemolytic == False:
        return "NO"
        
    

    