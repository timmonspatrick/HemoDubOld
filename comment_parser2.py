from hemo_phrases import YES_HEMOLYTIC, LITTLE_HEMOLYTIC, NO_HEMOLYTIC, YES_HEMOLYTIK, LITTLE_HEMOLYTIK, NO_HEMOLYTIK
import re

highly_hemolytic = {300 : 50, 200 : 30, 100 : 20, 50 : 15, 20 : 10, 10 : 5}
            
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

    
def hemolytik_parser_old(comment):
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




def ugml2um(ug, mw):
    #ug in units of ug / mL
    #mw in units of g/M == ug/uM
    uM = 1000 * float(ug) / float(mw) 
    return uM

def mgml2um(mg, mw):
    #mg in units of mg / mL
    #mw in units of g/M == mg/uM
    uM = float(mg) / float(mw) 
    return uM

def avg(a, b):
    a, b = float(a), float(b)
    return a+b/2

def exponent(n):
    if n.split("x")[1] == "103":
        return float(n.split("x")[0])* 1000
    elif n.split("x")[1] == "102":
        return float(n.split("x")[0])* 100
    elif n.split("x")[1] == "10-3":
        return float(n.split("x")[0])*0.001
    elif n.split("x")[1] == "10-5":
        return float(n.split("x")[0])*0.00001
    else:
        return n

def hemolytik_parser(comment, mw):
    
    hemolytic = False
    
    if "-" == comment[0]:
        comment = comment[1:]
    
    comment = comment.replace("hemolysis upto", "hemolysis at").replace("hemolysis up to", "hemolysis at")
    comment = comment.replace("LC50", "50% hemolysis at").replace("LD50", "50% hemolysis at").replace("IC50", "50% hemolysis at").replace("HC50", "50% hemolysis at").replace("EC50", "50% hemolysis at").replace("HD50", "50% hemolysis at")
    comment = comment.replace("Lethal concentration", "50% hemolysis at")
    comment = comment.replace("No hemolysis", "0% hemolysis at 100μM").replace("no hemolysis", "0% hemolysis at 100μM").replace("negligible hemolysis", "5% hemolysis")
    comment = comment.replace("Little", "5%").replace("No or little hemolysis", "0% hemolysis").replace("Inactive", "0% hemolysis at 100μM")
    comment = comment.replace("%  hemolysis","% hemolysis").replace("hemolysis  at","hemolysis at")
    comment = comment.replace("hemolytic", "hemolysis")
    comment = comment.replace("Increased hemolysis", "20% hemolysis")
    comment = comment.replace("% at", "% hemolysis at")
    comment = comment.replace("hemolysis in", "hemolysis at")
    comment = comment.replace("fracions Hemoysis","fractions hemolysis").replace("hemolyis","hemolysis")
    comment = comment.replace("hemolysis 100µM", "hemolysis at 100µM").replace("hemolysis 100 µM","hemolysis at 100 µM")
    comment = comment.replace("hemolysis  upto", "hemolysis at").replace("Hmolysis","hemolysis")
    comment = comment.replace("hemolysis 5.0μM", "hemolysis at 5.0μM").replace("%hemolysis", "% hemolysis")
    comment = comment.replace("μ ","µ")
    comment = comment.replace("hemolytoc","hemolysis").replace("low hemolysis activity up to", "5% hemolysis at")
    comment = comment.replace("hemolysis upto", "hemolysis at").replace("hemolysis from","hemolysis at").replace("0.1 hemolysis at", "10% hemolysis at")
    
    if "fractions" in comment:
        frac, rest = comment.split("fractions")
        comment = str(float(frac) * 100) + " %" +  rest
        
    if "(" in comment and ")" in comment:
        comment = comment.split("(")[0] + comment.split(")")[1]
    
    if "% hemolysis at" in comment:    
        comment = comment.replace("µ","µ").replace("µ","µ").replace("μ","µ").replace("µ","µ").replace("mL","ml").replace("=","")

        hemopercent = comment.split("hemolysis at")[0].replace("<","").replace("%","").replace("≥","").replace("~","").replace(">","")
        conc = comment.split("hemolysis at")[1].replace("≥","").replace("above","").replace(">","").replace("~","").replace("<","")
        
        
        if conc.replace(r"µM","") != conc:  ### if "μM" in conc:
            if "±" in conc:
                conc = conc.split("±")[0]
            else:
                conc = conc.split("µM")[0]       

            unit = "µM"
            
            if "-" in conc:
                conc = avg(conc.split("-")[0], conc.split("-")[1])
        
        elif conc.replace(r"mM","") != conc:  ### if "μM" in conc:
            if "±" in conc:
                conc = conc.split("±")[0]
            unit = "mM"
            conc = conc.split("mM")[0]
            if "-" in conc:
                conc = avg(conc.split("-")[0], conc.split("-")[1])
            conc = float(conc) / 1000
            unit = "µM"
            
        elif conc.replace("µg","") != conc and "ml" in conc:
            conc = conc.split(r"µg/ml")[0]
            print(conc)
            if "-" in conc:
                conc = avg(conc.split("-")[0], conc.split("-")[1])
            elif "–" in conc:
                conc = avg(conc.split("–")[0], conc.split("–")[1])                
            conc = ugml2um(conc, mw)
            unit = "µM"

        elif conc.replace("µmol","") != conc and "L" in conc:
            if "±" in conc:
                conc = conc.split("±")[0]
            conc = conc.split(r"µmol/L")[0]
            if "-" in conc:
                conc = avg(conc.split("-")[0], conc.split("-")[1])
            
            unit = "µM"

        elif "mg/ml" in conc:
            conc = conc.split(r"mg/ml")[0]
            if "-" in conc:
                conc = avg(conc.split("-")[0], conc.split("-")[1])
            
            conc = mgml2um(conc, mw)
            unit = "µM"

        else:
            conc = conc
            unit = "µM"

        
        if conc == "NA":
            conc = 0

        if type(conc) == str:
            if "x" in conc:               
                conc = exponent(conc)

        conc = float(conc)
        
        if "-" in hemopercent:
            h1 = hemopercent.split("-")[0]
            h2 = hemopercent.split("-")[1]
            hemopercent = avg(h1, h2)
        if type(hemopercent) == str:
            if "±" in hemopercent:
                hemopercent = hemopercent.split("±")[0]
        hemopercent = float(hemopercent)

        for h_conc in highly_hemolytic:
            if conc <= h_conc:
                if hemopercent >= highly_hemolytic[h_conc]:
                    hemolytic = True

    elif "MHC" in comment:
        if "μM" in comment:
            unit = "μM"
            conc = comment.replace("MHC","").replace(">","").replace("μM","").replace("<","").replace("=","").replace("≥","").replace("~","")
        elif r'μg/ml' in comment or r'μg/mL' in comment:
            comment = comment.replace("mL","ml")
            conc = comment.replace("MHC","").replace(">","").replace("μg/ml","").replace("<","").replace("=","").replace("≥","").replace("~","")
            if "-" in conc:
                conc = avg(float(conc.split("-")[0]), float(conc.split("-")[1]))
            
            conc = ugml2um(conc, mw)
            unit = "μM"
        
        elif r'g/ml' in comment:
            unit = "g/ml"
            conc = comment.replace("MHC","").replace(r'g/ml','').replace(">","").replace(" ","").replace("=","")
            conc = exponent(conc)
            
         
        if float(conc) <= 50:
            hemolytic = True

    elif "EC25" in comment:
        hemolytic = True
        
    elif "non-hemolytic" in comment.lower() or "non-hemolysis" in comment.lower() or "Not detected hemolysis" in comment:
        pass
    
    else:
        print("ELSE " + comment)

    if hemolytic == True:
        return "YES"
    elif hemolytic == False:
        return "NO"

def dbaasp_parser(comments, mw):
    ##Consistent with the boundaries set in the HemoPi Paper
    
    hemolytic = False
    
    comments = [comment for comment in comments if "erythro" in comment["targetCell"].lower() ]
    for comment in comments:
        commentx = comment["lysis"]

        conc = comment["concentration"].replace(">=","").replace("<=","").replace(">","").replace("<","")
        conc = conc.replace("up to","")
        conc = conc.split("±")[0]
        
        if "-" in conc:
            conc = avg(conc.split("-")[0], conc.split("-")[1])
        
        if conc == "NA":
            conc = 0
            
        conc = float(conc)
        
        
        ###Unit conversion########
        unit = comment["unit"]
        
        if unit == "µg/ml":
            conc = ugml2um(conc, mw)
            unit = "µM"
        
        assert unit == "µM"
        
        if "mhc" not in commentx.lower():
            
            hemopercent = float(commentx.lower().replace("hemolysis","").replace("%","").replace("<","").replace(">","").split("±")[0].replace("(","").replace(")","").replace("cell death","").replace("ic",""))
                    
            for h_conc in highly_hemolytic:
                if conc <= h_conc:
                    if hemopercent >= highly_hemolytic[h_conc]:
                        hemolytic = True
            
        elif "mhc" in commentx.lower():
            if conc <= 50:
                hemolytic = True
            
    if hemolytic == True:
        return "YES"
    elif hemolytic == False:
        return "NO"
        
    

    