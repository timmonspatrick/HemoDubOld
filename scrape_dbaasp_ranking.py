import requests
import json

activity_measures= [("373", "NotActive"), 
     ("14651", "0-10% Hemolysis"),
     ("14652", "10-20% Hemolysis"),
     ("14653", "20-30% Hemolysis"),
     ("14654", "30-40% Hemolysis"),
     ("14655", "40-50% Hemolysis"),
     ("14656", "50-60% Hemolysis"),
     ("14657", "60-70% Hemolysis"),
     ("14658", "70-80% Hemolysis"),
     ("14659", "80-90% Hemolysis"),
     ("14660", "90-100% Hemolysis"),
     ("7586", "50% Cell Death"),     
     ]

species = [("6542", "Human"),
           ("6548", "Rabbit"),
           ("6553", "Bovine"),
           ("12855", "Buffalo"),
           ("6629", "Cattle"),
           ("6639", "Chicken"),
           ("6655", "Cod"),
           ("6628", "Dog"),
           ("7206", "Erythrocytes"),
           ("13054", "Fish"),
           ("8616", "Goat",),
           ("6558", "Guinea Pig"),
           ("9102", "Hag Fish"),
           ("6551", "Horse"),
           ("6544", "Mouse"),
           ("6556", "Ovine"),
           ("6552", "Pig"),
           ("6888", "Pig2"),
           ("6631", "Rat"),
           ("10038", "Salmon"),
           ("6546", "Sheep"),
           ("6656", "Tilapia"),
           ("6601", "Trout"),
        ]

d = {n[0] : [] for n in activity_measures}
count = 0
for m in species:
    for n in activity_measures:
        url = "https://dbaasp.org/api/v1?query=ranking_search&target_cell_id=%s&activity_measure_id=%s&operation=<=&activity=1000000&format=json" % (m[0], n[0])
        
        data = requests.get(url).text
        data = json.loads(data)
        
        d[n[0]].extend(data["result"])
        count = count + len(data["result"])
        
print(d)
print(count)