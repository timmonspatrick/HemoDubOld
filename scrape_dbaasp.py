import requests
for n in range(1, 10893):
    try:
        url = "https://dbaasp.org/api/v1?query=peptide_card&peptide_id=%s&format=json" % str(n)
        
        data = requests.get(url).text
        with open("DBAASP/DBAASP%s.html" % (str(n)), "w", encoding="utf-8") as f:
            f.write(data)
        print(str(n) + " successful")
    except:
        print(str(n) + " failed")
