import requests
def scrape(start, stop):
    for n in range(start, stop):
    
        url = "http://dramp.cpu-bioinfor.org/browse/All_Information.php?id=DRAMP%s&dataset=General" % ("00000" + str(n))[-5:]
        data = requests.get(url)
        data.encoding = "utf-8"
        data = data.text
        with open("DRAMP%s.html" % ("00000" + str(n))[-5:], "w", encoding="utf-8") as f:
            f.write(data)
        print(str(n) + " successful")


scrape(1, 4719)
