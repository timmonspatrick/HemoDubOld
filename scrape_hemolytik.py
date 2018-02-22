import requests
for n in range(1001, 4202):
    try:
        url = "http://crdd.osdd.net/raghava/hemolytik/display.php?details=%s" %(str(n))
        data = requests.get(url)
        data.encoding = "utf-8"
        data = data.text
        with open("HEMOLYTIK%s.html" % (str(n)), "w", encoding="utf-8") as f:
            f.write(data)
        print(str(n) + " successful")
    except:
        assert 1==2
        print(str(n) + " failed")
