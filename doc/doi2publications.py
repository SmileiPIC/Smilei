from requests import get
from json import loads

DOIs = [
"10.1017/S0022377821001070",
"10.1063/5.0062503",
"10.1088/1361-6587/ac318d",
"10.3847/2041-8213/ac30e0",
"10.1016/j.physrep.2021.09.001",
"10.1088/1361-6587/ac2996",
"10.3847/2041-8213/ac1795",
"10.1088/1367-2630/ac1a97",
"10.1063/5.0055169",
"10.1063/5.0052599",
"10.1063/5.0052003",
"10.1364/OL.432817",
"10.1002/ctpp.202000219",
"10.24200/nst.2021.1197",
"10.1088/1367-2630/ac1975",
"10.1103/PhysRevApplied.15.054053",
"10.1103/PhysRevLett.126.234801",
"10.1103/PhysRevE.103.053202",
"10.1088/1361-6587/ac2614",
"10.1088/1361-6587/abf448",
"10.1038/s41467-021-24006-x",
"10.1038/s41567-021-01325-w",
"10.1051/0004-6361/202141049",
"10.1088/1361-6587/ac0352",
"10.1103/PhysRevA.103.053114",
"10.1063/5.0031313",
"10.1039/D0NR07025D",
"10.1016/j.cjph.2021.02.010",
"10.1063/5.0031555",
"10.1038/s41598-021-84264-z",
"10.1063/5.0037028",
"10.1063/5.0040374",
"10.1103/PhysRevLett.126.044801",
"10.1103/PhysRevLett.126.064801",
"10.1103/PhysRevE.103.L021201",


    # '10.1063/5.0037028',
    # ['10.1088/1367-2630/ac1a97', '1907.02621'],
    # ['10.1103/PhysRevLett.126.044801', '1902.05014'],
    # ['10.1088/1742-6596/1596/1/012052', '1912.04064'],
    # ['10.1088/1742-6596/1596/1/012055', '1912.04674'],
    # ['10.1103/PhysRevE.102.033204', '2006.04433'],
    # ['10.1007/s41115-020-0007-6', '2002.09411'],
    # ['10.1029/2019GL086546', '2002.02243'],
    # ['10.1017/S0022377820000264', '1911.09562'],
    # ['10.1103/PhysRevE.101.033204', '1906.05902'],
    # ['10.1017/S0022377819000898', '1905.11131'],
    # ['10.1103/PhysRevE.100.061201', '1911.03440'],
    # ['10.1088/1674-4527/19/12/182', '1908.08170'],
    # '10.3847/2041-8213/ab5b0a',
    # '10.1063/1.5122225',
    # ['10.1088/1361-6587/ab49cf', '1912.04127'],
    # ['10.1016/j.cpc.2019.05.001', '1810.03949'],
    # ['10.1103/PhysRevE.99.033307', '1809.04435'],
    # ['10.1103/PhysRevLett.122.104803', '1806.04976'],
    # ['10.1088/1361-6587/aace22', '1802.02927'],
    # ['10.1093/mnras/sty979', '1712.02883'],
    # ['10.1103/PhysRevE.97.043209', '1707.02618'],
    # ['10.1103/PhysRevE.96.033204', '1705.05402'],
    # '10.1088/1361-6587/aa8a54',
    # '10.1063/1.4996856',
    # '10.1002/2016JA023831',
    # '10.1103/PhysRevE.95.023203',
    # '10.1017/S002237781600057X',
    # '10.1063/1.4955322',
    # '10.1016/j.nima.2016.03.112',
    # '10.1103/PhysRevLett.116.075001',
]

# Some abbreviations
abbrv = {
    "Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment":"Nucl. Inst. Meth. in Phys. Res. A",
}

import unicodedata

def strip_accents(text):
    try:
        text = unicode(text, 'utf-8')
    except NameError: # unicode is a default on python 3 
        pass

    text = unicodedata.normalize('NFD', text)\
           .encode('ascii', 'ignore')\
           .decode("utf-8")

    return str(text)

# Get publication list from DOI numbers
keys = []
citations = []
for doi in DOIs:
    
    if type(doi) is str:
        url = "http://dx.doi.org/"+doi
        arxiv = ""
    elif type(doi) is list:
        url = "http://dx.doi.org/"+doi[0]
        arxiv = doi[1]
    
    try:
        a = loads( get(url, headers={"accept": "application/vnd.citationstyles.csl+json"}).content )
    except Exception as e:
        print("ERROR: could not resolve DOI "+doi)
        raise
    
    year = str(a["issued"]["date-parts"][0][0])
    key = a["author"][0]["family"] + year
    authors = ", ".join([A["given"]+" "+A["family"] for A in a["author"][:-1]]) + " and " + a["author"][-1]["given"]+" "+a["author"][-1]["family"]
    title = a["title"]
    if "container-title-short" in a:
        journal = a["container-title-short"]
    else:
        journal = a["container-title"]
    if journal in abbrv:
        journal = abbrv[journal]
    if "volume" in a:
        volume = a["volume"]
    else:
        volume = "???"
        print("WARNING missing volume "+url)
    if "page" in a:
        page = a["page"]
    elif "article-number" in a:
        page = a["article-number"]
    else:
        page = url.split("/")
        page = page[-1] or page[-2]
    year = year

    s = (
"""
  %s,
  `%s`,
  `%s %s, %s (%s) <%s>`_
"""
    ) % (authors, title, journal, volume, page, year, url)
    
    if arxiv:
        s += "  `arXiv:%s <https://arxiv.org/abs/%s>`_\n" % (arxiv, arxiv)
    
    keys += [key]
    citations += [s]
    
    print( key )

# Handle situations when the label appears several times
count = {k:((keys.count(k)-1) or -1)+1 for k in set(keys)}
abc = " abcdef"

for i in range(len(keys)):
    c = count[keys[i]]
    if c > 0:
        count[keys[i]] -= 1
        keys[i] += abc[c]

# Print results
print("---- copy-paste the following ----")
for k,c in zip(keys, citations):
    print(".. [%s]"% strip_accents(k))
    print(c)
