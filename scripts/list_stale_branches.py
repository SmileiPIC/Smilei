from subprocess import check_output

merged = check_output("git branch -r --merged", shell=True).decode().split()
not_merged = check_output("git branch -r --no-merged", shell=True).decode().split()

for branches in [merged, not_merged]:
    if branches == merged:
        print( "MERGED" )
    else:
        print( "NOT MERGED" )
    for b in branches:
        if b in ["origin/master", "origin/develop"] or "HEAD" in b or b[:7]!="origin/":
            continue
        authors = check_output("git log "+b+" -4 --pretty=format:'%an'", shell=True).decode().split("\n")
        dates = check_output("git log "+b+" -4 --pretty=format:'%ai'", shell=True).decode().split("\n")
        
        years = [int(d[:4]) for d in dates]
        if max(years) == 2022:
            continue
        
        author_date = {}
        for a,d in zip(authors, dates):
            dd = d.split()[0]
            author_date[a] = dd
        author_date = [a+" "+d for a,d in author_date.items()]
        
        print(b + ":\n\t" + "\n\t".join(author_date))
    print("")
