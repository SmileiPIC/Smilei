import re
import os
from tools import *

# DEFINE A CLASS FOR LOGGING DATA
class Log:
    pattern1 = re.compile(""
        +"[\n\t\s]+(Time[ _]in[ _]time[ _]loop)\s+([e+\-.0-9]+)\s+([e+\-<.0-9]+)\% coverage"
        +"([\n\t\s]+([\w ]+)\s+([e+\-.0-9]+)\s+([e+\-<.0-9]+)\%){2,15}"
    )
    pattern2 = re.compile(""
        +"[\t\s]+([\w ]+):?\s+([e+\-.0-9]+)\s+([e+\-<.0-9]+)\%"
    )

    def __init__(self, log_path, log_file):
        mkdir(log_path)
        self.log_file = log_file
        self.data = {}

    def scan(self, output):
        # Open file and find the pattern
        with open(output, 'r') as f:
            text = f.read()
        found = re.search(self.pattern1, text)
        if not found:
            print("WARNING: Unable to log data from "+output)
            return
        lines = found.string[found.start():found.end()].split("\n")[1:]
        matches = [re.search(self.pattern2, l).groups() for l in lines]
        # Get timers values and add to current timers
        for m in matches:
            key = m[0].replace(" ", "")
            if key == "Time_in_time_loop": key = "Timeintimeloop"
            value = float(m[1])
            if key in self.data:
                self.data[key] += value
            else:
                self.data[key] = value

    def append(self):
        if self.data:
            # Append commit and date to current data
            self.data["commit"] = gitversion
            self.data["date"] = strftime("%Y_%m_%d_%H:%M:%S")
            # Open previous database
            try:
                with open(self.log_file, 'r') as f:
                    db = json.load(f)
            except:
                db = {}
            maxlen = max([len(v) for v in db.values()] or [0])
            # Update the database
            for k,v in self.data.items():
                if k in db:
                    db[k] += [None]*(maxlen-len(db[k])) + [v]
                else:
                    db[k] = maxlen*[None] + [v]
            # Overwrite the file
            with open(self.log_file, 'w+') as f:
                json.dump(db, f)
