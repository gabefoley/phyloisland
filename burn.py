import re

s = "[131249:134744](+)"
m = re.search(r"\[[\w+]+\]", s)
loc = m.group(1)
print(loc)