import re, io, sys

s = io.open('../manuscript.tex', encoding='utf-8').read()

m = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', s, re.S)
ab = m.group(1)
ab = re.sub(r'%.*', '', ab)
ab = re.sub(r'\\[a-zA-Z]+\*?(\[[^\]]*\])?(\{[^{}]*\})?', ' ', ab)
ab = re.sub(r'[$\\{}]', ' ', ab)
words = [w for w in re.split(r'\s+', ab) if re.search(r'[A-Za-z0-9]', w)]
print('abstract words:', len(words))

print('--- headings ---')
for mm in re.finditer(r'\\(section|subsection)\*?\{([^}]*)\}', s):
    print(mm.group(1), '|', mm.group(2)[:70])
print('--- labels ---')
for mm in re.finditer(r'\\label\{([^}]*)\}', s):
    print(mm.group(1))
