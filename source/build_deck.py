"""Build PROTEUS_explained.html from deck_template.html and deck_data.json."""
import json
s = open('deck_template.html').read()
d = json.load(open('deck_data.json'))['bubble']
keep = [i for i, t in enumerate(d['t']) if t <= 4.5]
b = {'t': [round(d['t'][i], 4) for i in keep], 'p': [round(d['p'][i], 3) for i in keep],
     'R': [round(d['R'][i], 5) for i in keep], 'R0': d['R0']}
open('PROTEUS_explained.html', 'w').write(s.replace('__DATA__', json.dumps({'bubble': b})))
print('built')
